"""
    CIGAR <: AbstractCIGAR

A CIGAR string represents the sequence of insertions, matches and deletions
that comprise a pairwise alignment.
Construct a `CIGAR` from any object `x` where `MemoryView(x)` returns a
`MemoryView{UInt8}`, i.e. any memory-backed bytearray, or string.

Use [`CIGARStrings.try_parse`](@ref) to attempt to parse a `CIGAR` string
without throwing an exception if the data is invalid.

See also: [`CIGARElement`](@ref)

# Extended help
CIGAR strings are sequences of `CIGARElement`, from the 5' to the 3' of the
query (or N- to C-terminal for amino acids).
CIGAR strings comprise the entire query, i.e. the sum of lengths of elements
with the `XMI=SH` operations equals the length of the query.

For example, the query `AGCGTAGCACACC` that aligns from query base 5 and ref
base 1002, like this:

    Q: 5    TAG--CACACC   13
    R: 1002 TAGGACAC-CC 1011

Is summarized by the CIGAR `4S3M2D3M1I2M`. The operations `HX=PN` are more rarely
used, see [`CIGAROp`](@ref) for a description of the operations.
"""
struct CIGAR <: AbstractCIGAR
    mem::ImmutableMemoryView{UInt8}
    n_ops::UInt32
    aln_len::UInt32
    ref_len::UInt32
    query_len::UInt32

    function CIGAR(
            ::Unsafe,
            mem::ImmutableMemoryView{UInt8},
            n_ops::UInt32,
            aln_len::UInt32,
            ref_len::UInt32,
            query_len::UInt32,
        )
        return new(mem, n_ops, aln_len, ref_len, query_len)
    end
end

function CIGAR(x)
    y = try_parse(CIGAR, x)
    return y isa CIGARError ? throw(y) : y
end

MemoryViews.MemoryView(x::CIGAR) = x.mem

function Base.iterate(
        x::CIGAR,
        state::Int = 1
    )::Union{Tuple{CIGARElement, Int}, Nothing}
    mem = x.mem
    len = length(mem)
    state > len && return nothing
    n = 0
    @inbounds for i in state:len
        b = mem[i] - 0x30
        if b < 0x0a
            n = (10 * n) + b
        else
            b -= UInt8('=') - UInt8('0')
            enc = @inbounds OP_LUT[min(b + 1, 28)]
            op = reinterpret(CIGAROp, enc)
            return (CIGARElement(unsafe, op, n % UInt32), i + 1)
        end
    end
    return unreachable()
end

# TODO: Limit width for long cigars
function Base.show(io::IO, x::CIGAR)
    buf = IOBuffer()
    print(buf, summary(x), "(\"")
    write(buf, MemoryView(x))
    print(buf, "\")")
    return write(io, take!(buf))
end

Base.length(x::CIGAR) = x.n_ops % Int

"""
    try_parse(::Type{CIGAR}, x)::Union{CIGAR, CIGARError}

Cast `x` to a `MemoryView{UInt8}`, and try parsing a [`CIGAR`](@ref) from it.
If the parsing is unsuccessful, return a [`CIGARError`](@ref)

# Examples
```jldoctest
julia> c = CIGARStrings.try_parse(CIGAR, "2S1M9I");

julia> c isa CIGAR # success
true

julia> c = CIGARStrings.try_parse(CIGAR, "1S7H9M1S");

julia> c.kind
InvalidHardClip::CIGARErrorType = 0x03
```
"""
function try_parse(::Type{CIGAR}, x)::Union{CIGARError, CIGAR}
    mem = ImmutableMemoryView(x)::ImmutableMemoryView{UInt8}
    # H must be either first or last
    # S must be either first or last, or preceded or succeeded by H
    n_ops = 0
    aln_len = 0
    ref_len = 0
    query_len = 0
    last_was_num = false
    is_first = true
    last_was_H = false
    next_must_be_H = false
    n = 0
    for i in eachindex(mem)
        b = @inbounds mem[i]
        b -= 0x30
        last_was_num = b < 0x0a
        if last_was_num
            n = (10 * n) + b
            n > 0x0fffffff && return CIGARError(i, Errors.IntegerOverflow)
        else
            b -= UInt8('=') - UInt8('0')
            b > (UInt8('X') - UInt8('=')) && return CIGARError(i, Errors.InvalidOperation)
            enc = @inbounds OP_LUT[b + 1]
            enc == 0xff && return CIGARError(i, Errors.InvalidOperation)
            iszero(n) && return CIGARError(i, Errors.ZeroLength)
            op = reinterpret(CIGAROp, enc)
            if op === OP_H
                if !is_first && i != lastindex(mem)
                    return CIGARError(i, Errors.InvalidHardClip)
                end
            elseif next_must_be_H
                return CIGARError(i, Errors.InvalidSoftClip)
            elseif op === OP_S
                if !is_first && !last_was_H
                    next_must_be_H = true
                end
            end
            c = consumes(op)
            n_ops += 1
            ref_len += c.ref * n
            query_len += c.query * n
            aln_len += c.aln * n
            n = 0
            last_was_H = op === OP_H
            is_first = false
        end
    end
    last_was_num && return CIGARError(lastindex(mem), Errors.Truncated)
    if max(aln_len, n_ops, ref_len, query_len) > typemax(UInt32)
        return CIGARError(lastindex(mem), Errors.IntegerOverflow)
    end
    return CIGAR(unsafe, mem, n_ops % UInt32, aln_len % UInt32, ref_len % UInt32, query_len % UInt32)
end

function Base.print(out::IO, cigar::CIGAR)
    return write(out, cigar.mem)
end

function Base.:(==)(x::CIGAR, y::CIGAR)
    return length(x) == length(y) &&
        x.aln_len == y.aln_len &&
        x.ref_len == y.ref_len &&
        x.query_len == y.query_len &&
        x.mem == y.mem
end

query_length(x::CIGAR) = x.query_len % Int
ref_length(x::CIGAR) = x.ref_len % Int
aln_length(x::CIGAR) = x.aln_len % Int

"""
    unsafe_switch_memory(cigar::T, mem::ImmutableMemoryView{UInt8})::T where {T <: AbstractCIGAR}

Create a new instance of `typeof(cigar)` equal to `cigar`, but using the new memory
`mem` which must be equal to the existing memory backing `cigar`.
This operation does not do any validation.

This function is unsafe, because it assumes that `mem == MemoryView(cigar)`.
If this assumption is violated, any subsequent operation on the resulting `AbstractCIGAR`
may cause undefined behaviour.

# Examples
```jldoctest
julia> mem = MemoryView("5S12M1X8M10S");

julia> cigar_1 = CIGAR(mem);

julia> cigar_2 = unsafe_switch_memory(cigar_1, copy(mem));

julia> cigar_1 == cigar_2
true

julia> MemoryView(cigar_1) === MemoryView(cigar_2)
false
```
"""
function unsafe_switch_memory(x::CIGAR, mem::ImmutableMemoryView{UInt8})
    return CIGAR(unsafe, mem, x.n_ops, x.aln_len, x.ref_len, x.query_len)
end

function count_matches(x::CIGAR, mismatches::Integer)::Int
    mismatches = UInt(mismatches)::UInt
    n_M = UInt(0)
    n_X = UInt(0)
    n_Eq = UInt(0)
    for i in x
        n_Eq += (i.op === OP_Eq) * (i.len % UInt)
        n_X += (i.op === OP_X) * (i.len % UInt)
        n_M += (i.op === OP_M) * (i.len % UInt)
    end
    if mismatches > n_M + n_X
        throw(DomainError(mismatches, "Mismatches exceed number of possible mismatches in the CIGAR"))
    end
    if mismatches < n_X
        throw(DomainError(mismatches, "Mismatches is lower than minimum possible mismatches in CIGAR"))
    end
    return (n_Eq % Int) + (n_M - mismatches + n_X) % Int
end

# N.B: This function makes use of the iterator state of `cigar`
function normalize!(cigar::CIGAR, mem::MutableMemoryView{UInt8})::CIGAR
    state = 1
    it = iterate(cigar, state)
    # dummy value for type stability, initial value is not used
    last_element = CIGARElement(unsafe, UInt32(0x00_00_00_10))

    # When merging elements, we overwrite them. So, to overwrite at the right
    # location, we begin writing at this location. Dummy value to begin with.
    previous_element_write_offset = 0
    n_ops = 0
    write_offset = 0
    source_mem = MemoryView(cigar)
    while it !== nothing
        (element, new_state) = it
        if element.op === OP_P
            state = new_state
            it = iterate(cigar, state)
            continue
        end
        op = if element.op === OP_H
            OP_S
        elseif element.op ∈ (OP_Eq, OP_X)
            OP_M
        else
            element.op
        end
        shift = (reinterpret(UInt8, op) * 0x07) & 63
        op_byte = ((CIGAR_BYTE_LUT >> shift) % UInt8) & 0x7f
        if last_element.op == op && write_offset > 0
            # Same operation, overwrite the previous one.
            # Not applicable for the first operation
            len = element.len + last_element.len
            if len > 268435455
                throw(CIGARError(state, Errors.IntegerOverflow))
            end
            len = len % UInt32

            # Update last element
            u = (len << 4) | reinterpret(UInt8, op)
            last_element = CIGARElement(unsafe, u)

            # Write digits, backwards, e.g. write "123" as "321".
            # We write them backwards because that's easiest when using
            # the divrem approach
            n_digits = 0
            write_offset = previous_element_write_offset
            previous_element_write_offset = write_offset
            while len > 0
                n_digits += 1
                (len, rm) = divrem(len, UInt32(10))
                @boundscheck checkbounds(mem, write_offset + 1)
                @inbounds mem[(write_offset += 1)] = rm % UInt8 + 0x30
            end
            # Now, reverse digits in-place
            @inbounds for j in 1:(n_digits >> 1)
                a = write_offset - j + 1
                b = write_offset - n_digits + j
                (mem[a], mem[b]) = (mem[b], mem[a])
            end
            # Now write the operation itself
            @boundscheck checkbounds(mem, write_offset + 1)
            @inbounds mem[(write_offset += 1)] = op_byte
        else
            previous_element_write_offset = write_offset
            # New operation, so we simply copy the memory of the source to `mem`
            n_bytes = new_state - state
            @boundscheck if length(mem) < write_offset + n_bytes
                throw(BoundsError(mem, write_offset + n_bytes))
            end
            # Copy over digits
            @inbounds for i in state:(new_state - 2)
                mem[(write_offset += 1)] = source_mem[i]
            end
            # Copy over new op
            @inbounds mem[(write_offset += 1)] = op_byte

            n_ops += 1

            # Update last element
            u = (getfield(element, :x) & 0xff_ff_ff_f0) | reinterpret(UInt8, op)
            last_element = CIGARElement(unsafe, u)
        end
        state = new_state
        it = iterate(cigar, state)
    end
    mem = @inbounds ImmutableMemoryView(mem)[1:write_offset]
    return CIGAR(unsafe, mem, n_ops % UInt32, cigar.aln_len, cigar.ref_len, cigar.query_len)
end

function normalize(cigar::CIGAR)
    return @inbounds @inline normalize!(cigar, similar(MemoryView(cigar)))
end

function unsafe_normalize(cigar::CIGAR)
    mem = MemoryView(cigar)
    newmem = MemoryViews.unsafe_from_parts(mem.ref, length(mem))
    # This works even though the input and output arrays alias, because
    # they alias with the same offset. Since the normalization process happens
    # in a single pass over the input array, and it is never the case that
    # normalization writes to higher indices than it reads from, we never read
    # any data that has been modified in the normalization process.
    return @inbounds @inline normalize!(cigar, newmem)
end
