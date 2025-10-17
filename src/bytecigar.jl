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
            enc = ((LUT >>> ((4 * b) & 127)) % UInt8) & 0x0f
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
            enc = ((LUT >>> ((4 * b) & 127)) % UInt8) & 0x0f
            enc == 0x0f && return CIGARError(i, Errors.InvalidOperation)
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
    max(aln_len, n_ops) > typemax(UInt32) && return CIGARError(lastindex(mem), Errors.IntegerOverflow)
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
If this assumption is violated, any subsequent operation on the resuling `AbstractCIGAR`
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
