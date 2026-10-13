const CIGAR_BYTE_LUT = let
    u = UInt64(0)
    for (i, byte) in enumerate(b"MIDNSHP=X")
        u |= (byte % UInt64) << (7 * i - 7)
    end
    u
end

"""
    BAMCIGAR <: AbstractCIGAR

A `BAMCIGAR` is an alternative representation of a `CIGAR`,
stored compactly in 32-bit integers.
Semantically, a `BAMCIGAR` behaves much like a `CIGAR`.

Construct a `BAMCIGAR` either from a `CIGAR`, or using [`encode!`](@ref)
or [`encode_append!`](@ref)

# Examples
```jldoctest
julia> c = CIGAR("9S123M1=3I15M2H");

julia> b = BAMCIGAR(c);

julia> c == b
true

julia> CIGAR(b)
CIGAR("9S123M1=3I15M2H")
```
"""
struct BAMCIGAR <: AbstractCIGAR
    mem::ImmutableMemoryView{UInt8}
    aln_len::UInt32
    ref_len::UInt32
    query_len::UInt32

    global function unsafe_bamcigar(
            mem::ImmutableMemoryView{UInt8},
            aln_len::UInt32,
            ref_len::UInt32,
            query_len::UInt32,
        )
        return new(mem, aln_len, ref_len, query_len)
    end
end

function Base.show(io::IO, x::BAMCIGAR)
    mem = MemoryView(CIGAR(x))
    print(io, "BAMCIGAR(CIGAR(\"")
    write(io, mem)
    return print(io, "\"))")
end

Base.length(x::BAMCIGAR) = length(x.mem) >>> 2

function Base.iterate(
        x::BAMCIGAR,
        state::Int = 1
    )::Union{Tuple{CIGARElement, Int}, Nothing}
    (state - 1) % UInt >= (length(x.mem) % UInt) && return nothing
    u = load_le_u32(x.mem, state)
    return (unsafe_cigarelement(u), state + 4)
end

@inline function load_le_u32(mem::ImmutableMemoryView{UInt8}, i::Int)
    return @inbounds begin
        mem[i] % UInt32 |
            (mem[i + 1] % UInt32) << 8 |
            (mem[i + 2] % UInt32) << 16 |
            (mem[i + 3] % UInt32) << 24
    end
end

CIGAR(x::BAMCIGAR) = encode_append!(UInt8[], CIGAR, x)


MemoryViews.MemoryView(x::BAMCIGAR) = x.mem

BAMCIGAR(x::CIGAR) = encode_append!(UInt8[], BAMCIGAR, x)

function Base.print(io::IO, x::BAMCIGAR)
    v = UInt8[]
    encode_append!(v, CIGAR, x)
    write(io, v)
    return nothing
end

function try_parse(::Type{BAMCIGAR}, x::AbstractArray{UInt8})::Union{CIGARError, BAMCIGAR}
    mem = ImmutableMemoryView(x)::ImmutableMemoryView{UInt8}

    # Length must be divisible by 4
    iszero(length(mem) & 0x03) || return CIGARError(1, Errors.NotModFourLength)

    aln_len = 0
    ref_len = 0
    query_len = 0
    is_first = true
    last_was_H = false
    next_must_be_H = false

    n_ops = length(mem) >>> 2

    for u_offset in 0:(n_ops - 1)
        idx = 4 * u_offset + 1
        u = load_le_u32(mem, idx)
        n = u >> 4
        iszero(n) && return CIGARError(idx, Errors.ZeroLength)
        enc = (u % UInt8) & 0x0f
        enc > 0x08 && return CIGARError(idx, Errors.InvalidOperation)
        op = reinterpret(CIGAROp, enc)
        if op === OP_H
            if !is_first && u_offset < n_ops - 1
                return CIGARError(idx, Errors.InvalidHardClip)
            end
        elseif next_must_be_H
            return CIGARError(idx, Errors.InvalidSoftClip)
        elseif op === OP_S
            if !is_first && !last_was_H
                next_must_be_H = true
            end
        end

        c = consumes(op)
        ref_len += c.ref * n
        query_len += c.query * n
        aln_len += c.aln * n
        last_was_H = op === OP_H
        is_first = false
    end
    if max(aln_len, ref_len, query_len) > typemax(UInt32)
        return CIGARError(lastindex(mem), Errors.IntegerOverflow)
    end
    return unsafe_bamcigar(mem, aln_len % UInt32, ref_len % UInt32, query_len % UInt32)
end

function Base.:(==)(x::BAMCIGAR, y::BAMCIGAR)
    return length(x) == length(y) &&
        x.aln_len == y.aln_len &&
        x.ref_len == y.ref_len &&
        x.query_len == y.query_len &&
        x.mem == y.mem
end

Base.:(==)(x::BAMCIGAR, y::CIGAR) = y == x
function Base.:(==)(x::CIGAR, y::BAMCIGAR)
    if (length(x) != length(y)) |
            (x.aln_len != y.aln_len) |
            (x.ref_len != y.ref_len) |
            (x.query_len != y.query_len)
        return false
    end
    return eq_cigars(x, y)
end

function BAMCIGAR(x::AbstractArray{UInt8})
    y = try_parse(BAMCIGAR, x)
    return y isa CIGARError ? throw(y) : y
end

query_length(x::BAMCIGAR) = x.query_len % Int
ref_length(x::BAMCIGAR) = x.ref_len % Int
aln_length(x::BAMCIGAR) = x.aln_len % Int

function unsafe_switch_memory(x::BAMCIGAR, mem::ImmutableMemoryView{UInt8})
    return unsafe_bamcigar(mem, x.aln_len, x.ref_len, x.query_len)
end

function count_matches(x::BAMCIGAR, mismatches::Integer)::Union{Nothing, Int}
    n_M = UInt(0)
    n_X = UInt(0)
    n_Eq = UInt(0)
    mem = x.mem
    GC.@preserve mem begin
        @inbounds for i in 1:4:length(mem)
            u = ltoh(unsafe_load(Ptr{UInt32}(pointer(mem, i))))
            op = u & 0x0f
            len = (u >>> 4) % UInt
            n_Eq += (op == 0x07) * len
            n_X += (op == 0x08) * len
            n_M += (op == 0x00) * len
        end
    end
    return count_matches_from_counts(n_M % Int, n_X % Int, n_Eq % Int, mismatches)
end

function normalize!(cigar::BAMCIGAR, mem::MutableMemoryView{UInt8})::BAMCIGAR
    write_offset = 0
    # dummy value for type stability, initial value is not used
    last_element = unsafe_cigarelement(UInt32(0x00_00_00_10))
    GC.@preserve mem begin
        for (element_number, element) in enumerate(cigar)
            element.op ∈ (OP_P, OP_H) && continue
            op = if element.op ∈ (OP_X, OP_Eq)
                OP_M
            else
                element.op
            end
            u = if op === last_element.op && write_offset > 0
                # Overwrite the previous one
                len = element.len + last_element.len
                if len > 268435455
                    # Index of the first byte of `element` in the input
                    throw(CIGARError(4 * element_number - 3, Errors.IntegerOverflow))
                end
                write_offset -= 4
                ((len % UInt32) << 4) | reinterpret(UInt8, op)
            else
                # Append to end
                @boundscheck checkbounds_lightboundserror(mem, write_offset + 4)
                (getfield(element, :x) & 0xff_ff_ff_f0) | reinterpret(UInt8, op)
            end
            ptr = Ptr{UInt32}(pointer(mem) + write_offset)
            unsafe_store!(ptr, htol(u))
            write_offset += 4
            last_element = unsafe_cigarelement(u)
        end
    end
    mem = @inbounds ImmutableMemoryView(mem)[1:write_offset]
    return unsafe_bamcigar(mem, cigar.aln_len, cigar.ref_len, cigar.query_len)
end

function normalize(cigar::BAMCIGAR)
    return @inbounds @inline normalize!(cigar, similar(MemoryView(cigar)))
end

function unsafe_normalize(cigar::BAMCIGAR)
    mem = MemoryView(cigar)
    newmem = MemoryViews.unsafe_from_parts(MemoryViews.unsafe_memoryref(mem), length(mem))
    # See comment in `unsafe_normalize(::CIGAR)` for why aliasing is valid
    return @inbounds @inline normalize!(cigar, newmem)
end
