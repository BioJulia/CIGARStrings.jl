const CIGAR_BYTE_LUT = let
    u = UInt64(0)
    for (i, byte) in enumerate(b"MIDNSHP=X")
        u |= (byte % UInt64) << (7 * i - 7)
    end
    u
end

"""
    BAMCIGAR <: AbstractCIGAR

A BAMCIGAR is an alternative representation of a CIGAR,
stored compactly in 32-bit integers.
Semantically, a BAMCIGAR behaves much similar to a CIGAR.

Construct a BAMCIGAR either from a CIGAR, taking an optional `Vector{UInt8}`
to use as backing storage, or using [`CIGARStrings.try_parse`](@ref).

# Examples
```jldoctest
julia> c = CIGAR("9S123M1=3I15M2H");

julia> b = BAMCIGAR(c, UInt8[]); # use existing storage

julia> c == b
true

julia> CIGAR(b)
CIGAR("9S123M1=3I15M2H")
```

"""
struct BAMCIGAR
    mem::ImmutableMemoryView{UInt8}
    n_ops::UInt32
    aln_len::UInt32
    ref_len::UInt32
    query_len::UInt32

    function BAMCIGAR(
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
    return (CIGARElement(unsafe, u), state + 4)
end

@inline function load_le_u32(mem::ImmutableMemoryView{UInt8}, i::Int)
    return @inbounds begin
        mem[i] % UInt32 |
            (mem[i + 1] % UInt32) << 8 |
            (mem[i + 2] % UInt32) << 16 |
            (mem[i + 3] % UInt32) << 24
    end
end

CIGAR(x::BAMCIGAR) = CIGAR(x, UInt8[])
function CIGAR(x::BAMCIGAR, v::Vector{UInt8})
    empty!(v)
    for element in x
        n = element.len % UInt32
        n_digits = 0
        while !iszero(n)
            (n, r) = divrem(n, UInt32(10))
            push!(v, r + 0x30)
            n_digits += 1
        end
        if n_digits > 1
            @inbounds for i in 1:(n_digits >>> 1)
                a = lastindex(v) - i + 1
                b = lastindex(v) - n_digits + i
                v[a], v[b] = v[b], v[a]
            end
        end
        shift = (7 * (getfield(element, :x) & 0x0f)) & 63
        byte = ((CIGAR_BYTE_LUT >> shift) % UInt8) & 0x7f
        push!(v, byte)
    end
    return CIGAR(
        unsafe,
        ImmutableMemoryView(v),
        x.n_ops,
        x.aln_len,
        x.ref_len,
        x.query_len,
    )
end

BAMCIGAR(x::CIGAR) = BAMCIGAR(x, UInt8[])
function BAMCIGAR(x::CIGAR, v::Vector{UInt8})
    resize!(v, 4 * length(x))
    i = 1
    for element in x
        u = getfield(element, :x)
        for _ in 1:4
            @inbounds v[i] = u % UInt8
            i += 1
            u >>= 8
        end
    end
    return BAMCIGAR(
        unsafe,
        ImmutableMemoryView(v),
        x.n_ops,
        x.aln_len,
        x.ref_len,
        x.query_len
    )
end

function try_parse(::Type{BAMCIGAR}, x)::Union{CIGARError, BAMCIGAR}
    mem = ImmutableMemoryView(x)::ImmutableMemoryView{UInt8}

    # Length must be divisible by 4
    iszero(length(mem) & 0x02) || return CIGARError(1, Errors.NotModFourLength)

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
        if op == OP_H
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
    max(aln_len, n_ops) > typemax(UInt32) && return CIGARError(lastindex(mem) - 3, Errors.IntegerOverflow)
    return BAMCIGAR(unsafe, mem, n_ops % UInt32, aln_len % UInt32, ref_len % UInt32, query_len % UInt32)
end

Base.:(==)(x::BAMCIGAR, y::CIGAR) = y == x
function Base.:(==)(x::CIGAR, y::BAMCIGAR)
    if (x.n_ops != y.n_ops) |
            (x.aln_len != y.aln_len) |
            (x.ref_len != y.ref_len) |
            (x.query_len != y.query_len)
        return false
    end
    return eq_cigars(x, y)
end

@noinline function eq_cigars(x::CIGAR, y::BAMCIGAR)
    for (i, j) in zip(x, y)
        i === j || return false
    end
    return true
end

function BAMCIGAR(x)
    y = try_parse(CIGAR, x)
    return y isa CIGARError ? throw(y) : y
end

query_length(x::BAMCIGAR) = x.query_len % Int
ref_length(x::BAMCIGAR) = x.ref_len % Int
aln_length(x::BAMCIGAR) = x.aln_len % Int
