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
to use as backing storage, or using [`CIGARStrings.try_parse`](@ref),
or [`BAMCIGAR(::MutableMemoryView{UInt8}, ::CIGAR)`](@ref)

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
struct BAMCIGAR <: AbstractCIGAR
    mem::ImmutableMemoryView{UInt8}
    aln_len::UInt32
    ref_len::UInt32
    query_len::UInt32

    function BAMCIGAR(
            ::Unsafe,
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
    return CIGAR(
        unsafe,
        cigar_view!(v, x),
        length(x) % UInt32,
        x.aln_len,
        x.ref_len,
        x.query_len,
    )
end

MemoryViews.MemoryView(x::BAMCIGAR) = x.mem

"""
    cigar_view!(v::Vector{UInt8}, x::BAMCIGAR)::ImmutableMemoryView{UInt8}

Write the ASCII (i.e. `CIGAR`) representation `x` into `v`,
emptying `v`'s original content.
A memory view of `v` is returned:

# Examples
```jldoctest
julia> v = [0x01, 0x02, 0x03];

julia> bc = BAMCIGAR(CIGAR("151M3D20M"));

julia> mem_view = cigar_view!(v, bc);

julia> mem_view == v
true

julia> String(mem_view) == string(CIGAR(bc))
true
```
"""
function cigar_view!(v::Vector{UInt8}, x::BAMCIGAR)
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
    return ImmutableMemoryView(v)
end

BAMCIGAR(x::CIGAR) = BAMCIGAR(x, UInt8[])

function BAMCIGAR(x::CIGAR, v::Vector{UInt8})
    resize!(v, 4 * length(x))
    return @inbounds BAMCIGAR(MemoryView(v), x)
end

"""
    BAMCIGAR(mem::MutableMemoryView{UInt8}, x::CIGAR)::BAMCIGAR

Construct a `BAMCIGAR` equal to `x`, using the memory `mem`.
After calling this, `mem` may not be mutated, and is considered
owned by the resulting `BAMCIGAR`.

Throw a `BoundsError` if `length(mem) < 4 * length(x)`.

# Examples
```jldoctest
julia> x = CIGAR("150M3D9S");

julia> mem = MemoryView(zeros(UInt8, 15));

julia> cigar = BAMCIGAR(mem, x)
BAMCIGAR(CIGAR("150M3D9S"))

julia> parent(MemoryView(cigar)) === parent(mem)
true
```
"""
function BAMCIGAR(mem::MutableMemoryView{UInt8}, x::CIGAR)
    @boundscheck if length(mem) < 4 * length(x)
        throw(BoundsError(mem, 4 * length(x)))
    end
    mem = @inbounds mem[1:(4 * length(x))]
    i = 1
    for element in x
        u = getfield(element, :x)
        for _ in 1:4
            @inbounds mem[i] = u % UInt8
            i += 1
            u >>= 8
        end
    end
    return BAMCIGAR(
        unsafe,
        ImmutableMemoryView(mem),
        x.aln_len,
        x.ref_len,
        x.query_len
    )
end

Base.print(io::IO, x::BAMCIGAR) = (write(io, cigar_view!(UInt8[], x)); nothing)

function try_parse(::Type{BAMCIGAR}, x)::Union{CIGARError, BAMCIGAR}
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
    return BAMCIGAR(unsafe, mem, aln_len % UInt32, ref_len % UInt32, query_len % UInt32)
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
    if length(x) != length(y) |
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
    y = try_parse(BAMCIGAR, x)
    return y isa CIGARError ? throw(y) : y
end

query_length(x::BAMCIGAR) = x.query_len % Int
ref_length(x::BAMCIGAR) = x.ref_len % Int
aln_length(x::BAMCIGAR) = x.aln_len % Int

function unsafe_switch_memory(x::BAMCIGAR, mem::ImmutableMemoryView{UInt8})
    return BAMCIGAR(unsafe, mem, x.aln_len, x.ref_len, x.query_len)
end
