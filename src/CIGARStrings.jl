module CIGARStrings

#=
TODO:
* Print an AbstractAlignmentReport with seqs
* AlignmentReport
    Show
    findpos
    query,aln,reflength
    aln_identity
=#

export CIGAR,
    OP_M, OP_I, OP_D, OP_N, OP_S, OP_H, OP_P, OP_Eq, OP_X,
    ReportElement, ref_length, aln_length, query_length, aln_identity

using MemViews: ImmutableMemView

struct Unsafe end
const unsafe = Unsafe()

@noinline unreachable() = error("Unreachable reached")

"""
    CIGAROp

1-byte primitive type representing a CIGAR operation.
Module-level constants for each value of this type are defined, e.g. `OP_M`.

| Variable | Char | Description                                          |
|:---------|:-----|:-----------------------------------------------------|
| `OP_M`   | M    | Alignment Match or mismatch                          |
| `OP_I`   | I    | Insertion relative to the reference                  |
| `OP_D`   | D    | Deletion relative to the reference                   |
| `OP_N`   | N    | Region skipped from the reference (usually an intron)|
| `OP_S`   | S    | Soft clip (clipped sequence present in query)        |
| `OP_H`   | H    | Hard clip (clipped sequence not present in query)    |
| `OP_P`   | P    | Padding                                              |
| `OP_Eq`  | Eq   | Alignment match, not mismatch                        |
| `OP_X`   | X    | Alignment mismatch                                   |
"""
primitive type CIGAROp 8 end

const (STRING_LUT, BYTE_LUT) = let
    string_lut = Memory{String}(undef, 9)
    byte_lut = Memory{UInt8}(undef, 9)
    for (i, (name, doc)) in enumerate([
        ("M", "Alignment match or mismatch"),
        ("I", "Insertion relative to the reference"),
        ("D", "Deletion relative to the reference"),
        ("N", "Region skipped from the reference (usually an intron)"),
        ("S", "Soft clip (clipped sequence present in query)"),
        ("H", "Hard clip (clipped sequence not present in query)"),
        ("P", "Padding"),
        ("Eq", "Alignment match, not mismatch"),
        ("X", "Alignment mismatch"),
    ])
        string_lut[i] = name
        sym = Symbol("OP_", name)
        val = reinterpret(CIGAROp, UInt8(i - 1))
        nm = name == "Eq" ? "=" : name
        byte_lut[i] = UInt8(only(nm))
        @eval begin
            @doc $(string("`'", nm, "'`: ", doc)) const $(sym) = reinterpret(CIGAROp, $(UInt8(i - 1)))
        end
    end
    (Tuple(string_lut), Tuple(byte_lut))
end

Base.show(io::IO, op::CIGAROp) = print(io, string("OP_", STRING_LUT[reinterpret(UInt8, op) + 0x01]))

const LUT = let
    n = typemax(UInt128)
    for (char, op) in zip("MIDNSHP=X", [OP_M, OP_I, OP_D, OP_N, OP_S, OP_H, OP_P, OP_Eq, OP_X])
        shift = 4 * (char - '=')
        n &= ~(UInt128(0x0f) << shift)
        n |= UInt128(reinterpret(UInt8, op)) << shift
    end
    n
end

"""
    ReportElement(op::CIGAROp, len::Integer)

Type representing a single element in a CIGAR string, consisting of a
`CIGAROp` and a length.
Access the operation and the length with the properties `.op` and `.len`.
Note that currently, the largest supported length is 268435455.
Operations cannot have length zero.

# Examples
```jldoctest
julia> ReportElement(OP_X, 3)
ReportElement(OP_X, 3)
```

See also: [`CIGAR`](@ref), [`CIGAROp`](@ref)
"""
struct ReportElement
    # 4 upper bits: CIGAROp, 28 lower: Value.
    x::UInt32

    ReportElement(::Unsafe, x::UInt32) = new(x)
    function ReportElement(::Unsafe, op::CIGAROp, len::UInt32)
        new(((reinterpret(UInt8, op) % UInt32) << 28) | len)
    end
end

function ReportElement(op::CIGAROp, len::Integer)
    ulen = UInt32(len)::UInt32
    ulen > 0x0fffffff && error("CIGARStrings only support cigar operations of length 268435455")
    iszero(ulen) && error("CIGAR cannot have zero-length element")
    ReportElement(unsafe, op, ulen)
end

function Base.show(io::IO, x::ReportElement)
    print(io, typeof(x), '(', x.op, ", ", x.len, ')')
end

function Base.getproperty(x::ReportElement, sym::Symbol)
    if sym === :len
        (getfield(x, :x) & 0x0fffffff) % Int
    elseif sym === :op
        reinterpret(CIGAROp, (getfield(x, :x) >>> 28) % UInt8)
    else
        error("No such field in ReportElement: ", sym)
    end
end
Base.propertynames(x::ReportElement) = (:len, :op)

abstract type AbstractAlignmentReport end
Base.eltype(::Type{<:AbstractAlignmentReport}) = ReportElement

function Base.print(out::IO, report::AbstractAlignmentReport)
    for i in iter
        b = @inbounds BYTE_LUT[reinterpret(UInt8, i.op) + 0x01]
        print(out, i.len, b)
    end
end

const CONSUMES = let
    x = UInt32(0)
    for query in [OP_M, OP_I, OP_S, OP_Eq, OP_X]
        shift = reinterpret(UInt8, query) * 2
        x |= UInt32(1) << shift
    end
    for ref in [OP_M, OP_D, OP_N, OP_Eq, OP_X]
        shift = reinterpret(UInt8, ref) * 2 + 1
        x |= UInt32(1) << shift
    end
    x
end

function consumes(x::CIGAROp)::@NamedTuple{query::Bool, ref::Bool}
    n = CONSUMES >>> ((2 * reinterpret(UInt8, x)) & 31)
    (;query=isodd(n), ref=isodd(n >>> 1))
end

@enum TranslationKind::UInt8 Pos Gap Outside

struct Translation
    # 4 upper bits for Kind (in case more is needed)
    x::UInt64

    Translation(::Unsafe, x::UInt64) = new(x)
end

function Translation(kind::TranslationKind, pos::Int)
    pos > 0x0fffffffffffffff && error() # TODO
    Translation(unsafe, kind, pos)
end

function Translation(::Unsafe, kind::TranslationKind, pos::Int)
    n = (pos % UInt64) | ((reinterpret(UInt8, kind) % UInt64) << 60)
    Translation(unsafe, n)
end

function Base.show(io::IO, x::Translation)
    print(io, summary(x), '(', x.kind, ", ", x.pos, ')')
end

Base.propertynames(x::Translation) = (:kind, :pos)
function Base.getproperty(x::Translation, sym::Symbol)
    if sym === :pos
        (getfield(x, :x) & 0x0fffffffffffffff) % Int
    elseif sym === :kind
        reinterpret(TranslationKind, (getfield(x, :x) >>> 60) % UInt8)
    else
        error("No such property in Translation: ", sym)
    end
end

struct Anchor
    query::UInt32
    ref::UInt32
    aln::UInt32
end
Base.zero(::Type{Anchor}) = Anchor(0, 0, 0)

function advance(x::Anchor, e::ReportElement)
    c = consumes(e.op)
    Anchor(
        x.query + (e.len % UInt32) * c.query,
        x.ref   + (e.len % UInt32) * c.ref,
        x.aln   + (e.len % UInt32) * (c.query | c.ref)
    )
end

include("cigar.jl")
include("report.jl")
include("functionality.jl")

end # module CIGARStrings
