module CIGARStrings

export CIGAR,
    OP_M, OP_I, OP_D, OP_N, OP_S, OP_H, OP_P, OP_Eq, OP_X, CIGAROp,
    CIGARElement, ref_length, aln_length, query_length, aln_identity,
    query_to_ref, query_to_aln, ref_to_query, ref_to_aln,
    aln_to_query, aln_to_ref, Translation, count_matches,
    BAMCIGAR, AbstractCIGAR, cigar_view!

public CIGARError, CIGARErrorType, Errors, try_parse, outside, pos, gap, TranslationKind

using MemoryViews: MemoryViews, ImmutableMemoryView, MemoryView

struct Unsafe end
const unsafe = Unsafe()

@noinline unreachable() = error("Unreachable reached")

baremodule Errors

    using Base: @enum

    @enum CIGARErrorType::UInt8 begin
        IntegerOverflow
        InvalidOperation
        ZeroLength
        InvalidHardClip
        InvalidSoftClip
        Truncated
        NotModFourLength
    end

end # baremodule

using .Errors: CIGARErrorType

"""
    CIGARErrorType

Single-byte enum representing the different kind of errors returned or thrown by
the package CIGARStrings.
This currently implemented list of error kinds is not exhaustive, in that more
could be added in minor versions of this package.

See also: [`CIGARError`](@ref)
"""
CIGARErrorType

"""
    CIGARError

Exception kind thrown by the package CIGARStrings.
`CIGARError`s contain two properties: `.kind`, returning a [`CIGARErrorType`](@ref),
and `.index`, returning an `Int`, pointing to the approximate byte index where the
exception was encountered.

The `kind` and `pos` are stable API, in the sense that e.g. an integer overflow
at position 55 will always be represented by a
`CIGARError(55, CIGARStrings.Errors.IntegerOverflow)`.
However, the same invalid CIGAR string may produce multiple different errors, and which
error is produced in that case **is NOT** stable API, because that depends on specifics
parsing internals.

```jldoctest
julia> ce = CIGARStrings.try_parse(CIGAR, "15M9");

julia> ce.index
4

julia> ce.kind
Truncated::CIGARErrorType = 0x05
```
"""
struct CIGARError <: Exception
    index::Int
    kind::CIGARErrorType
end

function Base.showerror(io::IO, error::CIGARError)
    kind = error.kind
    s = if kind == Errors.IntegerOverflow
        "Integer overflow: No operation can be longer than 268435455 nucleotides, and \
        CIGAR alignments can have at most 4294967295 bases."
    elseif kind == Errors.InvalidOperation
        "Invalid operation. Possible values are \"MIDNSHP=X\"."
    elseif kind == Errors.ZeroLength
        "CIGAR element cannot have a length of zero (or the number part of the CIGAR element was missing)"
    elseif kind == Errors.InvalidHardClip
        "Hard clip (H) can only occur at ends of CIGAR."
    elseif kind == Errors.InvalidSoftClip
        "Soft clip (S) must be first or last in CIGAR, or preceded or succeeded by a hard clip (H)."
    elseif kind == Errors.Truncated
        "CIGAR string is truncated and ends too early."
    else
        unreachable()
    end
    return print(io, "Error around byte ", error.index, ": ", s)
end

"""
    CIGAROp

1-byte primitive type representing a CIGAR operation.
Module-level constants for each value of this type are exported, e.g. `OP_M`.

A `CIGAROp` is guaranteed to be a 1-byte primitive with the same bitwise
representation as a `UInt8` with the value given as `N` in the table below,
e.g. `reinterpret(UInt8, OP_I) === 0x01`.

The `Consumes` entry below signifies if the operation advances the position
in the query (`Q`), the reference (`R`) and the alignment (`A`).
E.g. if the `CIGARElement` `9D` begins at query, ref, aln positions `(A, B, C)`,
then the positions are `(A, B+9, C+9)` after. 

| N | Consumes |Variable | Char | Description                                          |
|:--|:---------|:--------|:-----|:-----------------------------------------------------|
| 0 | `Q R A`  |`OP_M`   | M    | Alignment match or mismatch                          |
| 1 | `Q   A`  |`OP_I`   | I    | Insertion relative to the reference                  |
| 2 | `  R A`  |`OP_D`   | D    | Deletion relative to the reference                   |
| 3 | `  R  `  |`OP_N`   | N    | Reference skipped from the alignment (usually an intron)|
| 4 | `Q    `  |`OP_S`   | S    | Soft clip (semantically identical to hard clip)      |
| 5 | `Q    `  |`OP_H`   | H    | Hard clip (semantically identical to soft clip)      |
| 6 | `     `  |`OP_P`   | P    | Padding, position not present in query or reference  |
| 7 | `Q R A`  |`OP_Eq`  | =    | Alignment match, not mismatch                        |
| 8 | `Q R A`  |`OP_X`   | X    | Alignment mismatch                                   |

See also: [`CIGARElement`](@ref)

# Extended help
* `M` means a match or a mismatch. By default, most programs will emit `M`
  instead of `X` or `=`, since the important part of the CIGAR is where the
  insertions and deletions are placed. To determine which bases in an `M`
  are matches and mismatches, the two sequences can be compared, base-wise.
* `N` means that the region is spans an intron, which means the query sequence
  is deleted, but not due to an actual deletion (which would be a `D` operation).
  It can also be used for other uses where the reference bases is missing
  for another reason than a deletion, if such a use case is found.
* `S` and `H` are semantically identical. They both signify the end(s) of the
  query sequence which are not part of the alignment. They differ only in whether
  the clipped bases are written in the SEQ field of a SAM record.
  Typically, hard-clipped bases are present as soft clipped bases in another
  record, such that writing them in both records would be wasteful.
* `P` is only used for a multiple sequence alignment. If a third sequence contains an
  insertion relative to both the query and the reference, the query has a
  `P` at this position, indicating it's present in neither the query nor reference.
* Most alignments only contain the operations `M`, `I`, `D`, `S`.
"""
primitive type CIGAROp 8 end

const (STRING_LUT, BYTE_LUT) = let
    string_lut = Memory{String}(undef, 9)
    byte_lut = Memory{UInt8}(undef, 9)
    for (i, (name, doc)) in enumerate(
            [
                ("M", "Alignment match or mismatch"),
                ("I", "Insertion relative to the reference"),
                ("D", "Deletion relative to the reference"),
                ("N", "Region skipped from the reference (usually an intron)"),
                ("S", "Soft clip (clipped sequence present in query)"),
                ("H", "Hard clip (clipped sequence not present in query)"),
                ("P", "Padding"),
                ("Eq", "Alignment match, not mismatch"),
                ("X", "Alignment mismatch"),
            ]
        )
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
    CIGARElement(op::CIGAROp, len::Integer)

Type representing a single element in a CIGAR string, consisting of a
`CIGAROp` and a length.
For example, in `15S198M1D15M`, the four elements are `15S`, `198M`, `1D`,
and `15M`.
Access the operation and the length with the properties `.op` and `.len`.
Note that currently, the largest supported length is 268435455.
Operations cannot have length zero.

See also: [`CIGAR`](@ref), [`CIGAROp`](@ref)

# Examples
```jldoctest
julia> e = CIGARElement(OP_X, 3)
CIGARElement(OP_X, 3)

julia> e.len
3

julia> e.op
OP_X
```
"""
struct CIGARElement
    # 28 upper bits: Value, lower 4 bits: CIGAROp
    x::UInt32

    CIGARElement(::Unsafe, x::UInt32) = new(x)
    function CIGARElement(::Unsafe, op::CIGAROp, len::UInt32)
        return new((reinterpret(UInt8, op) % UInt32) | (len << 4))
    end
end

function CIGARElement(op::CIGAROp, len::Integer)
    ulen = UInt32(len)::UInt32
    ulen > 0x0fffffff && error("CIGARStrings only support cigar operations of length 268435455")
    iszero(ulen) && error("CIGAR cannot have zero-length element")
    return CIGARElement(unsafe, op, ulen)
end

function Base.show(io::IO, x::CIGARElement)
    return print(io, typeof(x), '(', x.op, ", ", x.len, ')')
end

function Base.getproperty(x::CIGARElement, sym::Symbol)
    return if sym === :len
        (getfield(x, :x) >>> 4) % Int
    elseif sym === :op
        reinterpret(CIGAROp, (getfield(x, :x) % UInt8) & 0x0f)
    else
        error("No such field in CIGARElement: ", sym)
    end
end
Base.propertynames(::CIGARElement) = (:len, :op)

"""
    abstract type AbstractCIGAR

This abstract type is (not yet) a defined interface.
Its concrete subtypes are `CIGAR` and `BAMCIGAR`.
"""
abstract type AbstractCIGAR end

Base.eltype(::Type{<:AbstractCIGAR}) = CIGARElement

include("bytecigar.jl")
include("bamcigar.jl")

const CONSUMES = let
    x = UInt32(0)
    for query in [OP_M, OP_I, OP_S, OP_H, OP_Eq, OP_X] # not PDN
        shift = reinterpret(UInt8, query) * 3
        x |= UInt32(1) << shift
    end
    for ref in [OP_M, OP_D, OP_N, OP_Eq, OP_X] # not PHIS
        shift = reinterpret(UInt8, ref) * 3 + 1
        x |= UInt32(1) << shift
    end
    for aln in [OP_M, OP_I, OP_X, OP_Eq, OP_D] # not PSHN
        shift = reinterpret(UInt8, aln) * 3 + 2
        x |= UInt32(1) << shift
    end
    x
end

"""
    consumes(x::CIGAROp)::@NamedTuple{query::Bool, ref::Bool, aln::Bool}

Return whether the given CIGAROp consumes (i.e. uses up) bases of the query,
reference, and alignment, respectively.
Note that this DOES NOT exactly match the definition as stated in the SAM specs,
because in those, "consumes query" means whether it consumes bases of the query
as printed in the SEQ field, and the SEQ field does not include bases marked as
`OP_H`, despite these actually being part of the query.
"""
function consumes(x::CIGAROp)::@NamedTuple{query::Bool, ref::Bool, aln::Bool}
    n = CONSUMES >>> ((3 * reinterpret(UInt8, x)) & 31)
    query = isodd(n)
    ref = isodd(n >>> 1)
    aln = isodd(n >>> 2)
    return (; query, ref, aln)
end

"""
    TranslationKind

This enum has values `outside`, `pos` and `gap`. It represents the result of
translating a position between query, reference and alignment.
* If `outside`, the input position translates to a position outside the alignment
* If `pos`, the input position corresponds to a non-gap symbol in the alignment
* If `gap`, the input position maps to a gap.

See also: [`Translation`](@ref)
"""
@enum TranslationKind::UInt8 outside pos gap

"""
    Translation(kind::TranslationKind, pos::Int)

The result of translating from a position in the coordinate system of the 
query / reference / alignment to a position in one of the others.
This type contains two documented properies: `kind::TranslationKind` and `pos::Int`.

* If the resulting position is outside the target coordinate, `.kind == outside`
  and `pos == 0`
* If the position maps to a non-gap symbol in the other coodinate system,
  `.kind == pos`, and the position is the `pos`'d symbol in the target coordinate
  system.
* If the position maps to a gap, `pos` is the position of the symbol before the
  gap, and `kind == gap`. When translating to the `aln` coordinate system,
  the kind is never `gap`.

See also: [`CIGARStrings.TranslationKind`](@ref)

# Examples:
```jldoctest
julia> c = CIGAR("2M3D2M2I1M");

julia> for i in 0:8
           refpos = query_to_ref(c, i)
           println(refpos.pos, ' ', refpos.kind)
       end
0 outside
1 pos
2 pos
6 pos
7 pos
7 gap
7 gap
8 pos
0 outside
```
"""
struct Translation
    # 4 upper bits for Kind (in case more is needed)
    x::UInt64

    Translation(::Unsafe, x::UInt64) = new(x)
end

function Translation(kind::TranslationKind, pos::Int)
    unsigned(pos) > 0x0fffffffffffffff && error("Cannot translate position > 1152921504606846975")
    return Translation(unsafe, kind, pos)
end

function Translation(::Unsafe, kind::TranslationKind, pos::Int)
    n = (pos % UInt64) | ((reinterpret(UInt8, kind) % UInt64) << 60)
    return Translation(unsafe, n)
end

const outside_translation = Translation(outside, 0)

function Base.show(io::IO, x::Translation)
    return print(io, summary(x), '(', x.kind, ", ", x.pos, ')')
end

Base.propertynames(::Translation) = (:kind, :pos)
function Base.getproperty(x::Translation, sym::Symbol)
    return if sym === :pos
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

function advance(x::Anchor, e::CIGARElement)
    c = consumes(e.op)
    return Anchor(
        x.query + (e.len % UInt32) * c.query,
        x.ref + (e.len % UInt32) * c.ref,
        x.aln + (e.len % UInt32) * c.aln
    )
end

"""
    query_length(::AbstractCIGAR)::Int

Get the number of biosymbols in the query of the `CIGAR`. This is the same
as the lengths of all `CIGARElement`s of type `M`, `I`, `S`, `H`, `=` and `X`.

See also: [`ref_length`](@ref), [`aln_length`](@ref)

# Example
```jldoctest
julia> query_length(CIGAR("1S5M2D6M2I5M"))
19
```
"""
function query_length end


"""
    ref_length(::AbstractCIGAR)::Int

Get the number of biosymbols in the reference of the `CIGAR`. This is the same
as the lengths of all `CIGARElement`s of type `M`, `D`, `N`, `=` and `X`.

See also: [`query_length`](@ref), [`aln_length`](@ref)

# Example
```jldoctest
julia> ref_length(CIGAR("1S5M2D6M2I5M"))
18
```
"""
function ref_length end

"""
    aln_length(::AbstractCIGAR)::Int

Get the number of biosymbols spanned by the alignment of the `CIGAR`. Clips,
padding and skips are not part of the alignment, but still part of the CIGAR.
Therefore, the alignment length is the same as the lengths of all `CIGARElement`s
of type `M`, `I`, `D`, `=` and `X`.

See also: [`query_length`](@ref), [`ref_length`](@ref)

# Example
```jldoctest
julia> aln_length(CIGAR("1S5M2D6M2I5M"))
20
```
"""
function aln_length end

function count_matches(x::AbstractCIGAR, mismatches::Integer)::Int
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
        error("Mismatches exceed number of possible mismatches in the CIGAR")
    end
    if mismatches < n_X
        error("Mismatches is lower than minimum possible mismatches in CIGAR")
    end
    return (n_Eq % Int) + (n_M - mismatches + n_X) % Int
end

"""
    aln_identity(::AbstractCIGAR, mismatches::Int)::Float64

Compute the alignment identity of the `CIGAR`, computed as the number of matches
(not mismtches) divided by alignment length.
Since the CIGAR itself may not provide information about the precise number of
mismatches, the amount of mismatches is an argument.
The result is always in [0.0, 1.0].

# Example
```jldoctest
julia> aln_identity(CIGAR("3M1D3M1I2M"), 2)
0.6
```
"""
function aln_identity(x::AbstractCIGAR, mismatches::Int)::Float64
    return count_matches(x, mismatches) / aln_length(x)
end

struct PositionMapper{I, C, S, D}
    integers::I
    cigar::C
    get_src::S
    get_dst::D
end

Base.IteratorSize(::Type{<:PositionMapper{I}}) where {I} = Base.IteratorSize(I)
Base.length(x::PositionMapper) = length(x.integers)
Base.eltype(::Type{<:PositionMapper}) = Pair{Int, Translation} # TODO: Only emit translation?
Base.size(x::PositionMapper) = size(x.integers)

function Base.iterate(mapper::PositionMapper)
    (target, next_integer_state) = let
        it = iterate(mapper.integers)
        it === nothing && return nothing
        it
    end
    # Relies on 1 being the first cigar state
    return _iterate(mapper, Int(target)::Int, typemin(Int), next_integer_state, zero(Anchor), 1)
end


function Base.iterate(mapper::PositionMapper, state)
    (integer_state, last_integer, anchor, cigar_state) = state
    # If no more integers, done
    (target, next_integer_state) = let
        it = iterate(mapper.integers, integer_state)
        it === nothing && return nothing
        it
    end
    return _iterate(mapper, Int(target)::Int, last_integer, next_integer_state, anchor, cigar_state)
end

function _iterate(mapper::PositionMapper, target::Int, last_integer::Int, next_integer_state, anchor, cigar_state)
    # Integers must be ascending
    target ≥ last_integer || error("Integers must be ascending in order")
    new_state = (next_integer_state, target, anchor, cigar_state)
    while true
        # If no more cigar elements, emit outside_translation
        (element, next_cigar_state) = let
            it = iterate(mapper.cigar, cigar_state) # TODO: Should the element be part of the state? To avoid getting the same element many times
            it === nothing && return (target => outside_translation, new_state)
            it
        end
        in(element.op, (OP_S, OP_H)) && return (target => outside_translation, new_state)
        next_anchor = advance(anchor, element)

        # The first time we overshoot the target, we hit the right operation.
        if mapper.get_src(next_anchor) ≥ target
            @show anchor
            # Non-gap elements increment the source position.
            if mapper.get_dst(next_anchor) > mapper.get_dst(anchor)
                return (target => Translation(unsafe, pos, mapper.get_dst(anchor) + (target - mapper.get_src(anchor))), new_state)
            else
                return (target => Translation(unsafe, gap, Int(mapper.get_dst(anchor))), new_state)
            end
        end
        anchor = next_anchor
        cigar_state = next_cigar_state
    end
    return
end


# CIGARs don't support random indexing, so we need to do a linear search.
# In this function and the following, if you search for query_to_ref, then
# get_src is a function f(x::Anchor) = x.query, and get_dst f(x::Anchor) = x.ref.
# The comments will, for e.g. query_to_ref, refer to query as src and ref as dst.
function pos_to_pos(
        aln::AbstractCIGAR,
        target::Int,
        get_src::Function,
        get_dst::Function,
    )::Translation
    return target < 1 ? outside_translation : pos_to_pos_linear(aln, target, get_src, get_dst, zero(Anchor))
end

@inline function pos_to_pos_linear(
        CIGARElement_iter,
        target::Int,
        get_src::Function,
        get_dst::Function,
        prev_anchor::Anchor
    )::Translation
    kind = outside
    for element::CIGARElement in CIGARElement_iter
        anchor = advance(prev_anchor, element)
        if get_src(anchor) ≥ target
            in(element.op, (OP_S, OP_H)) && return outside_translation
            # We know the src must have been incremented to at or above target by this element.
            # If the dst was also incremented, we might have overshot.
            if get_dst(anchor) > get_dst(prev_anchor)
                kind = pos
                break
                # If on the other hand the ref was not incremented, then we have hit
                # an element that increments our src (query) but not dst (ref),
                # and so we just report a gap from our previous anchor's dst.
            else
                kind = gap
                break
            end
        end
        prev_anchor = anchor
    end
    return if kind == outside
        outside_translation
    elseif kind == pos
        Translation(unsafe, pos, get_dst(prev_anchor) + (target - get_src(prev_anchor)))
    else
        Translation(unsafe, gap, Int(get_dst(prev_anchor)))
    end
end

"""
    query_to_ref(x::AbstractCIGAR, pos::Int)::Int

Get the 1-based reference position aligning to query position `pos`.
See [`Translation`](@ref) for more details.

# Examples
```jldoctest
julia> c = CIGAR("3M2D1M4I2M");

julia> query_to_ref(c, 4)
Translation(pos, 6)
```
"""
function query_to_ref(x::AbstractCIGAR, pos::Int)
    return @inline pos_to_pos(x, pos, i -> i.query, i -> i.ref)
end

"""
    query_to_aln(x::AbstractCIGAR, pos::Int)::Int

Get the 1-based alignment position aligning to query position `pos`.
See [`Translation`](@ref) for more details.

# Examples
```jldoctest
julia> c = CIGAR("3M2D1M4I2M");

julia> query_to_aln(c, 8)
Translation(pos, 10)
```
"""
function query_to_aln(x::AbstractCIGAR, pos::Int)
    return @inline pos_to_pos(x, pos, i -> i.query, i -> i.aln)
end

"""
    ref_to_query(x::AbstractCIGAR, pos::Int)::Int

Get the 1-based query position aligning to reference position `pos`.
See [`Translation`](@ref) for more details.

# Examples
```jldoctest
julia> c = CIGAR("3M2D1M4I2M");

julia> ref_to_query(c, 7)
Translation(pos, 9)
```
"""
function ref_to_query(x::AbstractCIGAR, pos::Int)
    return @inline pos_to_pos(x, pos, i -> i.ref, i -> i.query)
end

"""
    ref_to_aln(x::AbstractCIGAR, pos::Int)::Int

Get the 1-based alignment position aligning to reference position `pos`.
See [`Translation`](@ref) for more details.

# Examples
```jldoctest
julia> c = CIGAR("3M2D1M4I2M");

julia> ref_to_aln(c, 7)
Translation(pos, 11)
```
"""
function ref_to_aln(x::AbstractCIGAR, pos::Int)
    return @inline pos_to_pos(x, pos, i -> i.ref, i -> i.aln)
end

"""
    aln_to_query(x::AbstractCIGAR, pos::Int)::Int

Get the 1-based query position aligning to alignment position `pos`.
See [`Translation`](@ref) for more details.

# Examples
```jldoctest
julia> c = CIGAR("3M2D1M4I2M");

julia> aln_to_query(c, 10)
Translation(pos, 8)
```
"""
function aln_to_query(x::AbstractCIGAR, pos::Int)
    return @inline pos_to_pos(x, pos, i -> i.aln, i -> i.query)
end

"""
    aln_to_ref(x::AbstractCIGAR, pos::Int)::Int

Get the 1-based reference position aligning to alignment position `pos`.
See [`Translation`](@ref) for more details.

# Examples
```jldoctest
julia> c = CIGAR("3M2D1M4I2M");

julia> aln_to_ref(c, 9)
Translation(gap, 6)
```
"""
function aln_to_ref(x::AbstractCIGAR, pos::Int)
    return @inline pos_to_pos(x, pos, i -> i.aln, i -> i.ref)
end

end # module CIGARStrings
