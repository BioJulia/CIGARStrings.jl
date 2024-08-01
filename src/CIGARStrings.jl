module CIGARStrings

export CIGAR,
    OP_M, OP_I, OP_D, OP_N, OP_S, OP_H, OP_P, OP_Eq, OP_X, CIGAROp,
    CIGARElement, ref_length, aln_length, query_length, aln_identity,
    query_to_ref, query_to_aln, ref_to_query, ref_to_aln,
    aln_to_query, aln_to_ref

public CIGARError, CIGARErrorType, Errors

using MemoryViews: ImmutableMemoryView, MemoryView

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
end

end # baremodule

using .Errors: CIGARErrorType

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
    print(io, "Error around byte ", error.index, ": ", s)
end

"""
    CIGAROp

1-byte primitive type representing a CIGAR operation.
Module-level constants for each value of this type are defined, e.g. `OP_M`.

A `CIGAROp` is guaranteed to be a 1-byte primitive with the same bitwise
representation as a `UInt8` with the value given as `N` in the table below,
e.g. `reinterpret(UInt8, OP_I) === 0x01`.

The `Consumes` entry below signifies if the operation advances the position
in the query (`Q`), the reference (`R`) and the alignment (`A`).
E.g. if the `CIGARElement` `9D` begins at query, ref, aln positions `(0, 0, 0)`,
then the positions are `(0, 9, 9)` after. 

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

# Extended help
* `M` means a match or a mismatch. By default, most programs will emit `M`
  instead of `X` or `=`, since the important part of the CIGAR is where the
  insertions and deletions are placed. To determine which bases in an `M`
  are matches and mismatches, the two sequences need to be compared.
* `N` means that the region is spans an intron, which means the sequence
  is deleted, but not due to an actual deletion.
  It can also be used for other uses where the reference bases is missing
  for another reason than a deletion, if such a use case is found.
* `S` and `H` are semantically identical. They both signify the end(s) of the
  sequence which are not part of the alignment. They differ only in whether
  the covered bases are written in the SEQ field of a SAM record.
  Typically, hard-clipped bases are present as soft clipped bases in another
  record, such that writing them in both records would be wasteful.
* `P` is only used for multiple alignment. If a third sequence contains an
  insertion relative to both the query and the reference, the query has a
  `P` at this position, indicating it's present in neither the query or ref.
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
    CIGARElement(op::CIGAROp, len::Integer)

Type representing a single element in a CIGAR string, consisting of a
`CIGAROp` and a length.
For example, in `15S198M1D15M`, the four elements are `15S`, `198M`, `1D`,
and `15M`.
Access the operation and the length with the properties `.op` and `.len`.
Note that currently, the largest supported length is 268435455.
Operations cannot have length zero.

# Examples
```jldoctest
julia> e = CIGARElement(OP_X, 3)
CIGARElement(OP_X, 3)

julia> e.len
3

julia> e.op
OP_X
```

See also: [`CIGAR`](@ref), [`CIGAROp`](@ref)
"""
struct CIGARElement
    # 4 upper bits: CIGAROp, 28 lower: Value.
    x::UInt32

    CIGARElement(::Unsafe, x::UInt32) = new(x)
    function CIGARElement(::Unsafe, op::CIGAROp, len::UInt32)
        new(((reinterpret(UInt8, op) % UInt32) << 28) | len)
    end
end

function CIGARElement(op::CIGAROp, len::Integer)
    ulen = UInt32(len)::UInt32
    ulen > 0x0fffffff && error("CIGARStrings only support cigar operations of length 268435455")
    iszero(ulen) && error("CIGAR cannot have zero-length element")
    CIGARElement(unsafe, op, ulen)
end

function Base.show(io::IO, x::CIGARElement)
    print(io, typeof(x), '(', x.op, ", ", x.len, ')')
end

function Base.getproperty(x::CIGARElement, sym::Symbol)
    if sym === :len
        (getfield(x, :x) & 0x0fffffff) % Int
    elseif sym === :op
        reinterpret(CIGAROp, (getfield(x, :x) >>> 28) % UInt8)
    else
        error("No such field in CIGARElement: ", sym)
    end
end
Base.propertynames(x::CIGARElement) = (:len, :op)

"""
    CIGAR

A CIGAR string represents the sequence of insertions, matches and deletions
that comprise a pairwise alignment.
Construct a `CIGAR` from any object `x` where `MemoryView(x)` returns a
`MemoryView{UInt8}`, i.e. any memory-backed bytearray, or string.

See also: [`CIGARStrings.try_parse`](@ref), [`CIGARElement`](@ref)

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
struct CIGAR
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
        new(mem, n_ops, aln_len, ref_len, query_len)
    end
end

function CIGAR(x)
    y = try_parse(CIGAR, x)
    y isa CIGARError ? throw(y) : y
end

MemoryView(x::CIGAR) = x.mem

function Base.iterate(
    x::CIGAR,
    state::Int=1
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
    unreachable()
end

# TODO: Limit width for long cigars
function Base.show(io::IO, x::CIGAR)
    buf = IOBuffer()
    print(buf, summary(x), "(\"")
    write(buf, MemoryView(x))
    print(buf, "\")")
    write(io, take!(buf))
end

Base.length(x::CIGAR) = x.n_ops % Int
Base.eltype(::Type{CIGAR}) = CIGARElement

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
    len = 0
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
    CIGAR(unsafe, mem, n_ops % UInt32, aln_len % UInt32, ref_len % UInt32, query_len % UInt32)
end

function Base.print(out::IO, cigar::CIGAR)
    write(out, cigar.mem)
end

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
    consumes(x::CIGAROp) -> @NamedTuple{query::Bool, ref::Bool, aln::Bool}

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
    (;query, ref, aln)
end

@enum TranslationKind::UInt8 outside pos gap 

struct Translation
    # 4 upper bits for Kind (in case more is needed)
    x::UInt64

    Translation(::Unsafe, x::UInt64) = new(x)
end

function Translation(kind::TranslationKind, pos::Int)
    unsigned(pos) > 0x0fffffffffffffff && error() # TODO
    Translation(unsafe, kind, pos)
end

function Translation(::Unsafe, kind::TranslationKind, pos::Int)
    n = (pos % UInt64) | ((reinterpret(UInt8, kind) % UInt64) << 60)
    Translation(unsafe, n)
end

const outside_translation = Translation(outside, 0)

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

function advance(x::Anchor, e::CIGARElement)
    c = consumes(e.op)
    Anchor(
        x.query + (e.len % UInt32) * c.query,
        x.ref   + (e.len % UInt32) * c.ref,
        x.aln   + (e.len % UInt32) * c.aln
    )
end

query_length(x::CIGAR) = x.query_len % Int
ref_length(x::CIGAR) = x.ref_len % Int
aln_length(x::CIGAR) = x.aln_len % Int

function count_matches(x::CIGAR, mismatches::Int)::Int
    n_M = UInt(0)
    n_X = UInt(0)
    n_Eq = UInt(0)
    for i in x
        n_Eq += (i.op === OP_Eq) * (i.len % UInt)
        n_X += (i.op === OP_X) * (i.len % UInt)
        n_M += (i.op === OP_M) * (i.len % UInt)
    end
    if mismatches > n_M + n_X
        error() # TODO
    end
    (n_Eq % Int) + (n_M % Int - mismatches + n_X % Int)
end

function aln_identity(x::CIGAR, mismatches::Int)::Float64
    count_matches(x, mismatches) / aln_length(x)
end

# CIGARs don't support random indexing, so we need to do a linear search.
# In this function and the following, if you search for query_to_ref, then
# get_src is a function f(x::Anchor) = x.query, and get_dst f(x::Anchor) = x.ref.
# The comments will, for e.g. query_to_ref, refer to query as src and ref as dst.
function pos_to_pos(
    aln::CIGAR,
    target::Int,
    get_src::Function,
    get_dst::Function,
)::Translation
    target < 1 ? outside_translation : pos_to_pos_linear(aln, target, get_src, get_dst, zero(Anchor))
end

@inline function pos_to_pos_linear(
    CIGARElement_iter,
    target::Int,
    get_src::Function,
    get_dst::Function,
    prev_anchor::Anchor
)::Translation
    for element::CIGARElement in CIGARElement_iter
        anchor = advance(prev_anchor, element)
        if get_src(anchor) ≥ target
            # We know the src must have been incremented to at or above target by this element.
            # If the dst was also incremented, we might have overshot.
            return if get_dst(anchor) > get_dst(prev_anchor)
                position = get_dst(prev_anchor) + (target - get_src(prev_anchor))
                iszero(position) ? outside_translation : Translation(unsafe, pos, position)
            # If on the other hand the ref was not incremented, then we have hit
            # an element that increments our src (query) but not dst (ref),
            # and so we just report a gap from our previous anchor's dst.
            else
                position = Int(get_dst(prev_anchor))
                iszero(position) ? outside_translation : Translation(unsafe, gap, position)
            end
        end
        prev_anchor = anchor
    end
    outside_translation
end

function query_to_ref(x::CIGAR, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.query, i -> i.ref)
end

function query_to_aln(x::CIGAR, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.query, i -> i.aln)
end

function ref_to_query(x::CIGAR, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.ref, i -> i.query)
end

function ref_to_aln(x::CIGAR, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.ref, i -> i.aln)
end

function aln_to_query(x::CIGAR, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.aln, i -> i.query)
end

function aln_to_ref(x::CIGAR, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.aln, i -> i.ref)
end

end # module CIGARStrings
