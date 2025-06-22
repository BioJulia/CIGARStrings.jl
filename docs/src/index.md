```@meta
CurrentModule = CIGARStrings
DocTestSetup = quote
    using CIGARStrings
    cigar = CIGAR("9H1S15M1D3X3=17M12I8S")
end
```

# CIGARStrings.jl
CIGARStrings.jl provide functionality for parsing and working with Concise Idiosyncratic Gapped Alignment Report - or CIGAR - strings.
CIGARs were popularized by the [SAM format](https://en.wikipedia.org/wiki/SAM_(file_format)), and are a compact run length encoding notation to represent pairwise alignments.
They can be found in the SAM, BAM, PAF, and GFA formats.

For example, the following pairwise alignment of a query to a reference:

```
    Q: TAGAT---TAGCTAC
       ||||    ||  | |
    R: TAGAACCATA--TGC
```
Can be represented by the CIGAR `5M3D2M2I3M`, representing:
1. 5 matches/mismatches
2. Then, 3 deletions
3. Then, 2 matches/mismatches
4. Then, 2 insertions
5. Finally, 3 matches/mismatches.

A CIGAR strings is always written in terms of the _query_, and not the reference. 

## Individual alignment operations
One run of identical alignment operations, e.g. "5 matches/mismatches" are represented
by a single `CIGARElement`.
Conceptually, a `CIGARElement` is an alignment operation (represented by a `CIGAROp`) and a length:

```@docs
CIGARElement
CIGAROp
```

## CIGARs
The CIGAR strings themselves are represented by the `CIGAR` type, which is backed by
a `MemoryView` (from the MemoryViews.jl package).
`CIGAR`s can be constructed from most memory-backed, contiguous byte vectors, and are validated upon construction:

```jldoctest
julia> CIGAR("2M1D3M")
CIGAR("2M1D3M")

julia> CIGAR("1M1W1S")
ERROR: Error around byte 4: Invalid operation. Possible values are "MIDNSHP=X".
[...]
```

```@docs
CIGAR
```

`CIGAR`s are iterable, and returns its `CIGARElement`s, in order:

```jldoctest
julia> collect(CIGAR("2M1D3M"))
3-element Vector{CIGARElement}:
 CIGARElement(OP_M, 2)
 CIGARElement(OP_D, 1)
 CIGARElement(OP_M, 3)
```

They can be converted back to strings using `string(::CIGAR)`, or printed into
an `IO`, in which case their normal string representation is recovered:

```jldoctest
julia> c = CIGAR("2M1D3M");

julia> string(c)
"2M1D3M"

julia> io = IOBuffer(); print(io, c); String(take!(io))
"2M1D3M"
```

### Basic information about a CIGAR
The reference, query and alignment length can be obtained with the functions
[`ref_length`](@ref), [`query_length`](@ref) and [`aln_length`](@ref).

In the alignment below, represented as `5M3D2M2I3M`, the query length is
12, since there are 12 query symbols, the reference length is 13, and the
alignment length is 15.

```
    Q: TAGAT---TAGCTAC
       ||||    ||  | |
    R: TAGAACCATA--TGC
```

We always have `aln_length(c) ≥ max(query_length(c), ref_length(c))`

```jldoctest
julia> c = CIGAR("5M3D2M2I3M");

julia> ref_length(c)
13

julia> query_length(c)
12

julia> aln_length(c)
15
```

The alignment identity (number of matches, not mismatches divided by alignment length)
can be obtained with [`aln_identity`](@ref).
Since the number of mismatches may not be known from the CIGAR itself (i.e. "100M"
could have anywhere from 0 to 100 mismatches), the mismatches is passed in as an argument:

```jldoctest
julia> aln_identity(CIGAR("3M1D3M1I2M"), 2)
0.6
```

## Position translation
Sometimes it may be necessary to answer questions of the form
"which reference position does query position 8 align to?".

As an example, consider the alignment below.
The query position (QP), reference position (RP) and alignment position (AP)
are also written in this alignment.

```
    QP:12345   6789012

    Q: TAGAT---TAGCTAC
       ||||    ||  | |
    R: TAGAACCATA--TGC

    RP:1234567890  123
    AP:123456789012345
```

We can see that query position 6 aligns to reference position 9, which is also
alignment position 9.

These position translations can be obtained using the functions [`query_to_ref`](@ref),
[`query_to_aln`](@ref), [`ref_to_query`](@ref), [`ref_to_aln`](@ref), [`aln_to_query`](@ref) and [`aln_to_ref`](@ref).
They return a `Translation` object that contains two properties: `.pos` and `.kind`.

When a position translation has a straightforward answer, the `.kind` property is
`CIGARStrings.pos`, and the `.pos` field is the corresponding position:

```jldoctest
julia> c = CIGAR("4M3D2M2I3M"); # alignment above

julia> query_to_ref(c, 6)
Translation(pos, 9)

julia> aln_to_query(c, 9)
Translation(pos, 6)
```

```@docs
Translation
CIGARStrings.TranslationKind
```

## Errors and error recovery
CIGARStrings.jl allows you to parse a poential CIGAR string without throwing an exception if the data is invalid, using the function [`CIGARStrings.try_parse`](@ref).

```@docs
CIGARStrings.CIGARError
CIGARStrings.CIGARErrorType
CIGARStrings.try_parse
```