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

## Individual alignment operations
The individual alignment operations, e.g. "5 matches/mismatches" are represented by
`CIGARElement`. Conceptually, a `CIGARElement` is an alignment operation (represented by a `CIGAROp`) and a length:

```@docs
CIGARElement
CIGAROp
OP_M
OP_I
OP_D
OP_N
OP_S
OP_H
OP_P
OP_Eq
OP_X
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

### Basic statistics of CIGARs
The reference, query and alignment length can be obtained with the functions
`ref_length`, `query_length` and `aln_length`.

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

```@docs
ref_length
query_length
aln_length
```

The alignment identity (number of matches, not mismatches divided by alignment length)
can be obtained with `aln_identity`.
Since the number of mismatches may not be known from the CIGAR itself (i.e. "100M"
could have anywhere from 0 to 100 mismatches), the mismatches is passed in as an argument:

```jldoctest
julia> aln_identity(CIGAR("3M1D3M1I2M"), 2)
0.6
```

```@docs
aln_identity
```

## Position translation

* query_to_ref, query_to_aln, ref_to_query, ref_to_aln,
    aln_to_query, aln_to_ref

```@docs
Translation
CIGARStrings.TranslationKind
```

## Errors and error recovery
```@docs
CIGARStrings.CIGARErrorType
CIGARStrings.CIGARError
CIGARStrings.try_parse
```

## BioAlignments interoperability

