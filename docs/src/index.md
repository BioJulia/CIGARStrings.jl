```@meta
CurrentModule = CIGARStrings
DocTestSetup = quote
    using CIGARStrings
    using MemoryViews
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
A CIGAR string is represented by an `AbstractCIGAR`, which currently has two subtypes: `CIGAR` and `BAMCIGAR`.
These types differ in their memory layout: The former stores the CIGAR as its ASCII representation (as used in the SAM format), and the latter stores it in a binary format (as used in the BAM format).
Both typs store its underlying data as an `ImmutableMemoryView{UInt8}`.

```@docs
AbstractCIGAR
```

The API for these two types are almost interchangeable, so examples below will use `CIGAR`, since its plaintext representation makes examples easier.
See [BAMCIGAR section](@ref bamcigar) for a list of all differences between the two types.

CIGAR strings are validated upon construction

```jldoctest
julia> CIGAR("2M1D3M")
CIGAR("2M1D3M")

julia> CIGAR("1M1W1S")
ERROR: Error around byte 4: Invalid operation. Possible values are "MIDNSHP=X".
[...]
```

Since CIGAR strings occur in various bioinformatics file formats, it is expected
that users of CIGARStrings.jl will construct `CIGAR`s from a view into a buffer storing a chunk of the file.

This is zero-copy, and _ought_ not to allocate (although it currently does, due to [a limitation of the Julia compiler](https://github.com/JuliaLang/julia/issues/53584)):

For example:
```jldoctest
julia> data = b"some format with CIGAR: 15M9D18M etc";

julia> buffer = collect(data);

julia> c = CIGAR(view(buffer, 25:32))
CIGAR("15M9D18M")
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

The memory underlying the CIGAR types can be obtained with `MemoryView(x)` using the MemoryViews.jl package:

```jldoctest
julia> c = CIGAR("19M");

julia> println(MemoryView(c))
UInt8[0x31, 0x39, 0x4d]

julia> println(MemoryView(BAMCIGAR(c)))
UInt8[0x30, 0x01, 0x00, 0x00]
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

## Comparing CIGARs
When comparing `CIGAR`s using `==`, it will check if the `CIGAR`s are literally identical, in the
sense that they are composed of the same bytes:

```jldoctest compare
julia> c1 = CIGAR("10M");

julia> c2 = CIGAR("4=1X5=");

julia> c3 = CIGAR("10M");

julia> c1 == c2
false

julia> c1 == c3
true
```
However, in the above example, since the CIGAR operation `M` signifies a match or a mismatch, all three
CIGARs are indeed compatible, since `10M` is also a valid CIGAR annotation for the same alignment
as `4=1X5=`.

This notion of compatibility tested with `is_compatible`:

```jldoctest compare
julia> is_compatible(c1, c2)
true

julia> is_compatible(CIGAR("1X1M"), CIGAR("1=1M"))
false
```

```@docs
is_compatible
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

These position translation can be obtained using the function [`pos_to_pos`](@ref),
specifying the source and destination coordinate systems [`query`](@ref), [`ref`](@ref)
or [`aln`](@ref).
When passed an integer, this function returns `Translation` object that contains two properties: `.pos` and `.kind`.

When a position translation has a straightforward answer, the `.kind` property is
`CIGARStrings.pos`, and the `.pos` field is the corresponding position:

```jldoctest
julia> c = CIGAR("4M3D2M2I3M"); # alignment above

julia> pos_to_pos(query, ref, c, 6)
Translation(pos, 9)

julia> pos_to_pos(aln, query, c, 9)
Translation(pos, 6)
```

Note that these operations are in __linear time__, as they scan the CIGAR string from the beginning.

To efficiently query multiple translations in the same scan of the CIGAR string, you can pass a sorted (ascending) iterator of integers.
In this case, `pos_to_pos` will return a lazy iterator of `Pair{Int, Translation}`, representing `source_index => destination_index`:

```jldoctest
julia> c = CIGAR("4M3D2M2I3M"); # alignment above

julia> it = pos_to_pos(query, ref, c, [1, 5, 6, 11]);

julia> length(it)
4

julia> collect(it)
4-element Vector{Pair{Int64, Translation}}:
  1 => Translation(pos, 1)
  5 => Translation(pos, 8)
  6 => Translation(pos, 9)
 11 => Translation(pos, 12)
```

```@docs
pos_to_pos
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

## [The BAMCIGAR type](@id bamcigar)
A CIGAR in the BAM format is stored in an array of 32-bit integers.
However, in order to make zero-copy CIGARs possible, the `BAMCIGAR` type is backed by an `ImmutableMemoryView{UInt8}` instead, with the same memory layout as its equivalent `Memory{UInt32}`.

```@docs
CIGARStrings.BAMCIGAR
CIGARStrings.BAMCIGAR(::MutableMemoryView{UInt8}, ::CIGAR)
```

A `BAMCIGAR` can be constructed from its binary representation, using any type which implements `MemoryViews.MemoryView`:

```jldoctest
julia> BAMCIGAR("\x54\4\0\0\x70\4\0\0")
BAMCIGAR(CIGAR("69S71M"))
```

This is not zero-cost: Like `CIGAR` the type contains some metadata and is validated upon construction.

Like `CIGAR`, the `try_parse` function can be used:

```jldoctest
julia> CIGARStrings.try_parse(BAMCIGAR, "\x5f\4\0\0\x70\4\0\0")
CIGARStrings.CIGARError(1, CIGARStrings.Errors.InvalidOperation)
```

`CIGAR` and `BAMCIGAR` can be converted ifallably to each other:

```jldoctest
julia> c = CIGAR("6H19S18M1I22=8I2S");

julia> b = BAMCIGAR(c);

julia> b == c
true

julia> CIGAR(b) == c
true
```

Note that printing a `BAMCIGAR` allocates, because it needs to allocate a new piece of memory to store its ASCII representation.
For high performance applications, the function `cigar_view!` may be used:

```@docs
CIGARStrings.cigar_view!
```
