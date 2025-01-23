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
```


* Individual CIGAR operations

* CIGARs

* Errors

* BioAlignments

## Reference
```@autodocs
Modules = [CIGARStrings]
Order   = [:type, :function]
```
