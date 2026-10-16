# Migrating from 0.2 to 0.3
```@meta
CurrentModule = CIGARStrings
DocTestSetup = quote
    using CIGARStrings
    using MemoryViews
end
```

This page lists the breaking changes in version 0.3, and how to update code written for 0.2.
For a full list of changes, including new features, see the changelog.

## `count_matches` and `aln_identity` return `nothing`
[`count_matches`](@ref) and [`aln_identity`](@ref) no longer throw a `DomainError`
when the number of mismatches is impossible for the CIGAR. Instead, they return `nothing`.

`aln_identity` also returns `nothing` instead of `NaN` when the alignment length is zero.
To distinguish the two cases, check whether the alignment length is zero.

To restore throwing behaviour, use `something`, e.g. `something(count_matches(c, n))`,
which will throw an `ArgumentError`.

## `BAMCIGAR` is only constructed from byte arrays
`BAMCIGAR(x)` and `CIGARStrings.try_parse(BAMCIGAR, x)` now require `x` to be an
`AbstractArray{UInt8}`. Previously, `BAMCIGAR("100M")` was accepted, and silently
interpreted the bytes of the string as binary BAM data.
If you really do have binary BAM data stored in a string, use `codeunits`.

To convert a textual CIGAR string to a `BAMCIGAR`, parse it as a `CIGAR` first:
`BAMCIGAR(CIGAR("100M"))`.

## `==` compares elements, not bytes
`==` on `AbstractCIGAR`s now checks whether the two CIGARs contain the same
`CIGARElement`s in the same order. Hence, a `CIGAR` may equal a `BAMCIGAR`, and leading
zeros in a `CIGAR` are ignored. `hash` is consistent with this.
To check whether the underlying memory of two CIGARs contain the same bytes
(which is what `==` used to do), compare their `MemoryView`s, as in
`MemoryView(a) == MemoryView(b)`.

## `CIGARElement` throws `ArgumentError`
`CIGARElement(op, len)` now throws an `ArgumentError` if `len` is not in `1:268435455`,
instead of an `ErrorException` or `InexactError`.
Code that catches these errors should be updated, or better, check the length in advance.
