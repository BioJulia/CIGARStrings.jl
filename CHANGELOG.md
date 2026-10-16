# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.3.0
### Breaking changes
* `count_matches` now returns `nothing` instead of throwing a `DomainError` when
  the number of mismatches is impossible for the CIGAR. Its return type is now
  `Union{Nothing, Int}`.
* `aln_identity` now returns `nothing` instead of throwing when the number of
  mismatches is impossible, and instead of returning `NaN` when the alignment
  length is zero. Its return type is now `Union{Nothing, Float64}`.
* `BAMCIGAR(x)` and `CIGARStrings.try_parse(BAMCIGAR, x)` now only accept
  `x::AbstractArray{UInt8}`. Previously, any `x` that could be converted to a
  `MemoryView` was accepted. This means that e.g. `BAMCIGAR("100M")` no longer
  silently interprets the string's bytes as binary BAM data.
* `CIGARElement(op, len)` now throws an `ArgumentError` when `len` is out of range,
  instead of an `ErrorException` or `InexactError`.
* `==` on `AbstractCIGAR`s now compares the elements, instead of the underlying bytes.
  As a result, CIGARs differing only in leading zeros, such as `CIGAR("01M")` and
  `CIGAR("1M")`, now compare equal.
* `consumes` is now exported.

### New features
* New function `reverse(::AbstractCIGAR)`, which returns a new CIGAR with the
  elements in reverse order.
* New function `unsafe_reverse!`, which reverses the elements of a CIGAR in-place.
* `CIGAR` and `BAMCIGAR` can now be constructed from an `AbstractVector{CIGARElement}`.
* New public function `CIGARStrings.try_from_elements`, a non-throwing version of
  constructing a CIGAR from an `AbstractVector{CIGARElement}`.
* `Char(::CIGAROp)` returns the character of the operation, e.g. `Char(OP_M) == 'M'`.
* `tryparse(CIGAROp, ::Char)` parses an operation from its character, returning
  `nothing` if the character is not a valid CIGAR operation.
* `print(io, ::CIGARElement)` writes the element in CIGAR format, e.g. `5M`, so that
  `string(CIGARElement(OP_M, 5)) == "5M"`.
* `aln_identity` now accepts any `Integer` for the number of mismatches,
  not only `Int`.

### Bug fixes
* `hash` is now consistent with `==` for `AbstractCIGAR`s. Previously, equal
  CIGARs could have different hashes, e.g. a `CIGAR` and an equal `BAMCIGAR`,
  breaking their use in `Dict` and `Set`.
* `normalize`, `normalize!` and `unsafe_normalize` on a `CIGAR` now remove
  leading zeros from element lengths, e.g. `01M` becomes `1M`.
* `BAMCIGAR` data is now consistently read and written as little-endian.
  This fixes `count_matches` and `normalize!` for `BAMCIGAR` on big-endian platforms.
* The index of the `CIGARError` thrown on length overflow by `normalize!` on a
  `BAMCIGAR` now correctly points to the element in the input.

### Performance
* Iterating a `PositionMapper` (returned by `pos_to_pos` with multiple positions) and
  an `AlignedBlocks` (returned by `aligned_blocks`) is significantly faster.

### Documentation
* Many docstring and manual fixes, including corrections to wrong examples and
  to the description of the `N`, `S` and `H` operations.
