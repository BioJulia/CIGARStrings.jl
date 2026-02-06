# AGENTS.md

This file provides guidance to AI agents when working with code in this repository.

## Build and Test Commands

```bash
# Run all tests
JULIA_TEST_FAILFAST=true julia --project -e 'using Pkg; Pkg.test()'

# Run tests with specific test file
JULIA_TEST_FAILFAST=true julia --project test/runtests.jl

# Run Julia REPL with project activated
julia --project

# Generate documentation locally
julia --project=docs docs/make.jl
```

## Architecture

CIGARStrings.jl provides types for representing CIGAR strings (Concise Idiosyncratic Gapped Alignment Reports) used in the SAM/BAM bioinformatics formats.

### Core Types

- **`AbstractCIGAR`** - Abstract supertype for CIGAR representations
- **`CIGAR`** - ASCII string representation (e.g., "15M1D10M")
- **`BAMCIGAR`** - Compact 32-bit integer representation used in BAM files
- **`CIGAROp`** - 1-byte primitive representing operations (M, I, D, N, S, H, P, =, X)
- **`CIGARElement`** - Operation + length pair (upper 28 bits = length, lower 4 bits = op)
- **`Translation`** - Result of position mapping between coordinate systems

### File Structure

- `src/CIGARStrings.jl` - Main module: exports, error handling, CIGAROp/CIGARElement definitions, position translation logic, `is_compatible`
- `src/bytecigar.jl` - `CIGAR` type: ASCII parsing/iteration/serialization
- `src/bamcigar.jl` - `BAMCIGAR` type: binary parsing/iteration, conversion to/from CIGAR

### Key Design Patterns

- Both CIGAR types store precomputed `aln_len`, `ref_len`, `query_len` for O(1) length queries
- Parsing uses `try_parse` returning `Union{T, CIGARError}` pattern (no exceptions on invalid input)
- Memory is backed by `ImmutableMemoryView{UInt8}` from MemoryViews.jl
- Internal `Unsafe` sentinel type for constructors that skip validation
- Position translation via `pos_to_pos(from, to, cigar, positions)` returns lazy `PositionMapper` iterator

### Coordinate Systems

Three coordinate systems: `query`, `ref`, `aln` (alignment with gaps). Position functions like `query_to_ref`, `ref_to_aln` etc. translate between them, returning `Translation` objects with `.kind` (outside/pos/gap) and `.pos`.
