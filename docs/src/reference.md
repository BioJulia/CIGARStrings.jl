```@meta
CurrentModule = CIGARStrings
DocTestSetup = quote
    using CIGARStrings
    cigar = CIGAR("9H1S15M1D3X3=17M12I8S")
end
```

# CIGARStrings.jl reference

```@docs; canonical = false
CIGAR
BAMCIGAR
CIGARElement
CIGAROp
Translation
CIGARStrings.TranslationKind
CIGARStrings.CIGARError
CIGARStrings.CIGARErrorType
CIGARStrings.try_parse
```

```@docs
unsafe_switch_memory
OP_M
OP_I
OP_D
OP_S
OP_H
OP_N
OP_P
OP_X
OP_Eq
CIGARStrings.PositionMapper
query
ref
aln
ref_length
query_length
aln_length
aln_identity
count_matches
query_to_ref
query_to_aln
ref_to_query
ref_to_aln
aln_to_query
aln_to_ref
```
