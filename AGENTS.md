# AGENTS.md
CIGARStrings.jl provides types for representing CIGAR strings (Concise Idiosyncratic Gapped Alignment Reports) used in the SAM/BAM bioinformatics formats.

Design decisions:
* Performance is a priority in this codebase. Care about eliminating even cheap branches.
* For invalid data, all function which throw should also have non-throwing equivalents.
  This does not apply to errors from API misuse which the caller can cheaply check for
  in advance (e.g. an out-of-range length passed to `CIGARElement`). These may simply throw.
  As a deliberate exception, `normalize`, `normalize!` and `unsafe_normalize` throw on
  length overflow when merging elements, without non-throwing variants, since this only
  happens with pathological input.

## When making changes:
Do not add yourself (you, the agent) as a co-author.
Format with `runic -i .`. If it fails, alert the user and do not attempt to install Runic.
