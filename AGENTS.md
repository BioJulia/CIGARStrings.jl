# AGENTS.md
CIGARStrings.jl provides types for representing CIGAR strings (Concise Idiosyncratic Gapped Alignment Reports) used in the SAM/BAM bioinformatics formats.

Design decisions:
* Performance is a priority in this codebase. Care about eliminating even cheap branches.
* For invalid data, all function which throw should also have non-throwing equivalents.

## When making changes:
Do not add yourself (you, the agent) as a co-author.
Format with `runic -i .`. If it fails, alert the user and do not attempt to install Runic.
