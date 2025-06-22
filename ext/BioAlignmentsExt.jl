module BioAlignmentsExt

# TODO: A bunch of stuff is subtly wrong or weird in BioAlignments, so I've
# decided to propone implementing this.
# For example: 
# * Alignment starting position is one-indexed when constructing, but shown
# zero-indexed.
# * Soft and hard clips advances the alignment in an `Alignment`, which is wrong.

using BioAlignments: BioAlignments as BA
import CIGARStrings: CIGAR, Anchor, advance, CIGAROp

function BA.AlignmentAnchor(anchor::Anchor, op::CIGAROp)
    return BA.AlignmentAnchor(anchor.query, anchor.ref, anchor.aln, reinterpret(BA.Operation, op))
end

"""
    Alignment(cigar::CIGAR, refpos::Int)::Alignment

Create an `Alignment` type from a CIGAR. The CIGAR's alignment starting position
in the reference is given as `refpos`.

!!! NOTE
    Due to limitations of the current BioAlignments

# Examples
```jldoctest
julia> c = CIGAR("19M2D5M1I22M5S");

julia> BioAlignments.Alignment(c, 1005)
Alignment:
  aligned range:
    seq: 0-52
    ref: 1004-1052
  CIGAR string: 19M2D5M1I22M
```
"""
function BA.Alignment(cigar::CIGAR, refpos::Int)
    if refpos < 1
        error("Alignment cannot start before ref position 1")
    end
    anchor = Anchor(0, refpos - 1, 0)
    anchors = Vector{BA.AlignmentAnchor}(undef, length(cigar) + 1)
    anchors[1] = (BA.AlignmentAnchor(0, anchor.ref, 0, BA.OP_START))
    for (i, element) in enumerate(cigar)
        anchor = advance(anchor, element)
        anchors[i + 1] = BA.AlignmentAnchor(anchor, element.op)
    end
    # Note: CIGAR constructor checks everything that Alignment needs.
    return BA.Alignment(anchors, false)
end

end # module BioAlignmentsExt
