module BioAlignmentsExt

using BioAlignments: BioAlignments as BA
import CIGARStrings: CIGAR, Anchor, advance, CIGAROp

function BA.AlignmentAnchor(anchor::Anchor, op::CIGAROp)
    BA.AlignmentAnchor(anchor.query, anchor.ref, anchor.aln, reinterpret(BA.Operation, op))
end

function BA.Alignment(cigar::CIGAR, refpos::Int)
    if refpos < 1
        error("Alignment cannot start before ref position 1")
    end
    anchor = Anchor(0, refpos - 1, 0)
    anchors = Vector{BA.AlignmentAnchor}(undef, length(cigar) + 1)
    anchors[1] = (BA.AlignmentAnchor(0, anchor.ref, 0, BA.OP_START))
    for (i, element) in enumerate(cigar)
        anchor = advance(anchor, element)
        anchors[i+1] = BA.AlignmentAnchor(anchor, element.op)
    end
    # Note: CIGAR constructor checks everything that Alignment needs.
    BA.Alignment(anchors, false)
end

end # module BioAlignmentsExt
