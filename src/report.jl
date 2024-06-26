const CHUNK_SIZE = 8

struct AlignmentReport <: AbstractAlignmentReport
    elements::Memory{ReportElement}
    anchors::Memory{Anchor}
    len::Int
end
Base.length(x::AlignmentReport) = x.len

function last_anchor(x::AlignmentReport)
    iszero(length(x)) && return nothing
    x.anchors[n_anchors(x)] # TODO inbounds
end
n_anchors(x::AlignmentReport) = cld(length(x) % UInt, CHUNK_SIZE) % Int

function Base.iterate(aln::AlignmentReport, i::Int=1)
    if (i - 1) % UInt < length(aln)
        e = @inbounds aln.elements[i]
        (e, i + 1)
    else
        nothing
    end
end

# TODO: Protect against overflow in anchor. How?
function AlignmentReport(it)
    len = Int(length(it))::Int
    (last_element::ReportElement, rest) = let
        y = Iterators.peel(it)
        if isnothing(y)
            return AlignmentReport(
                Memory{ReportElement}(),
                Memory{Anchor}(),
                0,
            )
        else
            y
        end
    end
    elements = Memory{ReportElement}(undef, len)
    anchors = Memory{Anchor}(undef, cld(len % UInt, CHUNK_SIZE) % Int)
    anchor = zero(Anchor)
    report_len = 0
    for element::ReportElement in rest
        if element.op == last_element.op
            last_element = ReportElement(element.op, element.len + last_element.len)
            continue
        end
        anchor = advance(anchor, last_element)
        report_len += 1
        elements[report_len] = last_element # this cannot be inbounds
        (anchor_index, rm) = divrem(report_len % UInt, CHUNK_SIZE)
        if iszero(rm)
            anchors[anchor_index] = anchor
        end
        last_element = element
    end
    report_len += 1
    elements[report_len] = last_element
    (anchor_offset, rm) = divrem(report_len % UInt, CHUNK_SIZE)
    anchors[anchor_offset + 1 - iszero(rm)] = advance(anchor, last_element)
    AlignmentReport(elements, anchors, report_len)
end