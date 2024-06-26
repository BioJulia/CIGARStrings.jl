query_length(x::CIGAR) = x.query_len % Int
ref_length(x::CIGAR) = x.ref_len % Int
aln_length(x::CIGAR) = x.aln_len % Int

function count_matches(x::AbstractAlignmentReport, mismatches::Int)::Int
    n_M = UInt(0)
    n_X = UInt(0)
    n_Eq = UInt(0)
    for i in x
        n_Eq += (i.op === OP_Eq) * (i.len % UInt)
        n_X += (i.op === OP_X) * (i.len % UInt)
        n_M += (i.op === OP_M) * (i.len % UInt)
    end
    if mismatches > n_M + n_X
        error() # TODO
    end
    (n_Eq % Int) + (n_M % Int - mismatches + n_X % Int)
end

function aln_identity(x::AbstractAlignmentReport, mismatches::Int)::Float64
    count_matches(x, mismatches) / aln_length(x)
end

# CIGARs don't support random indexing, so we need to do a linear search.
# In this function and the following, if you search for query_to_ref, then
# get_src is a function f(x::Anchor) = x.query, and get_dst f(x::Anchor) = x.ref.
# The comments will, for e.g. query_to_ref, refer to query as src and ref as dst.
function pos_to_pos(
    aln::CIGAR,
    target::Int,
    get_src::Function,
    get_dst::Function,
)::Translation
    outside = Translation(Outside, target)
    target < 1 && return outside
    result = pos_to_pos_linear(aln, target, get_src, get_dst, zero(Anchor))
    isnothing(result) ? outside : result
end

@inline function pos_to_pos_linear(
    reportelement_iter,
    target::Int,
    get_src::Function,
    get_dst::Function,
    prev_anchor::Anchor
)::Union{Nothing, Translation}
    for element::ReportElement in reportelement_iter
        anchor = advance(prev_anchor, element)
        if get_src(anchor) ≥ target
            # We know the src must have been incremented to at or above target by this element.
            # If the dst was also incremented, we might have overshot.
            return if get_dst(anchor) > get_dst(prev_anchor)
                Translation(unsafe, Pos, get_dst(prev_anchor) + (target - get_src(prev_anchor)))
            # If on the other hand the ref was not incremented, then we have hit
            # an element that increments our src (query) but not dst (ref),
            # and so we just report a gap from our previous anchor's dst.
            else
                Translation(unsafe, Gap, Int(get_dst(prev_anchor)))
            end
        end
        prev_anchor = anchor
    end
    nothing
end

function pos_to_pos(
    aln::AlignmentReport,
    target::Int,
    get_src::Function,
    get_dst::Function,
)::Translation
    outside = Translation(Outside, target)
    target < 1 && return outside
    # aln is empty or target is higher than highest src in the alignment
    target > (@something get_src(last_anchor(aln)) return outside) && return outside
    # Find the first anchor with a higher value that what we search for.
    # Since all anchors are sorted in the aln, we can use binary search.
    mn, mx = 1, n_anchors(aln)
    while mx > mn
        mid = (mn + mx) >>> 1
        if target > get_src(aln.anchors[mid])
            mn = mid + 1
        else
            mx = mid
        end
    end
    # Now if e.g. anchor 3 had higher value, we know that corresponds to
    # the elements in 17-24 (for a CHUNK_SIZE = 8).
    prev_anchor = isone(mn) ? zero(Anchor) : aln.anchors[mn - 1]
    span = ((mn - 1) * CHUNK_SIZE + 1):min(mn * CHUNK_SIZE, length(aln))
    pos_to_pos_linear(
        ImmutableMemView(aln.elements)[span],
        target,
        get_src,
        get_dst,
        prev_anchor
    )
end

function query_to_ref(x::AbstractAlignmentReport, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.query, i -> i.ref)
end

function query_to_aln(x::AbstractAlignmentReport, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.query, i -> i.aln)
end

function ref_to_query(x::AbstractAlignmentReport, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.ref, i -> i.query)
end

function ref_to_aln(x::AbstractAlignmentReport, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.ref, i -> i.aln)
end

function aln_to_query(x::AbstractAlignmentReport, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.aln, i -> i.query)
end

function aln_to_ref(x::AbstractAlignmentReport, pos::Int)
    @inline pos_to_pos(x, pos, i -> i.aln, i -> i.ref)
end
