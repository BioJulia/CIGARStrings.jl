struct CIGAR <: AbstractAlignmentReport
    mem::ImmutableMemView{UInt8}
    n_ops::UInt32
    aln_len::UInt32
    ref_len::UInt32
    query_len::UInt32

    function CIGAR(
        ::Unsafe,
        mem::ImmutableMemView{UInt8},
        n_ops::UInt32,
        aln_len::UInt32,
        ref_len::UInt32,
        query_len::UInt32,
    )
        new(mem, n_ops, aln_len, ref_len, query_len)
    end
end

function CIGAR(x)
    y = try_parse(CIGAR, x)
    isnothing(y) ? error("Invalid CIGAR string") : y
end

function Base.iterate(
    x::CIGAR,
    state::Int=1
)::Union{Tuple{ReportElement, Int}, Nothing}
    mem = x.mem
    len = length(mem)
    state > len && return nothing
    n = 0
    @inbounds for i in state:len
        b = mem[i] - 0x30
        if b < 0x0a
            n = (10 * n) + b
        else
            b -= UInt8('=') - UInt8('0')
            enc = ((LUT >>> ((4 * b) & 127)) % UInt8) & 0x0f
            op = reinterpret(CIGAROp, enc)
            return (ReportElement(unsafe, op, n % UInt32), i + 1)
        end
    end
    unreachable()
end

# TODO: Limit width for long cigars
Base.show(io::IO, x::CIGAR) = print(io, summary(x), "(\"", String(x.mem), "\")")
Base.length(x::CIGAR) = x.n_ops % Int
Base.eltype(::Type{<:CIGAR}) = ReportElement

function try_parse(::Type{CIGAR}, x)::Union{Nothing, CIGAR}
    mem = ImmutableMemView(x)::ImmutableMemView{UInt8}
    # H must be either first or last
    # S must be either first or last, or preceded or succeeded by H
    n_ops = 0
    aln_len = 0
    ref_len = 0
    query_len = 0
    is_first = true
    must_be_last = false
    last_was_H = false
    next_must_be_H = false
    last_was_num = false
    len = 0
    n = 0
    for b in mem
        b -= 0x30
        last_was_num = b < 0x0a
        if last_was_num
            n = (10 * n) + b
            n > 0x0fffffff && return nothing
        else
            b -= UInt8('=') - UInt8('0')
            b > (UInt8('X') - UInt8('=')) && return nothing
            enc = ((LUT >>> ((4 * b) & 127)) % UInt8) & 0x0f
            enc == 0x0f && return nothing
            iszero(n) && return nothing
            op = reinterpret(CIGAROp, enc)
            must_be_last && return nothing
            if op == OP_H
                must_be_last |= !is_first
            else
                next_must_be_H && return nothing
            end
            if op == OP_S
                next_must_be_H |= !(is_first | last_was_H)
            end
            last_was_H = op == OP_H
            is_first = false
            c = consumes(op)
            n_ops += 1
            ref_len += c.ref * n
            query_len += c.query * n
            aln_len += (c.ref | c.query) * n
            n = 0
        end
    end
    last_was_num && return nothing
    max(aln_len, n_ops) > typemax(UInt32) && return nothing
    CIGAR(unsafe, mem, n_ops % UInt32, aln_len % UInt32, ref_len % UInt32, query_len % UInt32)
end