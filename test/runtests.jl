module CIGARTests

using CIGARStrings
using CIGARStrings: Errors
using BioAlignments
using Test
using MemoryViews

# Instantiation and validation
@testset "Instantiation" begin
    for good in [
        "",
        "500S",
        "10H",
        "100M",
        "5H9S1D1D1D2I9S6H"
    ]
        n_elements = count(i -> !in(i, 0x30:0x39), codeunits(good))
        for x in Any[
            view(good, 1:ncodeunits(good)),
            codeunits(good),
            collect(codeunits(good))
        ]
            c = CIGAR(x)
            @test length(collect(c)) == n_elements
        end
    end

    for (bad, err) in Any[
        ("H", Errors.ZeroLength),
        ("1S1H1D", Errors.InvalidHardClip),
        ("1K", Errors.InvalidOperation),
        ("11111", Errors.Truncated),
        ("1HM", Errors.ZeroLength),
        ("1D1S1I1H", Errors.InvalidSoftClip),
        ("219382982739847498327912M", Errors.IntegerOverflow),
        ("1H~~", Errors.InvalidOperation),
        ("ÆÅA", Errors.InvalidOperation),
    ]
        @test CIGARStrings.try_parse(CIGAR, bad).kind == err
    end
end

@testset "Iteration" begin
    cig = CIGAR("1H4S19M33I191P9N1X22=3M1H")
    elements = collect(cig)
    @test elements == [
        CIGARElement(OP_H, 1),
        CIGARElement(OP_S, 4),
        CIGARElement(OP_M, 19),
        CIGARElement(OP_I, 33),
        CIGARElement(OP_P, 191),
        CIGARElement(OP_N, 9),
        CIGARElement(OP_X, 1),
        CIGARElement(OP_Eq, 22),
        CIGARElement(OP_M, 3),
        CIGARElement(OP_H, 1),
    ]
end

@testset "Writing" begin
    for s in [
        "",
        "15M1=9D",
        "5S1=3I9I41N13X5P16S4H",
    ]
        c = CIGAR(s)
        @test string(c) == s
        buf = IOBuffer()
        print(buf, c)
        v = take!(buf)
        m = MemoryView(c)
        @test m isa ImmutableMemoryView{UInt8}
        @test m == v
    end
end

@testset "Translate positions" begin
    is_outside(x) = x === CIGARStrings.outside_translation

    #       1234    5678 9 012    3456
    # R XXXXXXXX----XXXX-X-XXX----XXXX
    # Q HHSS|||......--|.|.-||.....SSS
    #   1234567890123  4567 8901234567
    # A     123456789 123456789 123
    c = CIGAR("2H2S4M4I1X2D1=1I1M1I1X1M4I3S")
    
    # Before alignment
    @test is_outside(query_to_aln(c, 0))
    @test is_outside(ref_to_aln(c, 0))
    @test is_outside(ref_to_query(c, 0))
    @test is_outside(query_to_ref(c, 0))
    @test is_outside(aln_to_query(c, 0))
    @test is_outside(aln_to_ref(c, 0))
    @test is_outside(query_to_aln(c, 1))
    @test is_outside(query_to_aln(c, 4))
    @test is_outside(query_to_ref(c, 1))
    @test is_outside(query_to_ref(c, 4))
    
    # After alignment 
    # TODO: Document how this works exactly
end

# TODO
@testset "Conversion to alignment" begin
    cig = CIGAR("")
    aln = Alignment(cig, 1)
end


end # module CIGARTests
