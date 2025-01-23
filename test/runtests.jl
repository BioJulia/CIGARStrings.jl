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
    pos(n) = CIGARStrings.Translation(CIGARStrings.pos, n)
    gap(n) = CIGARStrings.Translation(CIGARStrings.gap, n)

    #       1234    5678 9 012    3
    # R XXXXXXXX----XXXX-X-XXX----XXXX
    # Q HHSS|||......--|.|.-||.....SSS
    #   1234567890123  4567 8901234567
    # A     12345678901234567890123
    c = CIGAR("2H2S4M4I1X2D1=1I1M1I1D2M4I1X3S")
    
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
    @test is_outside(query_to_ref(c, 25))
    @test is_outside(query_to_aln(c, 25))
    @test is_outside(aln_to_query(c, 24))
    @test is_outside(aln_to_ref(c, 24))
    @test is_outside(ref_to_query(c, 14))
    @test is_outside(ref_to_aln(c, 14))
    
    # Within alignment
    @test query_to_ref(c, 5) == pos(1)
    @test query_to_aln(c, 5) == pos(1)
    @test ref_to_query(c, 1) == pos(5)
    @test ref_to_aln(c, 1) == pos(1)
    @test aln_to_query(c, 1) == pos(5)
    @test aln_to_ref(c, 1) == pos(1)

    # Various
    @test query_to_ref(c, 7) == pos(3)
    @test query_to_ref(c, 8) == pos(4)
    @test query_to_ref(c, 9) == gap(4)
    @test query_to_ref(c, 12) == gap(4)
    @test query_to_ref(c, 13) == pos(5)
    @test query_to_ref(c, 14) == pos(8)
    @test query_to_ref(c, 15) == gap(8)
    @test query_to_ref(c, 16) == pos(9)
    @test query_to_ref(c, 17) == gap(9)
    @test query_to_ref(c, 18) == pos(11)
    @test query_to_ref(c, 19) == pos(12)
    @test query_to_ref(c, 20) == gap(12)
    @test query_to_ref(c, 23) == gap(12)
    @test query_to_ref(c, 24) == pos(13)

    @test query_to_aln(c, 9) == pos(5)
    @test query_to_aln(c, 13) == pos(9)
    @test query_to_aln(c, 14) == pos(12)
    @test query_to_aln(c, 17) == pos(15)
    @test query_to_aln(c, 18) == pos(17)
    @test query_to_aln(c, 19) == pos(18)
    @test query_to_aln(c, 24) == pos(23)

    @test aln_to_query(c, 9) == pos(13)
    @test aln_to_query(c, 10) == gap(13)
    @test aln_to_query(c, 11) == gap(13)
    @test aln_to_query(c, 12) == pos(14)
    @test aln_to_query(c, 13) == pos(15)
    @test aln_to_query(c, 14) == pos(16)
    @test aln_to_query(c, 15) == pos(17)
    @test aln_to_query(c, 16) == gap(17)
    @test aln_to_query(c, 17) == pos(18)
    @test aln_to_query(c, 18) == pos(19)
    @test aln_to_query(c, 19) == pos(20)
    @test aln_to_query(c, 23) == pos(24)

    @test ref_to_query(c, 5) == pos(13)
    @test ref_to_query(c, 6) == gap(13)
    @test ref_to_query(c, 7) == gap(13)
    @test ref_to_query(c, 8) == pos(14)
    @test ref_to_query(c, 9) == pos(16)
    @test ref_to_query(c, 10) == gap(17)
    @test ref_to_query(c, 11) == pos(18)
    @test ref_to_query(c, 13) == pos(24)

    @test ref_to_aln(c, 1) == pos(1)
    @test ref_to_aln(c, 4) == pos(4)
    @test ref_to_aln(c, 5) == pos(9)
    @test ref_to_aln(c, 8) == pos(12)
    @test ref_to_aln(c, 9) == pos(14)
    @test ref_to_aln(c, 10) == pos(16)
    @test ref_to_aln(c, 12) == pos(18)
    @test ref_to_aln(c, 13) == pos(23)
end

@testset "Conversion to alignment" begin
    s = "3S1M2D5I5X1H"
    cig = CIGAR(s)
    aln = Alignment(cig, 5)
    @test BioAlignments.cigar(aln) == s[3:end-2]
end

end # module CIGARTests
