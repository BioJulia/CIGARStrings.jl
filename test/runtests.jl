module CIGARTests

using CIGARStrings
using CIGARStrings: Errors
using Test
using MemoryViews

# Instantiation and validation
@testset "Instantiation" begin
    for good in [
            "",
            "500S",
            "10H",
            "100M",
            "5H9S1D1D1D2I9S6H",
        ]
        n_elements = count(i -> !in(i, 0x30:0x39), codeunits(good))
        for x in Any[
                view(good, 1:ncodeunits(good)),
                codeunits(good),
                collect(codeunits(good)),
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
    s = "1H4S19M33I191P9N1X22=3M1H"
    for cig in [CIGAR(s), BAMCIGAR(CIGAR(s))]
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
end

@testset "Writing" begin
    for s in [
            "",
            "15M1=9D",
            "5S1=3I9I41N13X5P16S4H",
        ]
        c = CIGAR(s)
        @test string(c) == s
        @test string(BAMCIGAR(c)) == s
        buf = IOBuffer()
        print(buf, c)
        v = take!(buf)
        m = MemoryView(c)
        @test m isa ImmutableMemoryView{UInt8}
        @test m == v
    end
end

@testset "Query / ref / aln length" begin
    c = CIGAR("150M")
    @test query_length(c) == aln_length(c) == ref_length(c) == 150
    bc = BAMCIGAR(c)
    @test c == bc # test all lengths are the same

    c = CIGAR("5M1D5M")
    @test query_length(c) == 10
    @test ref_length(c) == 11
    @test aln_length(c) == 11
    bc = BAMCIGAR(c)
    @test c == bc # test all lengths are the same

    c = CIGAR("4S2M3D2M7I10M")
    @test query_length(c) == 25
    @test ref_length(c) == 17
    @test aln_length(c) == 24
    bc = BAMCIGAR(c)
    @test c == bc # test all lengths are the same
end

@testset "Count matches" begin
    # Maximum mismtaches: 15 + 12 == 27
    # Minimum: 12
    # 15 + 7 + 12 == 34
    s = "15M9I7=3D12X"
    for c in [CIGAR(s), BAMCIGAR(CIGAR(s))]
        @test_throws InexactError count_matches(c, -1)
        @test_throws InexactError count_matches(c, typemin(Int))

        @test_throws Exception count_matches(c, 28)
        @test_throws Exception count_matches(c, 11)

        @test count_matches(c, 27) == 7
        @test count_matches(c, 20) == 14
        @test count_matches(c, 19) == 15
        @test count_matches(c, 12) == 22

        c = CIGAR("100M")
        @test count_matches(c, 0) == 100
        @test count_matches(c, 1) == 99
        @test count_matches(c, 97) == 3
        @test count_matches(c, 100) == 0
        @test_throws Exception count_matches(c, 101)
    end
end

@testset "aln identity" begin
    # Alignment length is 9 + 4 + 3 + 3 + 9 + 12 == 40
    # Maximum mismatches: 9 + 4 + 3 + 12 == 28
    # Minimum mismatches: 4
    s = "11S9M4X3D3M9I12M"
    for c in [CIGAR(s), BAMCIGAR(CIGAR(s))]
        @test_throws Exception aln_identity(c, 29)
        @test_throws Exception aln_identity(c, 3)

        @test aln_identity(c, 4) == 0.6
        @test aln_identity(c, 5) == 0.575
        @test aln_identity(c, 10) == 0.45
        @test aln_identity(c, 15) == 0.325
        @test aln_identity(c, 20) == 0.2
        @test aln_identity(c, 25) == 0.075
        @test aln_identity(c, 28) == 0.0
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
    s = "2H2S4M4I1X2D1=1I1M1I1D2M4I1X3S"
    for c in [CIGAR(s), BAMCIGAR(CIGAR(s))]

        # Getting properties of Translation object
        t = query_to_aln(c, 12)
        @test t.pos == 8
        @test t.kind == CIGARStrings.pos

        t = pos(9)
        @test t.pos == 9
        @test t.kind == CIGARStrings.pos

        t = gap(304834283)
        @test t.pos == 304834283
        @test t.kind == CIGARStrings.gap

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

    @testset "Iterator" begin
        s = "2H2S4M4I1X2D1=1I1M1I1D2M4I1X3S"
        T = CIGARStrings.PositionMapper
        for c in [CIGAR(s), BAMCIGAR(CIGAR(s))]
            it = pos_to_pos(query, aln, c, (i for i in 1:3:25))
            @test it isa T
            @test length(it) == length(1:3:25)
            @test collect(it) == [query_to_aln(c, i) for i in 1:3:25]

            v = [-10, 0, 1, 3, 8, 8, 13, 15, 21]
            for coordinate in [query, ref, aln]
                it = pos_to_pos(coordinate, coordinate, c, v)
                @test length(it) == length(v)
                @test [i.pos for i in collect(it)] == map(v) do i
                    L = if coordinate === query
                        query_length(c)
                    elseif coordinate === ref
                        ref_length(c)
                    else
                        aln_length(c)
                    end
                    in(i, 1:L) ? i : 0
                end
            end

            @test_throws ArgumentError collect(pos_to_pos(query, ref, c, [5, 9, 8]))
        end
    end
end

@testset "BAMCIGAR specifics" begin
    @testset "Construction" begin
        s = "\xd0\0\0\0\x92\0\0\0"
        @test BAMCIGAR(s) == CIGAR("13M9D")
        c = BAMCIGAR(s)
        @test MemoryView(c) == MemoryView(s)

        @test CIGARStrings.try_parse(BAMCIGAR, "abc").kind == Errors.NotModFourLength

        for (bad, err) in Any[
                ("\0\0\0\0", Errors.ZeroLength),
                ("\x10\0\0\0\x15\0\0\0\x10\0\0\0", Errors.InvalidHardClip),
                ("\x1a\0\0\0", Errors.InvalidOperation),
                ("abcdefghi", Errors.NotModFourLength),
                ("\x10\0\0\0\x14\0\0\0\x10\0\0\0", Errors.InvalidSoftClip),
            ]
            @test CIGARStrings.try_parse(BAMCIGAR, bad).kind == err
        end
    end

    @testset "Contruction between the two" begin
        for s in [
                "",
                "500S",
                "10H",
                "100M",
                "5H9S1D1D1D2I9S6H",
            ]
            c = CIGAR(s)
            bc = BAMCIGAR(c)
            c2 = CIGAR(bc)
            bc2 = BAMCIGAR(c2)

            @test c == c2
            @test bc == bc2
            @test c == bc
            @test c2 == bc2
            @test c == bc2
        end
    end

    @testset "cigar_view" begin
        v = UInt8[]
        c = CIGAR("15M1D19S9H")
        m = cigar_view!(v, BAMCIGAR(c))
        @test m isa ImmutableMemoryView{UInt8}
        @test m == v
        @test String(m) == string(CIGAR(c))
    end
end

# TODO: More equality tests
# query length
# ref length
# aln length

end # module CIGARTests
