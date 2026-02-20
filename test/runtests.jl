module CIGARTests

using CIGARStrings
using CIGARStrings: Errors
using Test
using MemoryViews
using LightBoundsErrors: LightBoundsError

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

    @testset "try_parse catches total length overflow" begin
        maxop = "268435455"

        # Ref length overflows UInt32 while aln length stays tiny.
        ref_overflow = repeat(maxop * "N", 17) * "1M"
        err = CIGARStrings.try_parse(CIGAR, ref_overflow)
        @test err isa CIGARStrings.CIGARError
        @test err.kind == Errors.IntegerOverflow

        # Query length overflows UInt32 while aln length is still <= UInt32 max.
        query_overflow = maxop * "S" * repeat(maxop * "I", 14) * "16I" * maxop * "S"
        err = CIGARStrings.try_parse(CIGAR, query_overflow)
        @test err isa CIGARStrings.CIGARError
        @test err.kind == Errors.IntegerOverflow
    end

    @testset "unsafe_switch_memory" begin
        for str in [
                "",
                "500S",
                "10H",
                "100M",
                "5H9S1D1D1D2I9S6H",
            ]
            padded = collect(codeunits(str))
            pushfirst!(padded, 0x01, 0x02, 0x03)
            push!(padded, 0x61, 0x62, 0x63)
            vw = MemoryView(padded[4:(end - 3)])
            cigar = CIGAR(vw)
            for T in [CIGAR, BAMCIGAR]
                c = T === CIGAR ? cigar : BAMCIGAR(cigar)
                cp = collect(MemoryView(c))
                c2 = unsafe_switch_memory(c, ImmutableMemoryView(cp))
                @test MemoryView(c) == MemoryView(c2)
                @test c == c2
                if !isempty(MemoryView(c))
                    @test MemoryView(c) !== MemoryView(c2)
                end
            end
        end
    end

    s = "1H4S19M33I191P9N1X22=3M1H"
    for cig in [CIGAR(s), BAMCIGAR(CIGAR(s))]
        cig2 = copy(cig)
        @test cig == cig2
        @test MemoryView(cig) !== MemoryView(cig2)
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
    #     12345678901  2345 6789012345
    # A     12345678901234567890123
    s = "2H2S4M4I1X2D1=1I1M1I1D2M4I1X3S"
    for c in [CIGAR(s), BAMCIGAR(CIGAR(s))]

        # Getting properties of Translation object
        t = query_to_aln(c, 12)
        @test t.pos == 12
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
        @test is_outside(query_to_aln(c, 2))
        @test is_outside(query_to_ref(c, 1))
        @test is_outside(query_to_ref(c, 2))

        # After alignment
        @test is_outside(query_to_ref(c, 23))
        @test is_outside(query_to_aln(c, 23))
        @test is_outside(aln_to_query(c, 24))
        @test is_outside(aln_to_ref(c, 24))
        @test is_outside(ref_to_query(c, 14))
        @test is_outside(ref_to_aln(c, 14))

        # Within alignment
        @test query_to_ref(c, 3) == pos(1)
        @test query_to_aln(c, 3) == pos(1)
        @test ref_to_query(c, 1) == pos(3)
        @test ref_to_aln(c, 1) == pos(1)
        @test aln_to_query(c, 1) == pos(3)
        @test aln_to_ref(c, 1) == pos(1)

        # Various
        @test query_to_ref(c, 7) == gap(4)
        @test query_to_ref(c, 8) == gap(4)
        @test query_to_ref(c, 9) == gap(4)
        @test query_to_ref(c, 12) == pos(8)
        @test query_to_ref(c, 13) == gap(8)
        @test query_to_ref(c, 14) == pos(9)
        @test query_to_ref(c, 15) == gap(9)
        @test query_to_ref(c, 16) == pos(11)
        @test query_to_ref(c, 17) == pos(12)
        @test query_to_ref(c, 18) == gap(12)
        @test query_to_ref(c, 19) == gap(12)
        @test query_to_ref(c, 20) == gap(12)
        @test query_to_ref(c, 22) == pos(13)
        @test is_outside(query_to_ref(c, 23))

        @test query_to_aln(c, 9) == pos(7)
        @test query_to_aln(c, 13) == pos(13)
        @test query_to_aln(c, 14) == pos(14)
        @test query_to_aln(c, 17) == pos(18)
        @test query_to_aln(c, 18) == pos(19)
        @test query_to_aln(c, 19) == pos(20)
        @test query_to_aln(c, 22) == pos(23)

        @test aln_to_query(c, 9) == pos(11)
        @test aln_to_query(c, 10) == gap(11)
        @test aln_to_query(c, 11) == gap(11)
        @test aln_to_query(c, 12) == pos(12)
        @test aln_to_query(c, 13) == pos(13)
        @test aln_to_query(c, 14) == pos(14)
        @test aln_to_query(c, 15) == pos(15)
        @test aln_to_query(c, 16) == gap(15)
        @test aln_to_query(c, 17) == pos(16)
        @test aln_to_query(c, 18) == pos(17)
        @test aln_to_query(c, 19) == pos(18)
        @test aln_to_query(c, 23) == pos(22)

        @test ref_to_query(c, 5) == pos(11)
        @test ref_to_query(c, 6) == gap(11)
        @test ref_to_query(c, 7) == gap(11)
        @test ref_to_query(c, 8) == pos(12)
        @test ref_to_query(c, 9) == pos(14)
        @test ref_to_query(c, 10) == gap(15)
        @test ref_to_query(c, 11) == pos(16)
        @test ref_to_query(c, 13) == pos(22)

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
            @test collect(it) == [i => query_to_aln(c, i) for i in 1:3:25]

            v = [-10, 0, 1, 3, 8, 8, 13, 15, 21]
            for coordinate in [query, ref, aln]
                it = pos_to_pos(coordinate, coordinate, c, v)
                @test length(it) == length(v)
                @test [last(i).pos for i in it] == map(v) do i
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

        @testset "try_parse catches total length overflow" begin
            n = UInt32(268435455)
            push_u32le!(v::Vector{UInt8}, x::UInt32) = append!(v, reinterpret(UInt8, [htol(x)]))

            # Ref length overflows UInt32 while aln length stays tiny.
            ref_bytes = UInt8[]
            op_n = (n << 4) | UInt32(0x03)
            op_m = (UInt32(1) << 4) | UInt32(0x00)
            for _ in 1:17
                push_u32le!(ref_bytes, op_n)
            end
            push_u32le!(ref_bytes, op_m)
            err = CIGARStrings.try_parse(BAMCIGAR, ref_bytes)
            @test err isa CIGARStrings.CIGARError
            @test err.kind == Errors.IntegerOverflow

            # Query length overflows UInt32 while aln length is still <= UInt32 max.
            query_bytes = UInt8[]
            op_s = (n << 4) | UInt32(0x04)
            op_i = (n << 4) | UInt32(0x01)
            op_i16 = (UInt32(16) << 4) | UInt32(0x01)
            push_u32le!(query_bytes, op_s)
            for _ in 1:14
                push_u32le!(query_bytes, op_i)
            end
            push_u32le!(query_bytes, op_i16)
            push_u32le!(query_bytes, op_s)
            err = CIGARStrings.try_parse(BAMCIGAR, query_bytes)
            @test err isa CIGARStrings.CIGARError
            @test err.kind == Errors.IntegerOverflow
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


    @testset "BAMCIGAR actually creates BAMCIGAR" begin
        c = CIGAR("1S5M1D3M2S")
        c2 = BAMCIGAR(c)
        @test c2 isa BAMCIGAR
        @test c2 == c
    end
end

# TODO: More equality tests
# query length
# ref length
# aln length

@testset "is_compatible" begin
    @testset "CIGAR with CIGAR" begin
        # Empty CIGARs
        @test is_compatible(CIGAR(""), CIGAR(""))

        # Identical CIGARs
        @test is_compatible(CIGAR("10M"), CIGAR("10M"))
        @test is_compatible(CIGAR("5M3I2D"), CIGAR("5M3I2D"))

        # Collapsing consecutive operations
        @test is_compatible(CIGAR("1M1M"), CIGAR("2M"))
        @test is_compatible(CIGAR("1D1D1D"), CIGAR("3D"))
        @test is_compatible(CIGAR("5I"), CIGAR("2I3I"))
        @test is_compatible(CIGAR("1M1M1M1M"), CIGAR("2M2M"))

        # OP_M encompasses OP_X and OP_Eq
        @test is_compatible(CIGAR("3M"), CIGAR("1X1M1="))
        @test is_compatible(CIGAR("1D2M3I"), CIGAR("1D1X1M3I"))
        @test is_compatible(CIGAR("5M"), CIGAR("5="))
        @test is_compatible(CIGAR("5M"), CIGAR("5X"))
        @test is_compatible(CIGAR("10M"), CIGAR("3=4X3="))

        # OP_P and OP_H has no semantic meaning and is skipped
        @test is_compatible(CIGAR("1D1P2D1M3P"), CIGAR("3D1M"))
        @test is_compatible(CIGAR("5P"), CIGAR(""))
        @test is_compatible(CIGAR("1P1P1P"), CIGAR(""))
        @test is_compatible(CIGAR("2M1P3M"), CIGAR("2M3M"))
        @test is_compatible(CIGAR("1P2M1P3I1P"), CIGAR("2M3I"))
        @test is_compatible(CIGAR("3H5M"), CIGAR("3=2X"))
        @test is_compatible(CIGAR("3H5M2H"), CIGAR("3=2X"))
        @test is_compatible(CIGAR("5M2H"), CIGAR("2=3X"))

        # Complex combinations
        @test is_compatible(CIGAR("1H1S3M1P2I"), CIGAR("1S1=2M2I"))

        # Incompatible CIGARs
        @test !is_compatible(CIGAR("1="), CIGAR("1X"))
        @test !is_compatible(CIGAR("5M"), CIGAR("4M"))
        @test !is_compatible(CIGAR("5M"), CIGAR("6M"))
        @test !is_compatible(CIGAR("3M2I"), CIGAR("3M2D"))
        @test !is_compatible(CIGAR("1D1M"), CIGAR("1M1D"))
        @test !is_compatible(CIGAR("5="), CIGAR("4=1X"))
        @test !is_compatible(CIGAR("10M"), CIGAR(""))
        @test !is_compatible(CIGAR(""), CIGAR("10M"))
        @test !is_compatible(CIGAR("1S2M3I"), CIGAR("1S2M2I"))
    end

    @testset "BAMCIGAR with BAMCIGAR" begin
        # Empty CIGARs
        @test is_compatible(BAMCIGAR(CIGAR("")), BAMCIGAR(CIGAR("")))

        # Identical CIGARs
        @test is_compatible(BAMCIGAR(CIGAR("10M")), BAMCIGAR(CIGAR("10M")))

        # Collapsing consecutive operations
        @test is_compatible(BAMCIGAR(CIGAR("1M1M")), BAMCIGAR(CIGAR("2M")))
        @test is_compatible(BAMCIGAR(CIGAR("1D1D1D")), BAMCIGAR(CIGAR("3D")))

        # OP_M encompasses OP_X and OP_Eq
        @test is_compatible(BAMCIGAR(CIGAR("3M")), BAMCIGAR(CIGAR("1X1M1=")))
        @test is_compatible(BAMCIGAR(CIGAR("1D2M3I")), BAMCIGAR(CIGAR("1D1X1M3I")))

        # OP_P and OP_H have no semantic meaning
        @test is_compatible(BAMCIGAR(CIGAR("1D1P2D1M3P")), BAMCIGAR(CIGAR("3D1M")))
        @test is_compatible(BAMCIGAR(CIGAR("5P")), BAMCIGAR(CIGAR("")))
        @test is_compatible(BAMCIGAR(CIGAR("3H5M")), BAMCIGAR(CIGAR("3=2X")))
        @test is_compatible(BAMCIGAR(CIGAR("3H5M2H")), BAMCIGAR(CIGAR("3=2X")))

        # Incompatible CIGARs
        @test !is_compatible(BAMCIGAR(CIGAR("1=")), BAMCIGAR(CIGAR("1X")))
        @test !is_compatible(BAMCIGAR(CIGAR("5M")), BAMCIGAR(CIGAR("4M")))
        @test !is_compatible(BAMCIGAR(CIGAR("3M2I")), BAMCIGAR(CIGAR("3M2D")))
    end

    @testset "CIGAR with BAMCIGAR" begin
        # Compatible pairs
        @test is_compatible(CIGAR("10M"), BAMCIGAR(CIGAR("10M")))
        @test is_compatible(BAMCIGAR(CIGAR("10M")), CIGAR("10M"))

        @test is_compatible(CIGAR("1M1M"), BAMCIGAR(CIGAR("2M")))
        @test is_compatible(BAMCIGAR(CIGAR("1M1M")), CIGAR("2M"))

        @test is_compatible(CIGAR("3M"), BAMCIGAR(CIGAR("1X1=1M")))
        @test is_compatible(BAMCIGAR(CIGAR("5M")), CIGAR("2=3X"))

        @test is_compatible(CIGAR("3H5M"), BAMCIGAR(CIGAR("3=2X")))
        @test is_compatible(BAMCIGAR(CIGAR("3H5M2H")), CIGAR("3=2X"))

        @test is_compatible(CIGAR("2M1P3M"), BAMCIGAR(CIGAR("5M")))
        @test is_compatible(BAMCIGAR(CIGAR("1P2D1P")), CIGAR("2D"))

        # Incompatible pairs
        @test !is_compatible(CIGAR("1="), BAMCIGAR(CIGAR("1X")))
        @test !is_compatible(BAMCIGAR(CIGAR("5M")), CIGAR("6M"))
        @test !is_compatible(CIGAR("2M3I"), BAMCIGAR(CIGAR("2M3D")))
    end
end

@testset "Normalization" begin
    make(::Type{CIGAR}, s) = CIGAR(s)
    make(::Type{BAMCIGAR}, s) = BAMCIGAR(CIGAR(s))

    @testset "Consecutive ops collapsed after conversion" for T in [CIGAR, BAMCIGAR]
        # X, = and M all become M and merge together
        @test normalize(make(T, "1X1=3P2M1X")) == make(T, "5M")

        # Eq and X converted to M, merged with existing M
        @test normalize(make(T, "1=1X1M")) == make(T, "3M")

        # H is deleted
        @test normalize(make(T, "3H10M2H")) == make(T, "10M")

        # Already normalized: unchanged
        @test normalize(make(T, "5M1D8M1I7M")) == make(T, "5M1D8M1I7M")

        # All P removed, empty result
        @test normalize(make(T, "1P2P3P")) == make(T, "")

        # Multiple consecutive same-op elements collapsed
        @test normalize(make(T, "1M2M3M")) == make(T, "6M")
        @test normalize(make(T, "1D1D1D")) == make(T, "3D")
        @test normalize(make(T, "2I3I")) == make(T, "5I")

        # P between different ops does not cause merge
        @test normalize(make(T, "3M1P2I")) == make(T, "3M2I")

        # N is not converted
        @test normalize(make(T, "5M3N10M")) == make(T, "5M3N10M")

        # Single element conversions
        @test normalize(make(T, "5H")) == make(T, "")
        @test normalize(make(T, "5=")) == make(T, "5M")
        @test normalize(make(T, "5X")) == make(T, "5M")

        # P at start before M-compatible op
        @test normalize(make(T, "1P3M")) == make(T, "3M")
        @test normalize(make(T, "1P2P3=")) == make(T, "3M")
    end

    @testset "Preserves alignment properties" for T in [CIGAR, BAMCIGAR]
        for s in [
                "1X1=3P2M1X",
                "5H9S1D1D1D2I9S6H",
                "1M2M1=2=3M",
                "3H10M2H",
            ]
            c = make(T, s)
            cn = normalize(c)
            @test aln_length(cn) == aln_length(c)
            @test ref_length(cn) == ref_length(c)
            @test query_length(cn) == query_length(c)
        end
    end

    @testset "Integer overflow on merge" for T in [CIGAR, BAMCIGAR]
        # Merging 268435455 + 1 exceeds max element length (2^28 - 1)
        c = make(T, "268435455M1=")
        err = try
            normalize(c)
        catch e
            e
        end
        @test err isa CIGARStrings.CIGARError
        @test err.kind == Errors.IntegerOverflow

        # Also overflows when both operands are large
        c = make(T, "134217728M134217728=")
        err = try
            normalize(c)
        catch e
            e
        end
        @test err isa CIGARStrings.CIGARError
        @test err.kind == Errors.IntegerOverflow

        # Just at the boundary: should succeed
        c = make(T, "268435454M1=")
        cn = normalize(c)
        @test cn == make(T, "268435455M")
    end

    # normalize! requires no more memory than the final output size.
    @testset "normalize! memory requirement" for T in [CIGAR, BAMCIGAR]
        # (input, expected normalized string, number of ops after normalization)
        test_cases = [
            ("1X1=3P2M1X", "5M", 1),
            ("9M1=", "10M", 1), # merge causes digit growth
            ("1M2M1=2=3M", "9M", 1),
            ("3M1P2I", "3M2I", 2),
        ]
        for (input, expected_str, expected_n_ops) in test_cases
            c = make(T, input)
            needed = T === CIGAR ? sizeof(expected_str) : 4 * expected_n_ops

            # Buffer of exactly the output size works
            mem = MemoryView(zeros(UInt8, needed))
            cn = normalize!(c, mem)
            @test cn == make(T, expected_str)

            # One byte less is too small
            if needed > 0
                mem_small = MemoryView(zeros(UInt8, needed - 1))
                @test_throws LightBoundsError normalize!(c, mem_small)
            end
        end
    end

    @testset "unsafe_normalize" for T in [CIGAR, BAMCIGAR]
        for (input, expected) in [
                ("1X1=3P2M1X", "5M"),
                ("5M1D8M1I7M", "5M1D8M1I7M"),
                ("1P2P3P", ""),
            ]
            c = copy(make(T, input))
            cn = unsafe_normalize(c)
            @test cn == make(T, expected)
            @test cn isa T
        end
    end
end

@testset "encode! and encode_append!" begin
    cigars = [
        CIGAR(""),
        CIGAR("3S5M1D19X2I6M3H"),
        CIGAR("150M"),
        CIGAR("5H9S1D1D1D2I9S6H"),
    ]

    @testset "encode!" begin
        for c in cigars, T in [CIGAR, BAMCIGAR]
            for src in [c, BAMCIGAR(c)]
                mem = MemoryView(zeros(UInt8, 256))
                result = encode!(mem, T, src)
                @test result isa T
                @test result == src
                # Result memory is a prefix slice of mem
                nbytes = length(MemoryView(result))
                if nbytes > 0
                    @test MemoryView(result) === ImmutableMemoryView(mem)[1:nbytes]
                end
            end
        end

        # LightBoundsError with exactly one byte too little
        for s in ["10M5D3I", "150M", "5H9S1D1D1D2I9S6H"]
            c = CIGAR(s)
            bc = BAMCIGAR(c)
            for (T, src) in [(CIGAR, c), (BAMCIGAR, c), (CIGAR, bc), (BAMCIGAR, bc)]
                # First encode with large buffer to find exact size needed
                big = MemoryView(zeros(UInt8, 256))
                needed = length(MemoryView(encode!(big, T, src)))
                needed == 0 && continue
                # One byte too few
                too_small = MemoryView(zeros(UInt8, needed - 1))
                @test_throws LightBoundsError encode!(too_small, T, src)
                # Exact size works
                exact = MemoryView(zeros(UInt8, needed))
                result = encode!(exact, T, src)
                @test result == src
                @test MemoryView(result) === ImmutableMemoryView(exact)
            end
        end
    end

    @testset "encode_append!" begin
        for c in cigars, T in [CIGAR, BAMCIGAR]
            for src in [c, BAMCIGAR(c)]
                prefix = UInt8[0xaa, 0xbb, 0xcc]
                v = copy(prefix)
                result = encode_append!(v, T, src)
                @test result isa T
                @test result == src
                # Prefix bytes are preserved
                @test v[1:3] == prefix
                # Result is backed by the appended portion of v
                if !isempty(MemoryView(result))
                    @test MemoryView(result) === ImmutableMemoryView(v)[4:end]
                end
            end
        end

        # Appending to empty vector
        for (T, src) in [(BAMCIGAR, CIGAR("5M")), (CIGAR, BAMCIGAR(CIGAR("5M")))]
            v = UInt8[]
            result = encode_append!(v, T, src)
            @test result isa T
            @test result == src
            @test MemoryView(result) === ImmutableMemoryView(v)
        end

        # Multiple appends to the same vector
        c1 = CIGAR("10M")
        c2 = CIGAR("3S5M2I")
        v = UInt8[]
        r1 = encode_append!(v, BAMCIGAR, c1)
        r2 = encode_append!(v, BAMCIGAR, c2)
        @test r1 == c1
        @test r2 == c2
        @test length(v) == length(MemoryView(r1)) + length(MemoryView(r2))
    end
end

end # module CIGARTests
