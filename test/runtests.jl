using ContinuedFractionCoding
using Test, Random

@testset "ContinuedFractionCoding" begin

    @testset "gauss-kuzmin distribution" begin
        # partial sums telescope to the exact CDF 1 + log2((K+1)/(K+2))
        K = 10^6
        @test sum(gauss_kuzmin_prob, 1:K) ≈ 1 + log2((K + 1) / (K + 2)) atol = 1e-9
        @test gauss_kuzmin_prob(1) ≈ log2(4 / 3)
        @test GK_ENTROPY ≈ 3.4237 atol = 1e-4

        rng = MersenneTwister(1)
        s = [rand_gauss_kuzmin(rng) for _ in 1:200_000]
        for k in 1:5
            @test count(==(k), s) / length(s) ≈ gauss_kuzmin_prob(k) atol = 5e-3
        end
    end

    @testset "roundtrip: edge cases" begin
        @test cf_decode(cf_encode(Int[])) == Int[]
        @test cf_decode(cf_encode([1])) == [1]
        @test cf_decode(cf_encode([7])) == [7]
        @test cf_decode(cf_encode([1_000_000, 1, 2])) == [1_000_000, 1, 2]
        @test cf_decode(cf_encode(ones(Int, 500))) == ones(Int, 500)   # golden ratio
        @test cf_decode(cf_encode(fill(2, 300))) == fill(2, 300)       # sqrt(2)-ish
    end

    @testset "roundtrip: random messages" begin
        rng = MersenneTwister(42)
        for n in (1, 2, 10, 100, 1000)
            msg = [rand_gauss_kuzmin(rng) for _ in 1:n]
            @test cf_decode(cf_encode(msg)) == msg
            msg = rand(rng, 1:50, n)                    # non-GK source
            @test cf_decode(cf_encode(msg)) == msg
        end
    end

    @testset "rate on the gauss-kuzmin source" begin
        rng = MersenneTwister(7)
        n = 5_000
        msg = [rand_gauss_kuzmin(rng) for _ in 1:n]
        rate = length(cf_encode(msg)) / n
        @test GK_ENTROPY - 0.2 < rate < GK_ENTROPY + 0.2
    end

    @testset "strings" begin
        for s in ("", "a", "Hello, p-adic world!",
                  "Привет, мир! 🎲 ∑ garçon", "\0\xff" |> String)
            @test cf_decode_string(cf_encode(s)) == s
        end
        # a stream with quotients > 256 is rejected as a string
        @test_throws ArgumentError cf_decode_string(cf_encode([300, 1, 2]))
    end

    @testset "inspection" begin
        @test cf_string([3, 1, 2]) == "[0; 3, 1, 2]"
        @test cf_string(Int[]) == "[0]"
        @test cf_value([2]) == 1 // 2
        @test cf_value([3, 1, 2]) == 3 // 11          # 1/(3 + 1/(1 + 1/2))
        @test cf_value([1, 1, 1, 1, 1]) == 5 // 8     # Fibonacci convergent
        @test cf_interval(Int[]) == (0 // 1, 1 // 1)
        lo, hi = cf_interval([3, 1, 2])
        @test lo == 4 // 15 && hi == 3 // 11          # mediant (3+1)/(11+4)
        p = cf_point([3, 1, 2])
        @test lo < p < hi
        @test ispow2(denominator(p))
        out = sprint(cf_show, [3, 1, 2])
        @test occursin("[0; 3, 1, 2]", out) && occursin("3//11", out)
        @test cf_string("a") == "[0; 98]"
        @test cf_value("a") == 1 // 98
    end

    @testset "input validation" begin
        @test_throws ArgumentError cf_encode([1, 0, 2])
        @test_throws ArgumentError cf_encode([-3])
        # truncation moves the point out of the cylinder: the decoder either
        # detects it or returns wrong quotients, but never the original
        msg = [5, 4, 3, 2]
        bits = cf_encode(msg)
        @test try
            cf_decode(bits[1:end-3]) != msg
        catch e
            e isa ArgumentError
        end
        @test_throws ArgumentError cf_decode(falses(4))      # no gamma terminator
    end

end
