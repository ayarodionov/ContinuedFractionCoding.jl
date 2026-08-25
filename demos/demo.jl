# Demo: continued-fraction coding of a Gauss-Kuzmin source (its natural
# source) versus mismatched sources.  Run with:  julia --project demos/demo.jl
using ContinuedFractionCoding
using Random, Printf

rng = MersenneTwister(2026)
n = 20_000

# 1. The matched source: symbols ~ Gauss-Kuzmin.  CF coding is the
#    zero-redundancy arithmetic code for exactly this distribution.
msg = [rand_gauss_kuzmin(rng) for _ in 1:n]
bits = cf_encode(msg)
@assert cf_decode(bits) == msg
H = GK_ENTROPY
@printf("Gauss-Kuzmin source, n = %d\n", n)
@printf("  entropy rate : %.4f bits/sym\n", H)
@printf("  CF code rate : %.4f bits/sym  (overhead %.3f%%)\n\n",
        length(bits) / n, 100 * (length(bits) / n - H) / H)

# 2. Strings: code units become quotients b+1.  A quotient of size k costs
#    about 2*log2(k) bits, so ASCII text (~bytes 32..122) is expensive --
#    CF coding is spectacularly mismatched to text, which is the point.
s = "In the fields of Provence, a mathematician counts in trits."
sb = cf_encode(s)
@assert cf_decode_string(sb) == s
@printf("string demo: %d chars -> %d bits, %.2f bits/char (vs 8 raw)\n\n",
        ncodeunits(s), length(sb), length(sb) / ncodeunits(s))

# 3. Mismatched sources: CF coding is a fixed code, so it loses here.
for (name, gen, Hsrc) in (
        ("uniform 1..8", () -> rand(rng, 1:8), 3.0),
        ("all ones (golden ratio)", () -> 1, 0.0),
    )
    m = [gen() for _ in 1:n]
    b = cf_encode(m)
    @assert cf_decode(b) == m
    @printf("%-25s entropy %.3f, CF rate %.4f bits/sym\n",
            name, Hsrc, length(b) / n)
end
