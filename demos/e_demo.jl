# Fun: decode the "message" in the continued fraction of e -- the
# counterpart to pi_demo.jl.  Where pi's quotients look like Gauss-Kuzmin
# noise, e's follow Euler's exact pattern [2; 1, 2, 1, 1, 4, 1, 1, 6, ...]:
# e is as atypical a real as they come, and the coder's rate shows it.
#
# Run with:  julia --project demos/e_demo.jl
using ContinuedFractionCoding
using Printf

j = 50_000                                    # binary digits of e - 2 used
setprecision(BigFloat, j + 128)
m = BigInt(floor((BigFloat(MathConstants.e) - 2) * (big(1) << j)))  # m/2^j ~ e - 2

# ---- 1. decode e's bits as a CF code stream ----
n = 30
bits = BitVector()
ContinuedFractionCoding.push_gamma!(bits, n + 1)
for i in j-1:-1:0
    push!(bits, (m >> i) & 1 == 1)
end
q = cf_decode(bits)
println("e = [2; ", join(q, ", "), ", ...]")

v = 2 + cf_value(q)
@printf("value of that fraction: %.20f\n", Float64(v))
@printf("e                     : %.20f\n\n", Float64(BigFloat(MathConstants.e)))

# the "message" e spells, one byte per quotient (a - 1), printables only
nq = 5_000
function gauss_quotients(num::BigInt, den::BigInt, nq::Int)
    out = Vector{Int}(undef, nq)
    for i in 1:nq
        a, r = fldmod(den, num)
        out[i] = Int(a)
        num, den = r, num
    end
    return out
end
quot = gauss_quotients(m, big(1) << j, nq)
@assert quot[1:n] == q                        # decoder agrees with Gauss map
@assert all(k -> quot[3k-1] == 2k, 1:100)     # Euler's pattern: 2k at 3k-1

text = join(map(a -> 32 <= a - 1 <= 126 ? Char(a - 1) : '.', quot[1:150]))
println("as a byte message: \"", text, "\"")
println("-- once 2k - 1 > 31, e politely recites the odd ASCII characters\n")

# ---- 2. e is NOT a Gauss-Kuzmin source: the rate gives it away ----
code = cf_encode(quot)
@printf("%d quotients of e re-encoded : %.4f bits/quotient\n",
        nq, length(code) / nq)
@printf("Gauss-Kuzmin entropy rate    : %.4f bits/quotient (pi got 3.4124)\n",
        GK_ENTROPY)
@printf("true entropy of the pattern  : 0 bits/quotient -- it's deterministic\n")
@printf("largest quotient seen        : %d (at position %d)\n",
        maximum(quot), argmax(quot))
