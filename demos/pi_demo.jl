# Fun: decode the "message" hidden in the continued fraction of pi.
#
# The binary digits of pi - 3 are, to our decoder, a perfectly valid code
# stream: prepend a length header and cf_decode extracts the partial
# quotients [3; 7, 15, 1, 292, ...].  And since pi is (conjecturally) a
# typical real, its quotients follow the Gauss-Kuzmin law -- so re-encoding
# them should cost almost exactly GK_ENTROPY bits per quotient.
#
# Run with:  julia --project demos/pi_demo.jl
using ContinuedFractionCoding
using Printf

j = 20_000                                    # binary digits of pi - 3 used
setprecision(BigFloat, j + 128)
m = BigInt(floor((BigFloat(pi) - 3) * (big(1) << j)))   # m/2^j ~ pi - 3

# ---- 1. decode pi's bits as a CF code stream ----
n = 30
bits = BitVector()
ContinuedFractionCoding.push_gamma!(bits, n + 1)
for i in j-1:-1:0
    push!(bits, (m >> i) & 1 == 1)
end
q = cf_decode(bits)
println("pi = [3; ", join(q, ", "), ", ...]")

v = 3 + cf_value(q)
@printf("value of that fraction: %.20f\n", Float64(v))
@printf("pi                    : %.20f\n\n", Float64(BigFloat(pi)))

# the "message" pi spells, one byte per quotient (a - 1), printables only
text = join(map(a -> 32 <= a - 1 <= 126 ? Char(a - 1) : '.', q))
println("as a byte message: \"", text, "\"  -- alas, pi keeps its secrets\n")

# ---- 2. is pi a Gauss-Kuzmin source?  re-encode its quotients ----
function gauss_quotients(num::BigInt, den::BigInt, nq::Int)
    out = Vector{Int}(undef, nq)
    for i in 1:nq
        a, r = fldmod(den, num)
        out[i] = Int(a)
        num, den = r, num
    end
    return out
end
nq = 5_000
quot = gauss_quotients(m, big(1) << j, nq)
@assert quot[1:n] == q                        # decoder agrees with Gauss map

code = cf_encode(quot)
@printf("%d quotients of pi re-encoded: %.4f bits/quotient\n",
        nq, length(code) / nq)
@printf("Gauss-Kuzmin entropy rate    : %.4f bits/quotient\n", GK_ENTROPY)
@printf("largest quotient seen        : %d (at position %d)\n",
        maximum(quot), argmax(quot))
