# Fun: the continued fraction of Planck's constant.
#
# Since the 2019 SI redefinition, h = 6.62607015e-34 J*s *exactly* -- the
# constant is a rational number by decree, so unlike pi and e its continued
# fraction terminates: the quantum of action is a finite message.  The
# reduced constant hbar = h/2pi inherits pi's irrationality and babbles
# Gauss-Kuzmin noise forever.
#
# We code the mantissa's fractional part (the number's scale is bookkeeping;
# its digits are the physics).
#
# Run with:  julia --project demos/planck_demo.jl
using ContinuedFractionCoding
using Printf

as_text(q) = join(map(a -> 32 <= a - 1 <= 126 ? Char(a - 1) : '.', q))

function gauss_quotients(num::BigInt, den::BigInt, nmax::Int)
    out = Int[]
    while num > 0 && length(out) < nmax
        a, r = fldmod(den, num)
        push!(out, Int(a))
        num, den = r, num
    end
    return out
end

# ---- h: exact rational since 2019 ----
x = 62607015 // big(10)^8            # h = 6.62607015e-34: mantissa - 6
q = gauss_quotients(numerator(x), denominator(x), 10_000)
println("h = 6.62607015e-34 J*s exactly (SI, 2019)")
println("  0.62607015 = ", cf_string(q))
@assert cf_value(q) == x                       # the CF really terminates
bits = cf_encode(q)
@assert cf_decode(bits) == q
println("  a FINITE fraction: ", length(q), " quotients, ",
        length(bits), " bits of code")
println("  message : \"", as_text(q), "\"\n")

# ---- hbar = h/2pi: pi makes it irrational again ----
j = 20_000
setprecision(BigFloat, j + 128)
hbar = BigFloat(662607015 // big(10)^8) / (2 * BigFloat(pi))  # * 1e-34
frac = hbar - 1                                 # 1.054571817... - 1
m = BigInt(floor(frac * (big(1) << j)))
qh = gauss_quotients(m, big(1) << j, 5_000)
println("hbar = h/2pi = 1.054571817...e-34 J*s -- irrational, thanks to pi")
println("  0.0545718... = [0; ", join(qh[1:14], ", "), ", ...]")
println("  message : \"", as_text(qh[1:120]), "\"")
code = cf_encode(qh)
@printf("  %d quotients re-encoded : %.4f bits/quotient\n",
        length(qh), length(code) / length(qh))
@printf("  Gauss-Kuzmin entropy    : %.4f bits/quotient\n", GK_ENTROPY)
println("\nBy decree h is a finite message; hbar talks forever and says nothing.")
