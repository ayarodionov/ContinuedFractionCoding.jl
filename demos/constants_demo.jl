# Three more constants, three more morals.
#
#   1. Champernowne's 0.123456789101112... -- provably normal in base 10,
#      yet its CF hides quotients so large they overflow the decoder's
#      Int64.  Typical in one number system, pathological in another.
#   2. Khinchin's constant K: the geometric mean of the CF quotients of
#      almost every real is K = 2.6854...  We check the theorem on pi's
#      quotients -- and peek at K's own.
#   3. The fine-structure constant alpha ~ 1/137.036: measured, not
#      decreed, so its CF is only meaningful up to the error bar -- a
#      continued fraction with an experimental horizon.
#
# Run with:  julia --project demos/constants_demo.jl
using ContinuedFractionCoding
using Printf

function big_quotients(num::BigInt, den::BigInt, nmax::Int)
    out = BigInt[]
    while num > 0 && length(out) < nmax
        a, r = fldmod(den, num)
        push!(out, a)
        num, den = r, num
    end
    return out
end

shorten(a) = ndigits(a) > 12 ? "<$(ndigits(a))-digit number>" : string(a)

as_text(q) = join(map(a -> 32 <= a - 1 <= 126 ? Char(a - 1) : '.', q))

# ---- 1. Champernowne: normal in base 10, monstrous in CF ----
println("1. Champernowne's constant 0.12345678910111213...")
function champernowne_digits(n)
    io = IOBuffer()
    k = 1
    while io.size < n
        print(io, k)
        k += 1
    end
    return String(take!(io))[1:n]
end
m = parse(BigInt, champernowne_digits(600))
q = big_quotients(m, big(10)^600, 24)
println("   CF: [0; ", join(map(shorten, q[1:20]), ", "), ", ...]")
bits = cf_encode(q)
@printf("   %d quotients encode to %d bits -- the monster alone is worth ~%d\n",
        length(q), length(bits), 2 * ndigits(maximum(q)) * 4)
print("   decoding back: ")
try
    cf_decode(bits)
catch e
    println(typeof(e), " -- the ", ndigits(maximum(q)),
            "-digit quotient overflows the decoder's Int64.")
end
println("   message : \"", as_text(q), "\"")
println("   Base 10 sees perfectly uniform digits; CF sees structure so")
println("   extreme it breaks the coder.\n")

# ---- 2. Khinchin: almost every real's quotients have geometric mean K ----
println("2. Khinchin's constant K = 2.6854520010...")
j = 20_000
setprecision(BigFloat, j + 128)
mpi = BigInt(floor((BigFloat(pi) - 3) * (big(1) << j)))
qpi = big_quotients(mpi, big(1) << j, 5_000)
for n in (10, 100, 1000, 5000)
    gm = exp(sum(log, Float64.(qpi[1:n])) / n)
    @printf("   geometric mean of pi's first %4d quotients: %.4f\n", n, gm)
end
println("   -> crawling toward K = 2.6854, as Khinchin's theorem demands.")
println("   pi's quotients as text: \"", as_text(qpi[1:100]), "\"")
qk = big_quotients(big(26854520010653064), big(10)^16, 12)
println("   K's own CF starts [2; ", join(qk[2:end], ", "), ", ...] --")
println("   as text \"", as_text(qk), "\" -- whether K obeys its own")
println("   theorem is an open problem.\n")

# ---- 3. Fine-structure constant: a CF with an error bar ----
println("3. Fine-structure constant, 1/alpha = 137.035999177(21)  (CODATA 2022)")
den0 = big(137035999177)                       # 1/alpha * 10^9
lo  = big_quotients(big(10)^9, den0 + 21, 20)  # alpha at either end
mid = big_quotients(big(10)^9, den0,      20)  #   of the measured
hi  = big_quotients(big(10)^9, den0 - 21, 20)  #   uncertainty interval
horizon = findfirst(i -> !(lo[i] == mid[i] == hi[i]), 1:20) - 1
println("   central : [0; ", join(mid[1:horizon+2], ", "), ", ...]")
println("   the trustworthy part as text: \"", as_text(mid[1:horizon]), "\"")
println("   the three CFs (value, value+sigma, value-sigma) agree on ",
        horizon, " quotients,")
println("   then dissolve into the error bar.  h is exact by decree, pi by")
println("   mathematics -- alpha's continued fraction simply ends at the")
println("   precision of experiment: physics as a ", horizon, "-symbol message.")
