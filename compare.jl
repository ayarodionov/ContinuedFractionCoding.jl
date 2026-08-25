# Head-to-head: CF coding (fixed code, Gauss-Kuzmin-optimal) vs p-adic
# arithmetic coding (model-based, entropy-optimal for any static model).
using PadicCoding, ContinuedFractionCoding
using Random, Printf

rng = MersenneTwister(2026)
n = 20_000

entropy(freqs) = begin
    p = freqs ./ sum(freqs)
    -sum(x > 0 ? x * log2(x) : 0.0 for x in p)
end

function padic_rate(msg, nsyms; P = 2, N = 60)
    cnt = zeros(Int, nsyms)
    for s in msg; cnt[s] += 1; end
    m = PadicModel(max.(cnt, 1); P = P, N = N)
    bits = padic_encode(m, msg)
    @assert padic_decode(m, bits) == msg
    return length(bits) * log2(P) / length(msg), entropy(filter(>(0), cnt))
end

println("source                      entropy   CF rate   p-adic rate   (bits/sym)")
println("-"^74)

# 1. Gauss-Kuzmin source: CF's home turf.  For the p-adic coder we cap the
#    (unbounded) alphabet at the observed maximum and model observed counts.
msg = [rand_gauss_kuzmin(rng) for _ in 1:n]
cfr = length(cf_encode(msg)) / n
pr, He = padic_rate(msg, maximum(msg))
@printf("%-27s %7.4f  %8.4f  %10.4f\n", "Gauss-Kuzmin (unbounded)", He, cfr, pr)

# 2. Skewed 16-symbol source: p-adic's home turf (matched static model).
freqs = [max(1, round(Int, 10_000 * 0.7^k)) for k in 0:15]
cum = cumsum(freqs ./ sum(freqs))
msg = [searchsortedfirst(cum, rand(rng)) for _ in 1:n]
cfr = length(cf_encode(msg)) / n
pr, He = padic_rate(msg, 16)
@printf("%-27s %7.4f  %8.4f  %10.4f\n", "skewed 16-symbol", He, cfr, pr)

# 3. Uniform 1..8
msg = rand(rng, 1:8, n)
cfr = length(cf_encode(msg)) / n
pr, He = padic_rate(msg, 8)
@printf("%-27s %7.4f  %8.4f  %10.4f\n", "uniform 1..8", He, cfr, pr)

# 4. Text bytes (README of PadicCoding as sample text), symbols = byte + 1
txt = read(joinpath(dirname(dirname(pathof(PadicCoding))), "README.md"))
msg = Int.(txt) .+ 1
cfr = length(cf_encode(msg)) / length(msg)
pr, He = padic_rate(msg, 257)
@printf("%-27s %7.4f  %8.4f  %10.4f\n", "text bytes (README, $(length(msg)))", He, cfr, pr)

# 5. Quotients of pi (a GK-typical stream)
j = 20_000
setprecision(BigFloat, j + 128)
mpi = BigInt(floor((BigFloat(pi) - 3) * (big(1) << j)))
num, den = mpi, big(1) << j
quot = Vector{Int}(undef, 5000)
for i in eachindex(quot)
    a, r = fldmod(den, num)
    quot[i] = Int(a)
    global num, den = r, num
end
cfr = length(cf_encode(quot)) / length(quot)
pr, He = padic_rate(quot, maximum(quot))
@printf("%-27s %7.4f  %8.4f  %10.4f\n", "quotients of pi (5000)", He, cfr, pr)
