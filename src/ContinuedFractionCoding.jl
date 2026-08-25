"""
    ContinuedFractionCoding

Continued-fraction coding, just for fun.

A message of positive integers `a1, ..., an` is read as the partial
quotients of the continued fraction `x = [0; a1, a2, ..., an]`.  The set of
reals whose expansion starts with these quotients is an interval (a cylinder
of the Gauss map) with endpoints `h_n/k_n` and `(h_n + h_{n-1})/(k_n + k_{n-1})`
built from the convergents.  The encoder emits the binary expansion of a
shortest dyadic rational strictly inside the open cylinder; the decoder runs
the Gauss map `x -> 1/x - floor(1/x)` on that rational, in exact integer
arithmetic, and reads the quotients back.

This is the zero-redundancy arithmetic code for a memoryless source with the
Gauss-Kuzmin distribution `P(a = k) = log2(1 + 1/(k(k+2)))`, spending on
average `pi^2 / (6 ln^2 2) ~ 3.4237` bits per symbol (`GK_ENTROPY`).  For any
other source it is a fixed code, not an adaptive one.

The stream is self-delimiting: an Elias-gamma header carries the message
length, followed by the dyadic point's bits (trailing zeros dropped).
"""
module ContinuedFractionCoding

using Random

export cf_encode, cf_decode, cf_decode_string
export cf_show, cf_string, cf_value, cf_interval, cf_point
export gauss_kuzmin_prob, rand_gauss_kuzmin, GK_ENTROPY

"Entropy rate of the CF digit source, `pi^2 / (6 ln^2 2)` bits per symbol."
const GK_ENTROPY = pi^2 / (6 * log(2)^2)

"Gauss-Kuzmin probability of the partial quotient `k >= 1`."
gauss_kuzmin_prob(k::Integer) = log2(1 + 1 / (k * (k + 2)))

"""
    rand_gauss_kuzmin([rng]) -> Int

Sample a partial quotient from the Gauss-Kuzmin distribution by inverting
its CDF `P(a <= k) = 1 + log2((k+1)/(k+2))`.
"""
function rand_gauss_kuzmin(rng::AbstractRNG = Random.default_rng())
    t = exp2(rand(rng))                     # in [1, 2)
    return max(ceil(Int, 2 * (t - 1) / (2 - t)), 1)
end

# ---------------------------------------------------------------- encoding

"""
    cf_encode(msg) -> BitVector

Encode a vector of integers `>= 1`.  Output is an Elias-gamma header with
the message length followed by the binary expansion of a dyadic point
inside the cylinder `[0; msg...]`.
"""
function cf_encode(msg::AbstractVector{<:Integer})
    m, j = dyadic_point(msg)
    out = BitVector()
    push_gamma!(out, length(msg) + 1)
    for i in j-1:-1:0
        push!(out, (m >> i) & 1 == 1)
    end
    return out
end

# convergents h_n/k_n and h_{n-1}/k_{n-1} of [0; a1, ..., an]
function convergents(msg::AbstractVector{<:Integer})
    all(a -> a >= 1, msg) ||
        throw(ArgumentError("symbols must be integers >= 1"))
    hprev, kprev = big(1), big(0)
    h, k = big(0), big(1)
    for a in msg
        h, hprev = big(a) * h + hprev, h
        k, kprev = big(a) * k + kprev, k
    end
    return h, k, hprev, kprev
end

# the point m/2^j the encoder emits for this message
function dyadic_point(msg::AbstractVector{<:Integer})
    h, k, hprev, kprev = convergents(msg)

    # open cylinder interval: between h/k and the mediant (h+hprev)/(k+kprev)
    lon, lod = h, k
    hin, hid = h + hprev, k + kprev
    if lon * hid > hin * lod
        (lon, lod), (hin, hid) = (hin, hid), (lon, lod)
    end

    # smallest j with 2^-j < hi - lo, so that floor(lo*2^j) + 1 lands inside
    gapn = hin * lod - lon * hid
    gapd = hid * lod
    j = max(ndigits(gapd; base = 2) - ndigits(gapn; base = 2), 0)
    while gapn << j <= gapd
        j += 1
    end
    m = fld(lon << j, lod) + 1              # lo < m/2^j <= lo + 2^-j < hi
    while iseven(m) && j > 0                # drop trailing zeros: shortest
        m >>= 1                             # dyadic form of this point
        j -= 1
    end
    return m, j
end

"""
    cf_encode(s::AbstractString) -> BitVector

Encode a string: each UTF-8 code unit `b` becomes the partial quotient
`b + 1` (so the quotients are `1 .. 256`). Decode with
[`cf_decode_string`](@ref).
"""
cf_encode(s::AbstractString) = cf_encode(Int.(codeunits(s)) .+ 1)

"""
    cf_decode_string(bits) -> String

Inverse of `cf_encode(::AbstractString)`: decodes the quotients and maps
each back to the UTF-8 code unit `a - 1`.
"""
function cf_decode_string(bits::AbstractVector{Bool})
    q = cf_decode(bits)
    all(a -> 1 <= a <= 256, q) ||
        throw(ArgumentError("quotient out of byte range: not a string stream"))
    return String(UInt8.(q .- 1))
end

# -------------------------------------------------------------- inspection

"""
    cf_string(msg) -> String

The message written in continued-fraction notation, `"[0; a1, a2, ...]"`.
"""
cf_string(msg::AbstractVector{<:Integer}) =
    isempty(msg) ? "[0]" : "[0; " * join(msg, ", ") * "]"
cf_string(s::AbstractString) = cf_string(Int.(codeunits(s)) .+ 1)

"""
    cf_value(msg) -> Rational{BigInt}

The exact value of the finite continued fraction `[0; a1, ..., an]`.
"""
function cf_value(msg::AbstractVector{<:Integer})
    h, k, _, _ = convergents(msg)
    return h // k
end
cf_value(s::AbstractString) = cf_value(Int.(codeunits(s)) .+ 1)

"""
    cf_interval(msg) -> (lo, hi)

The open cylinder interval of reals whose continued fraction starts with
the message: between `h_n/k_n` and the mediant
`(h_n + h_{n-1})/(k_n + k_{n-1})`, as exact `Rational{BigInt}`s.
"""
function cf_interval(msg::AbstractVector{<:Integer})
    h, k, hprev, kprev = convergents(msg)
    lo, hi = minmax(h // k, (h + hprev) // (k + kprev))
    return lo, hi
end
cf_interval(s::AbstractString) = cf_interval(Int.(codeunits(s)) .+ 1)

"""
    cf_point(msg) -> Rational{BigInt}

The dyadic rational inside the cylinder whose binary expansion the encoder
emits; its numerator's bits (to denominator width) are the code body.
"""
function cf_point(msg::AbstractVector{<:Integer})
    m, j = dyadic_point(msg)
    return m // (big(1) << j)
end
cf_point(s::AbstractString) = cf_point(Int.(codeunits(s)) .+ 1)

"""
    cf_show([io], msg)

Print how coding sees the message (a vector of integers `>= 1` or a
string): its continued-fraction notation, exact value, cylinder interval,
and the emitted dyadic point with its code bits.
"""
cf_show(msg) = cf_show(stdout, msg)
cf_show(io::IO, s::AbstractString) = cf_show(io, Int.(codeunits(s)) .+ 1)
function cf_show(io::IO, msg::AbstractVector{<:Integer})
    lo, hi = cf_interval(msg)
    m, j = dyadic_point(msg)
    bits = reverse(digits(m; base = 2, pad = j))
    println(io, "message  : ", cf_string(msg))
    println(io, "value    : ", cf_value(msg), " = ", Float64(cf_value(msg)))
    println(io, "cylinder : (", lo, ", ", hi, "), width ", Float64(hi - lo))
    println(io, "point    : ", m, "/2^", j, " = 0.", join(bits), "_2")
    println(io, "code     : ", join(Int.(cf_encode(msg))),
            "  (", length(cf_encode(msg)), " bits: gamma header + point)")
    return nothing
end

# ---------------------------------------------------------------- decoding

"""
    cf_decode(bits) -> Vector{Int}

Inverse of [`cf_encode`](@ref): reads the length header, reconstructs the
dyadic rational and extracts its first `n` partial quotients with the Gauss
map.
"""
function cf_decode(bits::AbstractVector{Bool})
    n, pos = read_gamma(bits)
    n -= 1

    j = length(bits) - pos + 1              # x = m / 2^j from the tail bits
    num = big(0)
    for i in pos:lastindex(bits)
        num = num << 1 | (bits[i] ? 1 : 0)
    end
    den = big(1) << j

    out = Vector{Int}(undef, n)
    for i in 1:n
        num > 0 && num < den ||
            throw(ArgumentError("corrupt or truncated stream"))
        a, r = fldmod(den, num)             # 1/x = a + r/num
        r > 0 || i == n ||
            throw(ArgumentError("corrupt or truncated stream"))
        out[i] = Int(a)
        num, den = r, num                   # Gauss map: x <- 1/x - a
    end
    return out
end

# ------------------------------------------------------------- elias gamma

function push_gamma!(out::BitVector, N::Integer)
    L = ndigits(N; base = 2)
    for _ in 1:L-1
        push!(out, false)
    end
    for i in L-1:-1:0
        push!(out, (N >> i) & 1 == 1)
    end
    return out
end

function read_gamma(bits::AbstractVector{Bool})
    pos = firstindex(bits)
    z = 0
    while pos <= lastindex(bits) && !bits[pos]
        z += 1
        pos += 1
    end
    pos + z <= lastindex(bits) ||
        throw(ArgumentError("corrupt or truncated stream"))
    N = 0
    for i in pos:pos+z
        N = N << 1 | (bits[i] ? 1 : 0)
    end
    return N, pos + z + 1
end

end # module
