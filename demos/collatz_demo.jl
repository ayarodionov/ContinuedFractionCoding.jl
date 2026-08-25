# Fun: the 3n+1 (Collatz) trajectory of a number, read as the partial
# quotients of a continued fraction and decoded as a byte message.
#
# The trajectory of n is a perfectly valid message (all values >= 1), so
# [0; n, c1, c2, ..., 1] is a real number that *is* the trajectory.  Bytes
# are quotient - 1, printables shown, the rest as '.' -- whenever the orbit
# passes through 34..127 it spells characters.
#
# Run with:  julia --project demos/collatz_demo.jl
using ContinuedFractionCoding
using Printf

function collatz(n::Int)
    seq = [n]
    while n != 1
        n = iseven(n) ? n >> 1 : 3n + 1
        push!(seq, n)
    end
    return seq
end

as_text(q) = join(map(a -> 32 <= a - 1 <= 126 ? Char(a - 1) : '.', q))

for n0 in (27, 75, 88)
    seq = collatz(n0)
    bits = cf_encode(seq)
    @assert cf_decode(bits) == seq
    x = cf_value(seq)

    println("start ", n0, ": trajectory of ", length(seq),
            " steps, peak ", maximum(seq))
    head = seq[1:min(end, 12)]
    println("  CF      : [0; ", join(head, ", "),
            length(seq) > 12 ? ", ...]" : "]")
    @printf("  value   : %.15f...  (a rational with %d-digit denominator)\n",
            Float64(x), ndigits(denominator(x)))
    println("  message : \"", as_text(seq), "\"")
    @printf("  code    : %d bits = %.2f bits/step -- vs log2(%d) = %.1f bits\n",
            length(bits), length(bits) / length(seq), n0, log2(n0))
    println("            (the whole trajectory is determined by the start ",
            "number: true entropy ~", round(Int, log2(n0)), " bits)")
    println()
end
