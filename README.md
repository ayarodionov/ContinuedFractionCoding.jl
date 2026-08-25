# ContinuedFractionCoding.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-9558B2.svg)](https://julialang.org)

Continued-fraction coding — a just-for-fun companion to
[PadicCoding.jl](https://github.com/ayarodionov/PadicCoding.jl).

Arithmetic coding needs a number system whose representation prefixes
correspond to nested intervals. Base-P expansions have this property, but so
do continued fractions: the reals whose expansion starts with quotients
`a1, ..., an` form an interval — a cylinder of the Gauss map
`x -> 1/x - floor(1/x)`. This package codes in that number system:

- **Encoder**: reads the message symbols (integers `>= 1`) as the partial
  quotients of `x = [0; a1, ..., an]`, computes the cylinder interval from
  the convergents `h_n/k_n` and the mediant
  `(h_n + h_{n-1})/(k_n + k_{n-1})`, and emits the binary expansion of a
  shortest dyadic rational strictly inside it (the CF analogue of
  `selectPoint`), after an Elias-gamma length header.
- **Decoder**: reconstructs the dyadic rational and extracts the quotients
  by iterating the Gauss map in exact `BigInt` arithmetic.

## Why this is (only) a fixed code

By the Lévy and Rokhlin/Gauss–Kuzmin theorems, a CF digit of a
Lebesgue-random real both carries and costs `pi^2 / (6 ln^2 2) ~ 3.4237`
bits on average (`GK_ENTROPY`). CF coding is therefore the zero-redundancy
arithmetic code for exactly one source: memoryless positive integers with
the Gauss–Kuzmin distribution `P(a = k) = log2(1 + 1/(k(k+2)))` — a
Zipf-like law with a `1/k^2` tail. For that source it has no per-symbol
overhead, the way base-2 arithmetic coding is exact for dyadic models. For
any other source it is a fixed (non-adaptive) code. `demo.jl` shows both
regimes:

```
Gauss-Kuzmin source, n = 20000
  entropy rate : 3.4237 bits/sym
  CF code rate : 3.4356 bits/sym  (overhead 0.346%)

uniform 1..8              entropy 3.000, CF rate 4.0881 bits/sym
all ones (golden ratio)   entropy 0.000, CF rate 1.3900 bits/sym
```

(The all-ones rate is `2 log2(phi) ~ 1.3885` bits/symbol — Lévy's constant
for the slowest-converging continued fraction.)

## Usage

```julia
using ContinuedFractionCoding

msg  = [rand_gauss_kuzmin() for _ in 1:1000]   # its natural source
bits = cf_encode(msg)                          # BitVector
@assert cf_decode(bits) == msg
```

Strings work too — each UTF-8 code unit `b` becomes the quotient `b + 1`:

```julia
bits = cf_encode("garçon")            # BitVector
@assert cf_decode_string(bits) == "garçon"
```

(A quotient of size `k` costs about `2 log2(k)` bits, so ASCII text runs
at ~13 bits/char — CF coding is spectacularly mismatched to text, which
is part of the fun.)

To see the continued fraction a message becomes, use `cf_show` (or the
pieces: `cf_string`, `cf_value`, `cf_interval`, `cf_point`):

```julia
julia> cf_show([3, 1, 2])
message  : [0; 3, 1, 2]
value    : 3//11 = 0.2727272727272727
cylinder : (4//15, 3//11), width 0.006060606060606061
point    : 69/2^8 = 0.01000101_2
code     : 0010001000101  (13 bits: gamma header + point)
```

Exports: `cf_encode`, `cf_decode`, `cf_decode_string`, `cf_show`,
`cf_string`, `cf_value`, `cf_interval`, `cf_point`, `gauss_kuzmin_prob`,
`rand_gauss_kuzmin`, `GK_ENTROPY`.

## Decoding π and e

To a CF decoder, the binary digits of a real number *are* a code stream:
prepend a length header and `cf_decode` extracts the partial quotients.
`pi_demo.jl` does this to 20,000 binary digits of π − 3:

```
pi = [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84, ...]

as a byte message: "....................S........."  -- alas, pi keeps its secrets

5000 quotients of pi re-encoded: 3.4124 bits/quotient
Gauss-Kuzmin entropy rate    : 3.4237 bits/quotient
largest quotient seen        : 20776 (at position 431)
```

π is (conjecturally) a typical real: its quotients are statistically
Gauss–Kuzmin, so the coder finds them incompressible — the re-encoding
rate lands within 0.3% of `GK_ENTROPY`. Read as bytes, the message is
noise; only the quotient 84 surfaces as the letter "S".

`e_demo.jl` is the mirror image. By Euler's 1737 result
`e = [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, ...]` — the quotient `2k` at
every position `3k − 1` — so e is as *atypical* as a real can be:

```
e = [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, ...]

as a byte message: "...............!..#..%..'..)..+..-../..1..3..5..7..9..;..=..?..A..C..E..G..."
-- once 2k - 1 > 31, e politely recites the odd ASCII characters

5000 quotients of e re-encoded : 7.5180 bits/quotient
Gauss-Kuzmin entropy rate      : 3.4237 bits/quotient (pi got 3.4124)
true entropy of the pattern    : 0 bits/quotient -- it's deterministic
```

The growing quotients cost about `2 log2(2k)` bits each, so e re-encodes
at more than double the typical rate — while the true entropy of its
digit stream is zero. e is simultaneously the most expensive number for a
CF coder and infinitely compressible for anyone who knows the formula.
(It also needs 50,000 binary digits for 5,000 quotients where π needs
~17,000: Lévy's constant only holds for typical reals.)

`collatz_demo.jl` plays the same game with 3n+1 trajectories: the orbit
of `n` becomes the continued fraction `[0; n, c1, c2, ..., 1]` — a real
number that *is* the trajectory. Start 27 (112 steps, peak 9232) is a
rational with a 278-digit denominator whose byte message surfaces
fragments like `y<.[-.E"i4.O'`; start 75 ends in silence (`"J.p..T..?......"`)
because the final crash through the powers of two is all unprintable
bytes. The punchline: the CF code spends 8–17 bits per step on a
trajectory whose true entropy is just `log2(n)` — the start number
determines everything. Collatz orbits are even worse for a CF coder than
e: Kolmogorov complexity versus what any digit-based coder can see.

`planck_demo.jl` closes the trilogy of number types. Since the 2019 SI
redefinition, Planck's constant is *exact by decree* — `h = 6.62607015e-34`
is rational, so its continued fraction terminates: the quantum of action
is a finite message of 15 quotients, 58 bits. The reduced constant
`hbar = h/2pi` inherits pi's irrationality and re-encodes at 3.43
bits/quotient, indistinguishable from Gauss–Kuzmin noise. By decree h is
a finite message; hbar talks forever and says nothing.

`constants_demo.jl` adds three more constants, three more morals.
Champernowne's `0.123456789101112...` is provably normal in base 10, yet
its CF hides a 166-digit quotient that overflows the decoder's `Int64` —
typical in one number system, pathological in another. Khinchin's theorem
(the geometric mean of almost every real's quotients converges to
K = 2.6854...) is checked live on π's quotients — while whether K obeys
its own theorem remains open. And the fine-structure constant
`1/alpha = 137.035999177(21)`: recomputing the CF at both ends of the
error bar shows only 6 quotients survive the experimental uncertainty —
h is exact by decree, π by mathematics, but α's continued fraction ends
at the precision of experiment.

A measured head-to-head against PadicCoding.jl — where each coder wins
and why — is in [COMPARISON.md](COMPARISON.md) (`compare.jl` reproduces
it).

Run tests and demos with

```
julia --project -e 'using Pkg; Pkg.test()'
julia --project demos/demo.jl
julia --project demos/pi_demo.jl
julia --project demos/e_demo.jl
julia --project demos/collatz_demo.jl
julia --project demos/planck_demo.jl
julia --project demos/constants_demo.jl
```

## Notes

- The whole message is coded as one continued fraction with exact `BigInt`
  convergents (`k_n` grows like `e^(n pi^2 / 12)`), so this is a batch
  coder, not a streaming one — no PR/AR-style renormalization. Fine for
  fun-sized messages; encoding 20,000 symbols takes well under a second.
- The emitted dyadic point is within one bit of the shortest one in the
  cylinder: the encoder takes the first grid level finer than the cylinder
  width and then strips trailing zero bits.
- No external dependencies (only `Random`, `Printf`, `Test` from the
  standard library).
