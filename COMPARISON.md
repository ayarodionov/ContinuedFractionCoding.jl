# CF coding vs p-adic arithmetic coding

[![CI](https://github.com/ayarodionov/ContinuedFractionCoding.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ayarodionov/ContinuedFractionCoding.jl/actions/workflows/CI.yml)

A head-to-head of [ContinuedFractionCoding.jl](README.md) against
[PadicCoding.jl](https://github.com/ayarodionov/PadicCoding.jl), run by
`compare.jl` on five sources of 20,000 symbols each (5,000 for the π
quotients). The p-adic coder uses a static model built from the message's
own symbol counts, with unseen symbols in `1..max` padded to count 1;
rates are in bits per symbol.

| Source | Entropy | CF rate | p-adic rate |
|---|---|---|---|
| Gauss–Kuzmin (unbounded alphabet) | 3.39 | **3.44** | 10.79 |
| Skewed 16-symbol | 2.91 | 3.15 | **2.91** |
| Uniform 1..8 | 3.00 | 4.09 | **3.00** |
| Text bytes (PadicCoding README) | 5.06 | 12.58 | **5.10** |
| Partial quotients of π | 3.36 | **3.41** | 5.72 |

## Model-based vs fixed

The fundamental difference is that the p-adic coder separates *model* from
*coder*: give it any static frequency table and it codes at that model's
entropy, which is why it sits on the entropy line for every finite-alphabet
source above. CF coding has no model parameter at all — the code is frozen
into the number system, optimal for exactly one distribution, the
Gauss–Kuzmin law `P(k) = log2(1 + 1/(k(k+2)))`. On its one native source it
achieves 0.35% overhead with *zero* model — nothing to estimate, transmit,
or store — but on anything else it pays the full KL-divergence penalty: 8%
on the skewed source, 36% on uniform, 150% on text.

## Where CF genuinely wins: unbounded alphabets

The Gauss–Kuzmin and π rows are the striking ones. Those sources emit
integers with a `1/k^2` tail — the 20,000 samples contained values in the
tens of thousands. The p-adic coder (like any model-based arithmetic coder)
must enumerate its alphabet in advance, so the padded model hemorrhages
probability mass to thousands of never-used symbols and codes at 3x
entropy. CF coding handles any positive integer natively: a quotient of
size `k` intrinsically costs about `2 log2(k)` bits — a built-in universal
code for the integers, much like a Golomb code but tuned to the `1/k^2`
law. (A fair fight would give the p-adic coder an escape mechanism or
Golomb-style binning — a modeling fix, but engineering that CF gets for
free.)

## Everything else favors p-adic

- **Streaming.** The PR/AR rescaling machinery keeps the p-adic coder in
  fixed-width integers (`Int128`, or `UInt64` in `PadicCoding2`) with
  output emitted incrementally and interval width provably bounded below —
  it runs forever on constant memory. The CF coder is a batch coder: one
  giant continued fraction per message, `BigInt` convergents growing like
  `e^(n pi^2 / 12)`, the whole message in memory, output only at the end.
- **No CF renormalization exists.** There is no known CF analogue of PR/AR
  rescaling with a fixed word size — the Gauss-map cylinders have
  irrational-ratio widths and no `G(P^N)` grid to snap to. The closest
  thing (Stern–Brocot/mediant splitting) fails because the Gauss–Kuzmin
  distribution has infinite mean.
- **Speed.** `PadicCoding2` does ~23 MB/s on word operations; the CF
  coder's `O(n^2)`-bit arithmetic is fine for 20,000 symbols and hopeless
  beyond.
- **Flexibility.** The p-adic construction parameterizes the base (any
  prime P) and reproduces Huffman and Golomb–Rice codes as special cases;
  CF is one fixed system.

## The unifying view

Both are the same construction instantiated over different fibred number
systems, in Rényi's f-expansion sense. Both map a message to a nested
interval and pick a shortest dyadic point in it (`selectPoint` in the
p-adic paper and `dyadic_point` here are the same algorithm), and each is
the zero-overhead code for the source matching its number system's
invariant digit measure: dyadic sources for base 2, Gauss–Kuzmin/Zipf for
continued fractions. P-adic coding then adds the two things that make a
coder practical — a free model parameter and fixed-precision
renormalization — which is precisely what the CF system, for deep
number-theoretic reasons, cannot offer.

## Reproducing

`compare.jl` runs in the `compare/` environment, which declares
PadicCoding alongside this package. PadicCoding is not a dependency of
ContinuedFractionCoding itself, so it lives here rather than in the
top-level `Project.toml`.

```
julia --project=compare -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=compare compare.jl
```

The `develop` step points the environment at this checkout; it becomes
unnecessary once ContinuedFractionCoding is available from the General
registry, leaving just `Pkg.instantiate()`.
