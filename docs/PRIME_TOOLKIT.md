# The prime toolkit: every tool in `numtools`, and the structures it found

**Walk the primes with pure structure.** Deterministic logic settles the bulk of every range, and
a real primality test is spent only at the edge where the free lanes run out. This is the index
of the kit; [RAILS_AND_LANES.md](RAILS_AND_LANES.md) is the thinking, [COOKBOOK.md](COOKBOOK.md)
the recipes with their output.

```python
from sitoolz import numtools as nt
```

## The one idea

Every prime past 3 is `6k ± 1`. A small wheel of prime lanes `{2, 3, 5, 7, 11, 13}` divides out
most numbers for free. Below the horizon, the square of the next prime past the wheel, `289`,
surviving the wheel **is** being prime. The only numbers that slip past are semiprimes with both
factors at least 17, and a stack of chain-walks factors those deterministically too.

Why factoring is hard, in one line: a wheel that held every prime as a lane would catch every
composite by its smallest factor. The difficulty lives in the gap the wheel cannot afford to
fill, and the fold-wheel reaches into that gap by manufacturing the missing lanes from structure,
`a² − b²` and a multiplier, instead of storing them.

Measured on every number to 5,000: **87.8% deterministic logic, 0% heuristic, 12.2% one test**,
and that 12.2% is only the primes past the horizon.

## The atoms: primality and factoring

```python
nt.is_prime(n)         # six-wheel primality
nt.next_prime(n)       # the smallest prime above n
nt.prev_prime(n)       # the largest prime below n
nt.prime_pi(n)         # the count of primes to n
nt.nth_prime(k)        # the k-th prime, nth_prime(1) == 2
nt.factor(n)           # 360 -> [(2, 3), (3, 2), (5, 1)]
```

## The verdict: free where it can, a test only where forced

```python
nt.verdict(x)          # (verdict, how, cost): free lanes -> fermat -> fold-wheel -> test
nt.free_verdict(x)     # the deterministic part only
nt.lane_verdict(x)     # the wheel and the horizon
nt.lane_read(x)        # every deterministic lane for x
```

## The chain-walks: factoring the slip band

```python
nt.is_square(x)        # a perfect square on sight: the gap-0 factors
nt.fermat_walk(x)      # the near-square walk: close factor pairs
nt.rho_walk(x)         # Pollard rho: far pairs, kept as the deep fallback
nt.fold_wheel(x)       # the multiplier fold: a far pair folded close, then Fermat -> (p, q, k)
nt.fold_sharp(x)       # fold_wheel sharpened by the n mod 4 fork, about 30% fewer checks
```

## The rails: walkable DataFrames

```python
nt.master_rail(lo, hi, restart=True)   # everything: rails, kind, factors, both square views, the walk
nt.walk_rail(lo, hi, lanes=...)        # the walk with the settling method per row
nt.range_rail(lo, hi, restart=True)    # the left rail restarted from 0 on a range
nt.verdict_frame(hi)                   # the filled verdict table from 0
nt.factor_frame(lo, hi)                # the square-factoring walk: n = a² − b², squares at b = 0
nt.full_rail(hi, base10=True)          # the three rails, twin and sexy, the base-10 double-up
nt.five_rail(hi)                       # the prime rail through the lens of 5
```

Restartable rails reset `count / counn / railn` to 0 at `lo`, so any window can be studied on
its own:

```python
nt.range_rail(500, 512, restart=True)    # count starts at 0 at 500
nt.range_rail(500, 512, restart=False)   # the global tally from 0
nt.master_rail(96, 132)                  # everything, restarted, across the 2^7 neighbourhood
```

## The structure frames

```python
nt.twin_frame(lo, hi)          # twin pairs (p, p+2): the shared six-cell, the even front
nt.sexy_frame(lo, hi)          # sexy pairs (p, p+6): the lane they share
nt.gap_frame(lo, hi)           # consecutive gaps with running records
nt.cunningham_frame(seed)      # a Cunningham chain, 2p+1, rung by rung
nt.octave_frame(kmax)          # 2^k−1 / 2^k / 2^k+1: the Mersenne and Fermat shoulders
nt.fold_frame(lo, hi)          # every slip semiprime: factors, ratio, fold k, a, b, mod 4
nt.two_square_frame(lo, hi)    # every prime as a² + b², which exists iff p = 2 or p ≡ 1 mod 4
nt.sos_walk(lo, hi)            # the sum-of-squares walk: every a² + b² and the count r2
nt.sum_of_squares(n)           # every n = a² + b² with 1 <= a <= b
nt.slips(hi)                   # the slip-through semiprimes past the horizon
```

## The analysis tools

```python
nt.describe(n)                 # a full dossier on one number
nt.dashboard(lo, hi)           # the logic-versus-test accounting for a range
nt.wheel_compare(wheels)       # horizon, first slip, free share and fold recovery per wheel
```

## The structures found

**The mod-4 fork forks both signs.** One residue governs factoring and primality together:

| | `n ≡ 1 (mod 4)` | `n ≡ 3 (mod 4)` |
|---|---|---|
| composite | `a² − b²`, a odd, b even: factors | `a² − b²`, a even, b odd: factors |
| prime | `a² + b²`, unique (Fermat's two squares) | neither: irreducible both ways |

Difference of squares peels a composite open; sum of squares is the fingerprint a `1 mod 4`
prime carries. The sign is the fork. Zero exceptions to 4,000.

**Sum of two squares, the full theorem.** Any `n` is `a² + b²` exactly when every prime
`≡ 3 (mod 4)` in its factorization appears to an even power. Primes `≡ 1 mod 4` get one
representation each and products multiply the count: `65 = 5·13 = 1² + 8² = 4² + 7²`,
`325 = 5²·13 = 1² + 18² = 6² + 17² = 10² + 15²`. `nt.sos_walk` walks the whole pattern. Zero
exceptions to 500.

**The fold multiplier tracks the factor ratio in octaves.** `k = 8 → ratio ≈ 3`, `16 → 4.8`,
`32 → 8.6`, `64 → 13.8`. Powers of two make the fold a bit shift, and the power-of-two ladder
brackets any ratio by its logarithm. Any `k ≡ 2 (mod 4)` makes `k·n ≡ 2 (mod 4)`, never a
difference of squares, so those multipliers are barred for free.

**The six-lane predicts the factor lanes.** `n ≡ 1 (mod 6)` is `(6k−1)(6k−1)` or
`(6k+1)(6k+1)`; `n ≡ 5 (mod 6)` is always `(6k−1)(6k+1)`. Exact.

**The powers of two.** `2^k` never flips the rail. The primes live on the odd shoulders:
Mersenne `2^k − 1` (31, 127, …) and Fermat `2^k + 1` (17, 257, …). The Mersenne side is
divisible by 3 exactly when `k` is even, the Fermat side exactly when `k` is odd.

## Recipes

```python
nt.master_rail(120, 135)                     # the everything-walk across a power-of-two shoulder
nt.range_rail(10_000, 10_060, restart=True)  # a far window with a fresh local rail
nt.describe(1717)                            # 17×101, fold k=16, 59²−42², verdict free
nt.dashboard(2, 20_000)                      # the accounting for a range
nt.fold_frame(289, 3000)[['n', 'p', 'q', 'ratio', 'k']]   # the multiplier climbing with the octave
nt.two_square_frame(2, 200)                  # primes wearing their two squares
nt.wheel_compare([(2,3,5,7,11,13), (2,3,5,7,13), (2,3,5,7,11,13,17,19)])   # a weird wheel's early hole
```

## Where the certainties end

The wheel, square, Fermat and fold lanes are mod-class certainties, free and never wrong. Rho is
the one heuristic, kept as a deep fallback that now almost never fires. A number that survives
every walk is prime, and confirming it is the one real cost left: a thin band of primes past the
horizon, which shrinks every time the wheel grows.

*Free first. Walk the slip band. Pay one test only where the walks run dry. The gap the wheel
cannot fill is the whole reason there is a game to play.*
