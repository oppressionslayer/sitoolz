# The prime cookbook: worked recipes with their output

**Copy a block, run it, and see what it prints.** Every recipe is one of `numtools`' walks,
frames or verdicts; [PRIME_TOOLKIT.md](PRIME_TOOLKIT.md) is the full index and
[RAILS_AND_LANES.md](RAILS_AND_LANES.md) the ideas behind them.

```python
from sitoolz import numtools as nt
import pandas as pd
pd.set_option("display.width", 180)
```

## 1. Inspect one number

```python
nt.describe(1717)
```
```
{'n': 1717, 'prime': False, 'lane': '6k+1', 'mod4': 1, 'mod6': 1, 'mod8': 5,
 'kind': 'semiprime', 'factors': '17×101', 'diff_of_squares': '59^2-42^2',
 'fold_k': 16, 'verdict': ('composite', 'fold-wheel k=16: 17x101', 'free')}
```
`1717 = 17 × 101`, folded open at multiplier `k = 16` (a 4-bit shift), settled free with no test.

## 2. The mod-4 fork, right in the walk

```python
nt.master_rail(37, 50)[["x", "prime", "kind", "factors", "square_view"]]
```
```
 x  prime      kind factors square_view
37   True     prime      37       1²+6²      prime 1 mod 4:   a² + b²
38  False semiprime    2×19
39  False semiprime    3×13       8²−5²      composite:       a² − b²
41   True     prime      41       4²+5²
43   True     prime      43                  prime 3 mod 4:   neither
45  False composite   3^2×5       9²−6²
49  False    square     7^2       7²−0²      a square is just b = 0
```
The sign in `square_view` is `n mod 4`: composites subtract, `1 mod 4` primes add, and `3 mod 4`
primes have no split at all.

## 3. Restart the rail on any range

```python
nt.range_rail(10_000, 10_030, restart=True)    # count resets to 0 at 10000
nt.range_rail(10_000, 10_030, restart=False)   # count is the global tally from 0
```
`restart=True` gives a fresh local prime count for the window, so any stretch can be studied on
its own.

## 4. The fold multiplier grows with the factor ratio

```python
for n in (323, 1717, 2533):
    r = nt.fold_sharp(n)
    print(f"{n} = {r[0]}×{r[1]}   ratio {r[1]/r[0]:.1f}   k={r[2]}"
          + (f" (shift {int.bit_length(r[2])-1} bits)" if (r[2] > 1 and (r[2] & (r[2]-1)) == 0) else ""))
```
```
323 = 17×19     ratio 1.1   k=1                 close pair, no fold
1717 = 17×101   ratio 5.9   k=16 (shift 4 bits)
2533 = 17×149   ratio 8.8   k=32 (shift 5 bits)
```
The further apart the two factors, the bigger the multiplier that folds them close, and when `k`
is a power of two the fold is a bit shift. `fold_sharp` tries small `k` first and takes the
first hit, so the bands overlap as `n` grows; the trend is the structure.

## 5. The powers of two: Mersenne and Fermat shoulders

```python
nt.octave_frame(9)[["name", "n", "prime", "mod3"]]
```
```
 name   n  prime  mod3
2^k-1   7   True     1
  2^k   8  False     2
2^k+1   9  False     0
2^k-1  31   True     1     Mersenne prime
2^k+1  17   True     2     Fermat prime
2^k-1 127   True     1     Mersenne prime
2^k+1 257   True     2     Fermat prime
```
`2^k` never flips the rail. The Mersenne side is divisible by 3 exactly when `k` is even, the
Fermat side exactly when `k` is odd, read straight off the `mod3` column.

## 6. Twins, gaps, chains

```python
nt.twin_frame(0, 60)       # pairs (p, p+2) and the shared six-cell k = (p+1)/6
nt.gap_frame(2, 60)        # consecutive gaps with running records
nt.cunningham_frame(89)    # 89 -> 179 -> 359 -> 719 -> 1439 -> 2879, five doublings
```
```
# gap records to 60:  gap 1 at 2, gap 2 at 3, gap 4 at 7, gap 6 at 23
# cunningham(89): all six rungs 6k-1, all safe to double, chain breaks at 5759 = 13×443
```

## 7. The accounting for any range

```python
nt.dashboard(2, 10000)
```
```
               category  count  pct          kind
   wheel (small factor)   8076 80.8 deterministic
  below horizon (prime)     61  0.6 deterministic
square/fermat/fold walk    693  6.9 deterministic
        heuristic (rho)      0  0.0     heuristic
     brute test (prime)   1168 11.7    brute test
```
**88.3% deterministic logic, 0% heuristic, 11.7% one test**, and the test slice is only the
primes past the horizon.

## 8. Weird wheels

```python
nt.wheel_compare([(2,3,5,7,11,13), (2,3,5,7,13), (2,3,5,7,11,13,17,19)])
```
```
              wheel  horizon  first_slip  free_pct  slips fold_recovers
      2,3,5,7,11,13      289         289      81.7    373       373/373
         2,3,5,7,13      289         121      79.3    486       473/486    skip 11: a hole opens at 121
2,3,5,7,11,13,17,19      529         529      84.4    249       249/249
```
Skip a prime and the first hole jumps to that prime squared, long before the naive horizon,
and the fold lane still recovers nearly all of it. Grow the wheel and the free band climbs.
A wheel that held every prime cheaply would leave nothing hard, and that gap is exactly where
factoring lives.

## 9. Primes as a sum of two squares

```python
nt.two_square_frame(2, 40)
```
```
 p  m4  sum_sq   representable
 5   1  1^2+2^2      True
 7   3  —            False
13   1  2^2+3^2      True
17   1  1^2+4^2      True
```
`p = 2` or `p ≡ 1 (mod 4)` gives a unique `a² + b²`; `p ≡ 3 (mod 4)` gives none. Fermat's
two-square theorem, the prime side of the same fork.

## 10. Counting and navigation

```python
nt.prime_pi(1000)      # 168 primes to 1000
nt.nth_prime(100)      # 541
nt.next_prime(10**6)   # 1000003
nt.factor(360)         # [(2, 3), (3, 2), (5, 1)]
```

*Free first. Walk the slip band. Pay one test only where the walks run dry.*
