# Rails, lanes, gates and folds: the primes as structure you can walk

**Read the primes off a grid of free logic, settle most numbers without a single primality
test, and pay only at the edge.** The rails hang every prime on a power of two. The lanes are
residue classes that cost a couple of XORs and never lie. The gates forbid composites by
construction. The folds factor what slips through. And the base-prime bootstrap turns a small
set of primes into every prime below its square.

```python
from sitoolz import numtools as nt
```

## Every-time lanes and some-of-the-time lanes

Two kinds of prime logic, and telling them apart is the whole game:

- **Every-time lanes** are pure functions of the low bits of `x`: `x mod 8`, `(x+3)^5 < x`,
  which side of a twin `(x+1)^5` lands on. They cost a couple of XORs, they are never wrong,
  and they answer every single time.
- **Some-of-the-time lanes** are the tests: `is_prime(x)`, `is_prime(2x+1)`. They pay off
  sometimes, and finding out costs work.

Grade each on its own exam. "How often is `(x+3)^5 == x` prime?" gives 29%, which looks like
noise, and it is the wrong question. `(x+3)^5 == x` fires every time `x ≡ 1, 3 (mod 8)`, for
nothing. That is free, exact structure.

### The horizon, where free runs out

Take a wheel of lanes `[2, 3, 5, 7, 11, 13]`. A number no lane divides is a survivor. Below the
square of the next prime past the wheel, here `17² = 289`, survival **is** primality. In that
band the cheap mod-logic reproduces the prime walk with zero tests. The first number to slip
past is `289 = 17 × 17`, and every slip past the horizon is a semiprime `p·q` with both factors
at least 17. Add a lane and the horizon jumps to the next square.

| the free lanes say | verdict | cost |
|---|---|---|
| a lane divides `x` | composite | free, certain |
| survives all lanes, `x < horizon` | prime | free, certain |
| survives, `x ≥ horizon` | a slip-through semiprime could reach here | one walk, or one test |

With six lanes, 82% of all numbers to 5,000 get a free, certain verdict.

```python
nt.lane_read(x)        # every deterministic lane for x: mod classes and fold directions
nt.lane_verdict(x)     # (verdict, reason, cost): 'free' wherever it can, 'test' only where forced
nt.logic_walk(hi)      # the prime walk from free lanes alone, exact below the horizon
nt.slips(hi)           # the slip-through semiprimes, exactly where more work is needed
```

## The three chain-walks past the horizon

Past the horizon the only survivors are semiprimes, and each walk reaches a different gap
between the two factors. Sort the slips by factor gap and the catcher marches straight across:

```python
nt.is_square(x)        # gap 0: equal factors p², one isqrt
nt.fermat_walk(323)    # (17, 19)   close pairs: walk a down from √n, a step or two
nt.rho_walk(1717)      # (17, 101)  far pairs: a pseudo-random loop and a gcd
```

The tiling from 289 to 5,000, 295 slip-throughs, every one factored, no test:

```
gap 0        13   all square
gap 2–8      21   all fermat
gap 9–20     35   all fermat
gap 21–50    64   fermat 63 · rho 1
gap 51–120   92   fermat 17 · rho 75      the handoff
gap 121+     70   all rho
```

### The fold-wheel retires the heuristic

Rho is a walk that almost always wins fast. The fold-wheel replaces it with certain logic:
fold a far pair close by multiplying `n` by `k` so that `k·p ≈ q`, then Fermat lands it. The
multipliers are the wheel, and powers of two make the fold a bit shift.

```python
nt.fold_wheel(n)       # deterministic multiplier fold -> (p, q, k)
nt.fold_sharp(1717)    # (17, 101, 16)   k = 16 is a shift left by 4
nt.fold_sharp(2533)    # (17, 149, 32)   k = 32 is a shift left by 5
```

Two exact `mod 4` facts drive `fold_sharp`. A number `≡ 2 (mod 4)` is never a difference of two
squares, so any multiplier `k ≡ 2 (mod 4)` can never land and is skipped. For odd `k·n` only one
parity of `a` gives an integer `b`: `a` odd when `k·n ≡ 1`, `a` even when `k·n ≡ 3`, which halves
the checks. The multiplier tracks the factor ratio `q/p` in octaves: `k = 8 → ratio ≈ 3`,
`16 → 4.8`, `32 → 8.6`, `64 → 13.8`.

With the fold-wheel in the stack, every number to 5,000 settles as **87.8% deterministic logic,
0% heuristic, 12.2% one test**, and that last band is exactly the primes past the horizon.
Grow the wheel to shrink it.

```python
nt.verdict(91)         # ('composite', '7 divides it', 'free')
nt.verdict(323)        # ('composite', 'near-square walk: 17×19', 'free')
nt.verdict(1717)       # ('composite', 'fold-wheel k=16: 17×101', 'free')
nt.verdict(4999)       # ('prime', 'confirmed by test', 'test')
nt.free_verdict(x)     # the same, stopping before any walk: 'free' or 'test' only
nt.verdict_frame(hi)   # the filled table from 0: rails, wheels, structure, the settling method
```

## The walkable rail

```python
nt.walk_rail(lo=0, hi=120, lanes=(2, 3, 5, 7, 11, 13))   # a DataFrame, one row per x
```

Three rails roll down the left, `count / counn / railn`, lagged 0, 1 and 2, each flipping only
when a prime lands. `count` is `π(<x) − 1`. Beside them: `x`, `m5`, `lane` (the six-wheel lane,
`6k−1`, `6k+1` or `·`), `prime`, `twin`, `sexy`, `walk` (which walk settled it: `wheel`, `below`,
`square`, `fermat`, `rho`, `prime?`), `fac`, and `cost` (`free`, `cheap`, `test`).
`df.attrs['horizon']` carries the horizon for the wheel you chose.

```python
nt.walk_rail(284, 332)                                  # walk straight across the horizon 289
nt.walk_rail(122, 135)                                  # the 2^7 neighbourhood: 127 lands 'below'
nt.walk_rail(286, 296, lanes=(2, 3, 5, 7, 11, 13, 17, 19))   # horizon 529; 289 is now plain 'wheel'
```

| wheel | horizon | free verdicts to 3,000 |
|---|---:|---:|
| `2, 3, 5, 7` | 121 | 82.0% |
| `2, 3, 5, 7, 11, 13` | 289 | 85.9% |
| `2, 3, 5, 7, 11, 13, 17, 19` | 529 | 88.5% |

```python
# the powers of two: Mersenne and Fermat primes on the odd shoulders
for k in (5, 6, 7, 8):
    P = 1 << k
    print(nt.walk_rail(P - 6, P + 7)[['x', 'lane', 'prime', 'walk', 'fac', 'cost']])

nt.factor_frame(lo, hi)   # every number as a² − b²: a square is b = 0, a prime has no split
```

## The two rails: mod 4

Every odd prime is `1` or `3 (mod 4)`. `rail_climb_tree` hangs each on a side of a power of two:
`2^k − 1` catches the `3 mod 4` lane (**L**), `2^k + 1` the `1 mod 4` lane (**R**), and the
lone even prime rides **e**.

```python
df, info = nt.rail_climb_tree(0, 512, split=True)
info['rails']        # {'L': 51, 'R': 45, 'e': 1}
info['coverage']     # 1.0: every prime captured
```

Completeness is guaranteed: below 512 there are 96 odd primes, 51 on L, 45 on R, plus the 2.
L leads R and stays ahead, which is Chebyshev's bias in the prime race.

Each row carries the prime's structural coordinates: its `parent` anchor `2^k ∓ 1` and `wrap`
(parent XOR child), its complement `other = 2^k − child` and whether that is prime (a Goldbach
split of a power of two: `31 + 97 = 128`, `19 + 109 = 128`, `61 + 67 = 128`), `up = 2^(k+1) −
child`, the invariant `sums_p2` (always true), the descent seed `under = (child−1)/2`, the
ascent `over = 2·child + 1`, the Cunningham chain lengths `chain_up` and `chain_down`, and the
`gap` to the next prime.

### Chains: build by doubling

`over = 2·child + 1`, chained while prime, is a Cunningham chain, and exactly the Pocklington
rollball with `R = 1`: each rung is proven off the last with witness `a = 2`.

```
2 → 5 → 11 → 23 → 47                     length 5, breaks at 95 = 5·19
89 → 179 → 359 → 719 → 1439 → 2879       length 6, breaks at 5759 = 13·443
```

```python
df[(df.chain_down == 0) & (df.chain_up > 0)].sort_values('chain_up', ascending=False)   # chain bottoms
```

The general rollball frees `R`, so the ladder never stalls; see
[BUILT_NOT_FOUND.md](BUILT_NOT_FOUND.md).

### Gates: composites forbidden by residue

Each construction has a `mod p` gate. Mod 3 is the cleanest:

```
c ≡ 1 (mod 3)  ->  2c+1 and (c−1)/2 are always composite,  'up' can be prime
c ≡ 2 (mod 3)  ->  2c+1 and (c−1)/2 can be prime,          'up' is always composite
```

It generalises: `2c + 1 ≡ 0 (mod p)` exactly when `c ≡ (p−1)/2 (mod p)`, so every small prime
forbids one residue of `c`. Stack the gates and only the candidates that can be prime remain,
before any expensive step. That is the speed lever under the rollball.

### Folds: the gate-survivors climb and fold

`xor_fold_climb(k)` takes the octave `[0, 2^k)`, keeps the gate-survivors (coprime to the small
primes), climbs them by power-of-two XOR, and folds each node onto its `2^k − child` partner.

```python
df = nt.xor_fold_climb(8)             # octave 256, gates {3, 5, 7, 11, 13, 17}
print(nt.climb_tree_ascii(df))        # the fold tree
df[df.fold_survives]                  # 23 + 233 = 256, 83 + 173 = 256, 29 + 227 = 256, ...
df.sums_p2.all()                      # True: child + fold = 2^k on every row
nt.xor_fold_climb(8, fold_only=True)  # only the mutual pairs
```

Below the squaring bound the gate-survivors are exactly the primes. Gates to 17 remove
everything whose smallest factor is under 19, and the smallest composite left is `19² = 361`,
so every survivor under 256 is prime with no test at all: a tree of primes created by
elimination. Pick `k`, set the gates to about `2^(k/2)`, and the whole octave sits under the
bound, nesting octave inside octave: 256 ⊂ 65536 ⊂ 2³².

## The base-prime bootstrap

A base of primes to `P` finds every prime to `P²` with no new test: each base prime strikes
its multiples, any composite to `P²` has a factor at most `P`, and the survivors are prime.

```
256 primes (to 1619)        reach 1619² = 2,621,161
primes to 2.6 million       reach about 6.8 trillion
each round squares the reach
```

```python
base  = nt.ruler_sieve(2000)[:256].tolist()      # the base, largest 1619
reach = base[-1] ** 2                            # 2,621,161
known = 1_000_003
wraps = [2, 4, 6, 8, 10, 12, 16, 64, 128, 20, 24, 36, 60, 124]
cands = sorted({known ^ w for w in wraps if 2 < (known ^ w) < reach})
[c for c in cands if nt.is_prime_by_index(c, base)]   # [1000039], confirmed by the base alone
```

That is the division of labour: **XOR organises, addition discovers, and they compose.** A
wrap `w = p ^ q` connects primes into a rail tree, and the lean wrap set `roll8`, the multiples
of 4 to `2·hi` plus `{1, 2, 6}`, reaches every prime in every range tested through 65,536:
`{2, 6}` stitch the twins and cousins (`5^7 = 2`, `11^13 = 6`) and `1` catches the 2. A wrap
picks candidates; a small even wrap lands on a prime 8 to 18% of the time near a million, which
is twin and cousin density. The base set decides. `sitoolz.primebase` is this bootstrap carried
to 268 million base primes and a reach of 72 quadrillion; see [PRIME_NET.md](PRIME_NET.md).

```python
nt.rail_climb_tree(lo, hi, rail="roll8")   # the wrap tree
nt.roll_ball(lo, hi, base="double")        # wraps against a base
nt.ruler_sieve(n); nt.ruler_sieve_128(n)   # the segmented sieve
nt.prime_stream(limit); nt.prime_set(limit)
nt.is_prime_by_index(c, base); nt.primes_by_index(N)
nt.xor_grid(a, b); nt.wide()
```

## The thread in one line

**Rails** (mod 4) → **complements** (sum to `2^k`) → **chains** (double and prove) →
**gates** (mod p forbids) → **folds** (survivors climb and fold) → **the bootstrap**
(a base to `P` finds everything to `P²`). Structure plus construction gives primes created by
logic: sieve by elimination, build by Pocklington, and pay for a test only on the thin band of
primes past the horizon.

*Free first. Walk the slip band. Pay one test only where the walks run dry.*
