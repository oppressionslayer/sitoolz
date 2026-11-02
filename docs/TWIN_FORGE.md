# The twin forge: twin primes you can prove

**Find twin primes at any size, and build the big ones so both sides carry their own proof.**
Three ways up: sieve and read them off, search the `6k ± 1` lane, and forge `k·2ⁿ ± 1`, the
shape that makes both members provable because `k·2ⁿ` is factored by construction. The forge
runs in parallel, streams what it proves, and shipped an 82-digit pair proven on both sides,
a 623-digit and a 1,240-digit pair with the `+1` side proven, and 187 proven Proth primes along
the way.

## One: sieve and read them off

Twins are consecutive primes two apart. Build the primes, look for the 2s.

```python
from sitoolz import gap_forge
gap_forge.GapForge(1_000_000).twins       # lower members of every twin pair below a million
# [3, 5, 11, 17, 29, 41, 59, 71, ...]

from sitoolz import primebase as pb
pb.twins(0, 100)                          # [[3, 5], [5, 7], [11, 13], [17, 19], ...]
```

## Two: the 6k lane, any size in a second

Every twin pair above `(3, 5)` straddles a multiple of 6: it is `(6k − 1, 6k + 1)`. So test
random multiples of 6 and look both ways.

```python
import random
from sitoolz import numtools as nt

def find_twin(digits, rng=random):
    lo, hi = 10**(digits - 1) // 6, 10**digits // 6
    while True:
        m = 6 * rng.randrange(lo, hi)
        if nt.is_prime(m - 1) and nt.is_prime(m + 1):
            return m - 1, m + 1

find_twin(20)
# (70990165992303087377, 70990165992303087379)
```

Below 24 digits the test is a deterministic proof. Above that it is a strong probable-prime
screen, and a pair found this way is a candidate for the forge.

## Three: the forge, `k·2ⁿ ± 1`

Fix `n`, search a small `k` so that both `k·2ⁿ ± 1` are prime, and both members become
provable, because `k·2ⁿ` is fully factored the moment you write it down:

- `k·2ⁿ + 1` is a Proth number. `N − 1 = k·2ⁿ` is factored, so one Pocklington witness proves it.
- `k·2ⁿ − 1` has `N + 1 = k·2ⁿ` factored, so the Lucas / Brillhart–Lehmer–Selfridge test proves it.

Two engines meeting on one twin. This is the shape the largest known twins take; the
record pair `2996863034895·2¹²⁹⁰⁰⁰⁰ ± 1` has 388,342 digits.

`sitoolz.twin_forge` sieves `k` hard on both sides, screens the survivors with a probable-prime
test, and on a hit proves both sides outright, streaming each result as a packet:

```python
from sitoolz.twin_forge import estimate, hunt
estimate(256)                        # how big is this search at n = 256?
hunt(256, out="twins_256.jsonl")     # sieve, screen, prove, stream
```

### The 82-digit pair, proven both ways

`k = 19902 = 2·3·31·107`, `n = 256`. The two members:

```
k·2²⁵⁶ − 1 = 2304494160001066921319909743642906742095779234815578505513284836925487106094006271
k·2²⁵⁶ + 1 = 2304494160001066921319909743642906742095779234815578505513284836925487106094006273
```

The midpoint `k·2²⁵⁶` is divisible by 3, as every twin midpoint above `(3, 5)` must be. The
`+1` side is proven by Pocklington with the witnesses `2 → 5, 3 → 2, 31 → 2, 107 → 2`, one base
per prime factor of `N − 1 = 2²⁵⁷·3·31·107`, and the factored part `2²⁵⁶` already exceeds `√N`.
The `−1` side is proven by Lucas with a sequence `(P, Q)` whose discriminant has Jacobi symbol
`−1` modulo `N`: `U_{N+1} ≡ 0` and `gcd(U_{(N+1)/q}, N) = 1` for each prime `q` dividing `N + 1`.

Run the whole proof yourself in under a second. The Lucas prover self-tests against every
prime and composite below 4,000 before it is trusted:

```python
from sitoolz import twin_k2n256
twin_k2n256.verify()
# {'k': 19902, 'n': 256, 'digits': 82, 'difference': 2, 'midpoint_divisible_by_3': True,
#  'lucas_selftest': True, 'upper_proven': True, 'lower_proven': True,
#  'witnesses': {2: 5, 3: 2, 31: 2, 107: 2}, 'seconds': 0.08}
```

### The hunts at 2048, 4096 and 8192

Three runs of the forge, each proving every `+1` side it reached and screening the `−1` side
for a partner. The records ship with the package.

| n | digits | k tested up to | `+1` sides proven | pairs found |
|---:|---:|---:|---:|---:|
| 2,048 | 623 | 598,578 | 61 | 1, `k = 442533` |
| 4,096 | 1,240 | 2,390,088 | 121 | 1, `k = 2130693` |
| 8,192 | 2,472 | 435,690 | 5 | the hunt continues |

Both pairs have the `+1` member proven by Pocklington and the `−1` member a strong probable
prime. The Lucas proof for the `−1` side is the same call the 82-digit pair used, one call away
from a fully certified 1,240-digit twin.

```python
from sitoolz import showcase
showcase.twin_pairs(4096)
# [{'k': 2130693, 'n': 4096, 'digits': 1240, 'upper': ..., 'upper_proven_pocklington': True,
#   'lower': ..., 'lower_status': 'probable_prime'}]
len(showcase.twins(4096))            # 121 proven Proth primes of 1,238 digits
```

A fourth run aimed at `n = 332193`, a 100,001-digit twin, sieved 2,000,000 values of `k` down to
11,215 survivors and tested them at one every five minutes on four cores. That is the size
where the arithmetic wants a faster big-integer backend, and the forge is ready for one.

## Where each method belongs

| size | method | what you get |
|---|---|---|
| to a million | sieve | every pair, instantly |
| to 24 digits | the 6k lane | a pair with a deterministic verdict |
| any size | the 6k lane | a strong probable-prime pair |
| any size | the forge | a pair with proofs on both sides |

*Build the primes, watch for the 2s, and prove both sides.*
