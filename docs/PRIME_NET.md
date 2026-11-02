# The prime net: prove, factor and generate to 72 quadrillion from 15 MB

**Hold every prime to 268,499,993 in a 15 MB net and answer about any number to
72,092,246,241,000,049 from it.** `sitoolz.primebase` proves a sixteen-digit number prime in
an instant, factors it, hands you every prime in a window past a trillion, counts primes
exactly, maps the factor anatomy around a number, and rolls primes up and down the
Pocklington ladder.

```python
from sitoolz import primebase as pb
```

Nothing loads at import. The net is read on first call, and small queries decode only the
prefix they need.

## What it answers

```python
pb.first(31)                       # the first 31 primes -> 2, 3, 5, ..., 127
pb.nth(1_000_000)                  # the millionth prime -> 15485863
pb.primes(10**12, 10**12 + 100)    # every prime just past a trillion
pb.twins(0, 100)                   # [[3, 5], [5, 7], [11, 13], [17, 19], ...]

pb.is_prime(9_999_999_999_999_917) # True, and proven: no factor in the net means prime
pb.factor(72_092_246_240_999_997)  # 3, the smallest factor of a composite
pb.factorize(1_009_732_533_765_288)# {2: 3, 3: 1, 61: 1, 197207: 1, 3497381: 1}

pb.pi(10**9)                       # 50847534, the exact count of primes to a billion
pb.count_range(10**12, 10**12 + 10**6)   # exact primes in a window
pb.octaves(20, 24)                 # primes per [2^k, 2^(k+1)), beside the Li estimate
pb.li(10**12)                      # the logarithmic-integral estimate

pb.neighborhood(N, width=5)        # factor anatomy around N, with a proof anchor
pb.descend(9_999_999_999_999_917)  # roll a prime down its Pocklington ladder to a seed
pb.roll_up(73)                     # build larger proven primes N = 2·R·73 + 1
```

## How a 15 MB file proves a 17-digit number

A composite `n` has a prime factor at most `√n`. The net holds every prime to `√reach`, so a
number in range with no factor in the net is prime, by elimination rather than by a test.
The net ships as `uint8` half-gaps between consecutive primes, eight times smaller than the
primes themselves, and decodes in about 150 ms. `primes(lo, hi)` runs a segmented sieve
seeded from the net, so a window far out is as quick as one near zero.

## Pick your reach

Reach grows as the square of the net's top prime, so storage grows only as the square root of
the reach. The shipped row covers most work; the rows below it are one build call away.

| net holds primes to | provable reach | primes | on disk | RAM to build |
|---|---|---:|---:|---:|
| **2.68 × 10⁸** | **7.2 × 10¹⁶**, 72 quadrillion | 14.6 M | **15 MB** | 0.3 GB, ships |
| 1.8 × 10⁹ | 3.3 × 10¹⁸ | 85 M | 85 MB | 2 GB |
| 1.0 × 10¹⁰ | 1.0 × 10²⁰ | 434 M | 434 MB | 10 GB |
| 3.2 × 10¹⁰ | 1.0 × 10²¹ | 1.35 B | 1.35 GB | 32 GB |

`is_prime` and `factor` answer within the reach and say so past it. Counting goes further:
`pi(x)` is exact through the Lucy_Hedgehog method, `O(x^0.75)` with no prime list, into the
10¹² range, and `li(x)` estimates any `x` to well under 0.01%.

## Build a bigger net

```python
from sitoolz.primebase.build import build_base, pack
p = build_base(1_800_000_000)      # sieve to 1.8 × 10⁹: about 2 GB of RAM and 30 s
pack(p, "big_base")                # big_base.npy and big_base.meta.json
```

Point the engine at it before importing:

```bash
export PRIMEBASE_GAPS=/path/big_base.npy
export PRIMEBASE_META=/path/big_base.meta.json
```

For memory-mapped access to a raw array instead of the packed form:

```python
from sitoolz.primebase.build import export_raw
export_raw(pb.base(), "raw_base")                  # raw_base.npy, int64
import numpy as np
primes = np.load("raw_base.npy", mmap_mode="r")    # only the pages you touch hit RAM
```

## The ladder, both directions

`descend` factors `N − 1` to recover the prime a rollball prime was built from, and keeps
going down to a seed. `roll_up` builds the next proven primes above `q` as `N = 2·R·q + 1`.
`sitoolz.primebase.prove` carries the full construction: `build_up(seed, bits)`,
`certificate(N)`, `verify_certificate(steps)`. The 19,538-digit prime in
[BUILT_NOT_FOUND.md](BUILT_NOT_FOUND.md) is the ladder at full stretch.

## Tests

The engine ships its own 67 tests under `tests/primebase/`: the net decodes to the sieve,
counting agrees three ways, every proof and factorization replays, and streaming builds
match the packed form.
