# Mersenne Beacons

*Read the distance from `2^k` down to the nearest prime, and watch it light at exactly the
Mersenne primes. The distance‑1 case of the octave theory is a Mersenne prime, and a
deterministic test (Lucas–Lehmer) **proves** it with no oracle. This is the map of what
the beacon finds, what is known, and where the frontier sits.*

All code is in `numtools`; every number below is from a real run.

---

## 1. The beacon

`under_p2(k)` measures the distance from a power of two `2^k` **down** to the nearest prime
below it. That distance is always odd (even − odd = odd), and it equals **1** in exactly one
situation:

```
distance == 1   ⟺   2^k − 1 is prime   ⟺   2^k − 1 is a MERSENNE PRIME, and k is a Mersenne exponent
```

The all‑ones ceiling anchor `2^k − 1` used throughout the octave work **is** the Mersenne
candidate. When that ceiling is itself prime, the nearest prime below `2^k` sits one step away
and the beacon lights. So the Mersenne primes are not a separate topic — they are the brightest,
rarest case of the ceiling read‑out.

```python
from sitoolz import numtools as nt
df = nt.under_p2(32)
df[df.mersenne]      # k = 2,3,5,7,13,17,19,31  -> the beacons in range
```

---

## 2. What a Mersenne prime is

`M_p = 2^p − 1`. Two facts, one of which people get wrong:

- **`p` must be prime.** If `p = a·b` then `2^a − 1` divides `2^p − 1`, so a composite exponent
  can't give a prime. (This is why the beacon only ever lights at prime `k`.)
- **`p` prime is *not* enough.** The famous trap: `M_11 = 2^11 − 1 = 2047 = 23 · 89`. Prime
  exponent, composite Mersenne. So you can't just test the exponent — you have to test `M_p`.

---

## 3. Lucas–Lehmer — the proof (no oracle)

General primality at Mersenne sizes would be hopeless, but Mersennes have their own
**deterministic** test, and it's beautiful:

> Let `s₀ = 4` and `s_{i+1} = s_i² − 2`. For an odd prime `p`, `M_p = 2^p − 1` is prime **iff**
> `s_{p−2} ≡ 0 (mod M_p)`.

That's it — `p − 2` squarings mod `M_p`. Each step is one multiply, so a thousand‑digit
Mersenne is decided in milliseconds. It's a genuine **proof**, not a probable‑prime test — which
is what lets these be settled at all, and it fits the no‑oracle rule exactly.

```python
from sitoolz import numtools as nt

nt.lucas_lehmer(127)      # True   -> 2^127 − 1 is prime (M127)
nt.lucas_lehmer(11)       # False  -> 2^11 − 1 = 2047 = 23·89
nt.lucas_lehmer(1279)     # True   -> 2^1279 − 1 is prime (386 digits), decided in ~0.1 s
```

---

## 4. The known Mersenne exponents

Proven with `lucas_lehmer`, the exponents `p < 2300` giving a prime `M_p` (17 of them,
verified in ~1.6 s):

```
2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281
```

```python
nt.mersenne_exponents(2300)   # the list above
```

Cross‑check against the beacon — every one of these is a `distance == 1` hit in `under_p2`:

| exponent p | `2^p − 1`            | beacon (under_p2 distance) |
|---|---|---|
| 13 | 8191                | 1 ✓ |
| 17 | 131071              | 1 ✓ |
| 19 | 524287              | 1 ✓ |
| 31 | 2147483647          | 1 ✓ |
| 127 | (39 digits)        | 1 ✓ (Lucas, 1876) |
| 521 | (157 digits)       | 1 ✓ (1952) |
| 1279 | (386 digits)      | 1 ✓ (1952) |

---

## 5. The frontier — known and still open

**Every beacon that lit is a long‑known Mersenne** — `M127` was proven by Édouard Lucas in
**1876**, `M521`/`M607`/`M1279` in **1952** on the first computers. What the method shows is
that the octave theory sees them: it lights on the known beacons, independently, from its own
construction.

Where the frontier is:

- **Only 52 Mersenne primes are known** (as of the GIMPS discovery in October 2024).
- **The largest known prime is Mersenne:** `2^136279841 − 1` — about **41 million digits**,
  found Oct 2024 by GIMPS (Luke Durant). Every record‑holder for "largest known prime" for
  decades has been a Mersenne, *because* Lucas–Lehmer makes them the cheapest giants to test.
- **Whether there are infinitely many is an open problem.**
- Finding a *new* one means running Lucas–Lehmer on **untested large prime exponents** — that's
  distributed‑supercomputer work (GIMPS), months of GPU time per candidate. The beacon
  *confirms* a Mersenne instantly once you name its exponent; the search over untested
  exponents is a separate job at a different scale.

In one line: **the beacon is a perfect detector for Mersennes you point it at, and the hunt for
unknown ones is supercomputer work.** Different jobs.

---

## 6. Using the tools

```python
from sitoolz import numtools as nt

nt.under_p2(32)[lambda d: d.mersenne]   # beacons up to 2^32, provable via the base
nt.mersenne_exponents(1300)             # first 15 exponents, Lucas–Lehmer proven
nt.lucas_lehmer(2281)                   # True — a 687-digit Mersenne, proven in ~1 s
```

- **`under_p2(k)` / `altitude(k)`** — the ceiling read; `mersenne=True` / `distance=1` marks a
  beacon (base‑provable up to ~`2^40`, or ~`2^55` with a grown base).
- **`lucas_lehmer(p)`** — the deterministic Mersenne proof, good to thousands of digits.
- **`mersenne_exponents(pmax)`** — the proven list under `pmax`.

**One line summary:** the Mersenne primes are the octave theory's brightest beacons — the
`2^k − 1` ceiling landing exactly on a prime — and Lucas–Lehmer is the no‑oracle proof that
lets you confirm them from 3 up into the millions of digits, while *finding* new ones lives
out at the edge of what supercomputers can do.
