# Built, not found: a 19,538-digit prime in nine rolls

**Build a prime of any size from a seed you already trust, and carry its proof up with it.**
Lars Rocha's rollball starts at M127, rolls nine times, and lands on a 64,903-bit prime whose
certificate is the way it came up. No probable-prime test anywhere in the build. Pure Python.

```python
from sitoolz import showcase
p = showcase.big_prime()          # the number itself, 19,538 digits
str(p)[:30], str(p)[-30:]         # 488722292701599022162589797030 … 178852201085316611321153499213
```

## The move

Start from a prime `q` you can prove. Pick a multiplier `R` and form

```
N = 2·R·q + 1
```

With `R < q/2` the known factor `q` exceeds `√N`, so one witness `a` settles `N` outright by
Pocklington's theorem (1914):

```
a^(N−1) ≡ 1 (mod N)   and   gcd(a^((N−1)/q) − 1, N) = 1   ⟹   N is prime.
```

One modular exponentiation, one gcd, and `N` is proven. Set `q = N` and roll again. Each roll
about doubles the bit length, and `N − 1 = 2·R·q` stays fully factored at every rung, because
you chose every factor on the way up. Maurer (1995) and Shawe-Taylor formalised this as
recursive provable-prime construction; the rollball reaches it by feel and keeps it simple.

## The nine rolls from M127

M127 = 2¹²⁷ − 1 is proven by Lucas–Lehmer. Witness `a = 3` closed every rung.

| rung | bits | digits |
|---|---:|---:|
| seed M127 | 127 | 39 |
| roll 1 | 254 | 77 |
| roll 2 | 508 | 153 |
| roll 3 | 1,015 | 306 |
| roll 4 | 2,029 | 611 |
| roll 5 | 4,058 | 1,222 |
| roll 6 | 8,114 | 2,443 |
| roll 7 | 16,226 | 4,885 |
| roll 8 | 32,452 | 9,769 |
| roll 9 | **64,903** | **19,538** |

Build time: 14.8 hours of pure-Python big-integer arithmetic on one machine. Afterwards two
checks that share no machinery with the construction both say prime: Fermat base 2 and a
strong Baillie–PSW test.

```python
import gmpy2
N = showcase.big_prime()
pow(2, N - 1, N) == 1                    # True
bool(gmpy2.is_prime(gmpy2.mpz(N)))       # True
```

## The certificate is the distance you walked

For a prime you find by searching, factoring `N − 1` is the hard problem. For a rollball prime
it is free: you paced the distance out yourself.

```
distance from seed to top = 2·R₁ · 2·R₂ · … · 2·R₉       every Rᵢ you chose
```

So `N − 1` arrives fully factored. That one fact proves the prime, lets you walk it back down
to its seed, and makes it coprime to everything it was built from.

## Roll your own

`sitoolz.primebase.prove` rolls and verifies. From M31 to 1,024 bits is seven rungs in a fifth of
a second, each rung about doubling the bits with `R` just under `q/2`; every rung is a dict with
`N`, the previous prime `q`, and the witness `a`. The same ladder in Rust, `rusttoolz`, rolls
M127 to 4,096 bits in four seconds and to 8,192 in seven rungs.

```python
from sitoolz.primebase import prove
steps = prove.build_up(2**31 - 1, 1024)   # roll M31 up to 1024 bits: 7 rungs, 0.2 s
prove.verify_certificate(steps)           # True: every rung re-proved from the one below
steps[-1]["N"].bit_length()               # 1024
```

`numtools` carries the session's own rollers, including the pure-Python one that built the
big prime's method at smaller sizes:

```python
from sitoolz import numtools as nt
prime, cert = nt.roll_prime(700)          # a 211-digit proven prime with its certificate
nt.verify_certificate(cert)               # replay the chain from the seed
nt.pocklington_walk(seed=2**31 - 1, steps=10, rbits=12)   # watch the rungs land
```

A shipped 4,096-bit prime shows the chain shape at a glance: nine rungs from the seed 112543,
17 → 33 → 65 → 129 → 257 → 513 → 1,025 → 2,049 → 4,096 bits.

```python
c = showcase.chain_4096()
[rung[0].bit_length() for rung in c["chain"]]
```

## Split down, build up

Split any number by `(n − 1)/2`. On primes that halves the gaps and turns twin primes into
adjacent integers. It is the descent operator, and `m → 2m + 1` is the ascent:

```
build UP   (m → 2m+1):   2 → 5 → 11 → 23 → 47       a Cunningham chain
split DOWN ((n−1)/2):    47 → 23 → 11 → 5 → 2       the same chain, reversed
```

| direction | name | move |
|---|---|---|
| up: `q` with `2q + 1` prime | Sophie Germain | `×2 + 1`, the `R = 1` rollball |
| down: `p` with `(p − 1)/2` prime | safe prime | `(n − 1)/2`, the descent |

The general rollball is `R > 1`: when `2q + 1` is composite, pick another `R` and keep
climbing, which is why the ladder never stalls. And for any rollball prime, `(N − 1)/2 = R·q`
drops you straight onto its known factorization.

```python
nt.cunningham(89)                          # 89 → 179 → 359 → 719 → 1439 → 2879
from sitoolz import primebase as pb
pb.descend(9_999_999_999_999_917)          # roll a prime down its ladder to a seed
pb.roll_up(73)                             # the next proven primes N = 2·R·73 + 1
```

## Lineage

Euclid's `+1` → Pocklington's witness (1914) → Maurer and Shawe-Taylor's recursive
construction (1990s) → the rollball in sitoolz. Conceived and driven by Lars Rocha.
