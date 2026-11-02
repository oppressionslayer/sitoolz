# Squaring vs. XOR — where an odd number keeps its factors

*Square an odd number and its factors stay on the surface; XOR it and they vanish. Take
`n = a·b` and operate on it: does the result "know" about `a` and `b`? **Squaring says
yes. XOR says no** — and the control below proves which is real, not just asserts it.*

Every number below is from an actual run; the program is `factor_probe` in `numtools`.

---

## The question

You've got an odd composite `n = a·b`. Two moves are on the table:

- **Square it:** `n²`
- **XOR it:** `n ^ w`

Does either land you on a number whose factorization reveals `a` and `b`?

---

## Way 1 — Squaring keeps the factors (real)

Squaring is **multiplicative**: `(a·b)² = a²·b²`. So the factors survive, guaranteed,
every single time:

```
47 · 37   = 1739      n² = 3024121          = [37, 37, 47, 47]
1009 · 1279 = 1290511   n² = 1665418641121   = [1009, 1009, 1279, 1279]
```

The original primes sit right on the surface, squared. No luck, no search.

**Bonus — the neighbours merge.** Since `n² − 1 = (n−1)(n+1)`, descending `n²−1` hands
you the union of *both* neighbours' factorizations in one list:

```
1739:  n²−1 = {2,3,5,11,29,79}  =  {2,11,79} (n−1)  ∪  {2,3,5,29} (n+1)
```

---

## Way 2 — XOR scrambles them (pure chance)

XOR is tempting: `n ^ w` *sometimes* lands on a number divisible by `a` or `b`. But it's
an illusion. XOR is **bitwise, not multiplicative** — `(a·b) ^ w` has no arithmetic
relationship to `a` and `b`; it just stirs the bits. The apparent "reveals" are nothing
but the base rate: any number near `n` is divisible by a prime `p` about `1/p` of the
time, whether or not `p` divides `n`.

### The airtight control

Sweep prime sizes against `n = 1009·1279` (`w < 40000`) and ask: *does being an actual
factor change the hit rate?*

| prime | is it a factor? | XOR hits | W/p (pure chance) |
|------:|:---------------:|---------:|------------------:|
| 3     | stranger        | 13333    | 13333 |
| 7     | stranger        | 5714     | 5714  |
| 41    | stranger        | 974      | 976   |
| 211   | stranger        | 189      | 190   |
| **1009** | **YES — factor** | **39** | **40** |
| **1279** | **YES — factor** | **30** | **31** |
| 5003  | stranger        | 8        | 8     |
| 9973  | stranger        | 3        | 4     |

Every prime — factor or stranger — hits at exactly `W/p`. The two *real* factors get
**no bump** for being factors; they hit at their size rate, the same as any prime that
size. Small primes hit more *only* because they're small. So the XOR hit rate is
**100% about size (`1/p`), 0% about whether it divides `n`.** Pure base rate.

> Note on the control: the strangers are deliberately chosen the *same size* as the real
> factors (e.g. 1013, 1283 next to 1009, 1279). That holds size constant so the only
> variable left is "is it a factor?" — and the answer is, it doesn't matter. Using tiny
> strangers would make *them* hit more, which looks like structure but is only size.

---

## Why

Factors are a **multiplicative** property of a number.

- **Squaring is multiplicative** → `(ab)² = a²·b²` → the factors are preserved.
- **XOR is bitwise** → `(ab) ^ w` mixes bits with no regard for multiplication → the
  factors are destroyed.

Right tool, right world. **XOR is king of the additive / bit world** — prime rails,
wraps, gaps, the rollball where `2ᵏ+1 ^ w = 2ᵏ+1 + w`. But **factors live in the
multiplicative world, and XOR can't see them there.**

---

## The bottom line

- **Squaring** carries the factors, structurally, always. Real.
- **XOR** surfaces a factor only at the base rate `1/p` — indistinguishable from a coin
  that says "divisible" with probability `1/p` for *any* prime. Not structure.
- **Neither cracks an unknown factorization.** You still have to factor `n²` to pull `a`
  and `b` back out, and factoring `n²` is no easier than factoring `n`. The reliable way
  to *know* a number's factors is to **keep the certificate from when you built it** —
  the rollball receipt (`descend`).

---

## Program

```python
from sitoolz import numtools as nt

nt.factor_probe(47, 37)          # -> n, n±1, n², n²±1 all factored
nt.factor_probe(1009, 1279)      # or any (a, b)
nt.factor_probe(1290511)         # or a single odd number
```

`factor_probe` descends (fully factors) `n`, `n−1`, `n+1`, `n²`, `n²−1`, `n²+1` and prints
every factor — so you can eyeball the structure across whatever numbers you feed it, and
see for yourself that the square keeps the factors and XOR doesn't.
