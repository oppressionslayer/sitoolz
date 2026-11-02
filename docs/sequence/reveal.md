# The reveal

XOR each root with 5 and every perfect square from 0² up falls out in order. This is the
answer for anyone who tried the sequence and wants the mechanism, top to bottom: the move
first, then the proof it earns.

---

## The terms

```
1, 16, 25, 144, 441, 784, 1089, 2704, 4761, 7056, 9409, 15376, …
```

All perfect squares. Their roots:

```
1, 4, 5, 12, 21, 28, 33, 52, 69, 84, 97, 124, …
```

Those roots look like nothing. Differences are 3, 1, 7, 9, 7, 5, 19, 17, 15,
13, 27 — no pattern, no growth law you'd recognize.

---

## The move

XOR each root with 5.

```
  1 ⊕ 5 =   4        33 ⊕ 5 =  36
  4 ⊕ 5 =   1        52 ⊕ 5 =  49
  5 ⊕ 5 =   0        69 ⊕ 5 =  64
 12 ⊕ 5 =   9        84 ⊕ 5 =  81
 21 ⊕ 5 =  16        97 ⊕ 5 = 100
 28 ⊕ 5 =  25       124 ⊕ 5 = 121
```

Sorted, that is

```
0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121
```

**Every square from 0² to 11². No gaps, nothing extra, nothing missing.**

That's the whole thing.

---

## Why it was hard to see

`^` in the original code is **bitwise XOR**, not exponentiation. Almost
everyone reads `x^5` as "x to the fifth" and never recovers. The operation
is hiding in plain sight as a typo.

And there is no lookup route — the sequence isn't in the OEIS. Five searches
confirm it: the full list, the first eight terms, the x-form, the generation
ordering, and the tail. All `null`. So the usual first move fails silently.

---

## Why it works

XOR by a constant is an **involution**: apply it twice and you're back where
you started. So `n ↦ n ⊕ 5` is a *bijection* on the non-negative integers.

That means the condition "n ⊕ 5 is a square" isn't a filter at all. Solve it:

$$n \oplus 5 = m^2 \quad\Longleftrightarrow\quad n = m^2 \oplus 5$$

Every `m` gives exactly one `n`. Nothing is being searched or rejected. The
original loop scanned 16,430 candidates and kept 12 — which *looks* sparse,
about 0.07% — but the sparseness is just a fourth root:

$$\bigl|\{x < X\}\bigr| \;=\; X^{1/4}$$

| X | count | X^(1/4) |
|---|---|---|
| 10⁸ | 100 | 100.00 |
| 10¹² | 1,000 | 1,000.0 |
| 10²⁰ | 100,000 | 100,000 |
| 10¹⁰⁰ | **10²⁵** | 10²⁵ |

Exactly, with no error term. And you can name any term you like. The
10²⁵-th root is

```
100000000000000000000000000000000000000000000000005
```

computed in microseconds, because it's just `(10²⁵)² ⊕ 5`.

---

## The three offsets, and the missing fourth

Since `5 = 101₂`, XOR-ing with it flips only bits 0 and 2. So `n − m²` can
only take four values:

| bit 2 | bit 0 | change |
|---|---|---|
| 0 | 0 | +5 |
| 0 | 1 | +3 |
| 1 | 0 | −3 |
| 1 | 1 | −5 |

But only **three** ever occur:

```
m mod 4:   0 → +5,   1 → +3,   2 → −3,   3 → +3
```

The `−5` case never happens, and that absence is the prettiest part.
A perfect square is always `0`, `1`, or `4` mod 8 — in binary `000`, `001`,
`100`. **Bits 0 and 2 are never both set in a square.** So the fourth case
is arithmetically impossible.

---

## Bonus: this is Kummer's theorem

XOR agrees with ordinary addition exactly when adding produces **no carries**.
And Kummer's theorem says the number of base-2 carries in `a + b` equals the
exponent of 2 dividing `C(a+b, b)`.

So `m² ⊕ 5 = m² + 5` — the `+5` branch — happens precisely when
`C(m²+5, 5)` is **odd**. Verified for m up to 3000. Equivalently by Lucas'
theorem: `5`'s bits are a subset of `m²+5`'s bits.

The bit-flipping story and the binomial-coefficient story are the same story.

---

## The sting in the tail: try to index them

Given the terms, assign each one its index `m`. The obvious answer — sorted
position — gets **11 of 13 right**.

```
   x       sorted position     true m
   1              0              2      ← wrong
  16              1              1      ← correct, by luck
  25              2              0      ← wrong
 144              3              3
 441              4              4
   …              …              …
```

It misses exactly **two**, and they are **swapped** — positions 0 and 2 trade
places while position 1 sits still. Because the map is V-shaped:

```
m:  0    1    2    3     4     5
x:  25   16   1    144   441   784
```

It descends to a minimum at `m = 2`, then climbs forever after. Sorted order
and index order agree everywhere except across that dip, and nothing in the
terms tells you where the dip is. You need `m = isqrt(isqrt(x) ⊕ 5)`, which
needs the 5.

So even after you know it's a bijection, you can still get the indexing wrong
— you also have to know the bijection isn't monotonic.

---

## Why it isn't in the OEIS

Not an oversight. For any fixed `k`, the correction `(m² ⊕ k) − m²` is
periodic in `m` with period `2^(bitlen(k)−1)` — because a fixed `k` touches
only bounded low bits, and `m² mod 2^j` is periodic in `m`.

So every sequence in this family is a **quasi-polynomial**: degree 2, and a
finite table of offsets. Once you have the period there is nothing left to
find. Every term, forever, from four cases.

Degenerate case worth knowing: `k = 2` has period **1**. Since bit 1 is never
set in a square, `m² ⊕ 2 = m² + 2` for every `m` — XOR is *identically*
addition there. And `k = 2` is the only positive constant for which that
holds, because bit 1 is the only bit no square ever sets.

---

## Reproduce any of it

```python
from sitoolz import numtoolsd as ntd

ntd.xor_k_squares(5)           # the table
ntd.count_below(10**100)       # 10**25, instantly
ix = ntd.xor_k_index(5)
ix.unrank(10**25)              # the 51-digit root
ntd.xor_k_period(5)            # (4, [5, 3, -3, 3])
ntd.kummer_check(16, 5)        # carries / binomial / XOR, all agreeing
```

And `rulefinder` (`from sitoolz.sequence import rulefinder`) cracks it from bare
terms in 0.11 ms — by searching over *rules* rather than over terms:

```python
for k in range(1024):
    if all(is_square(root ^ k) for root in roots):
        return k
```

---

## The ladder's count

The session's `coun` column is the same object as `master_rail.count` in `numtools`.
In `master_rail`:

```
count = π(<x) − 1
counn = count − 1
railn = count − 2
```

A three-deep offset ladder off the prime-counting function: `count` sits one below
the number of primes under `x`, `counn` one below that, `railn` one below again. The
`coun` column of the session is identical to `master_rail.count`.
