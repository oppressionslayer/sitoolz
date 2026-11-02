# The bits of powers

### Read the bits a square cannot have, and use them as a gate

**Lars Rocha**

---

## For anyone reading

Write out some squares in binary and read the 2s column: it is zero every time. That one
fact decides exactly when XOR is addition, and it generalises to every power.

```
1^2 =  1 = 00000001
2^2 =  4 = 00000100
3^2 =  9 = 00001001
4^2 = 16 = 00010000
5^2 = 25 = 00011001
6^2 = 36 = 00100100
7^2 = 49 = 00110001
8^2 = 64 = 01000000
```

Look at the second column from the right — the **2s** place. It is zero every
time. Not usually. **Always**, for every square there has ever been.

The reason is short. Every whole number is either even or odd. An even number
is `2m`, and `(2m)^2 = 4m^2`, which is a multiple of 4. An odd number is
`2m+1`, and `(2m+1)^2 = 4m^2 + 4m + 1`, which is a multiple of 4 plus one. So a
square is always **0 or 1 more than a multiple of 4** — and neither of those
ever lights up the 2s bit.

That one fact is the seed of everything below.

## Why it matters: XOR is addition that forgot to carry

`XOR` compares two numbers bit by bit and keeps a 1 where they differ. It looks
nothing like addition, and then sometimes it is:

```
5 XOR 2 = 7     and     5 + 2 = 7      same
6 XOR 5 = 3     but     6 + 5 = 11     different
```

The rule is exact: **`a XOR b` equals `a + b` precisely when `a` and `b` share
no 1-bits.** In binary, `5 = 101` and `2 = 010` never collide, so nothing has to
carry and XOR *is* the sum. But `6 = 110` and `5 = 101` both own the 4s bit;
that collision needs a carry, XOR throws the carry away, and the answer comes
out wrong.

So "when does XOR act like plus?" is really "when do the bits stay out of each
other's way?" — and for squares, we already know one place they always do.

## The question

Take a square `s` and a mask `m`, and ask how often

```
s XOR m == s + m
```

Since that happens exactly when `s AND m == 0`, it is a question about which
bits of squares are ever set. Squares are rigid down there, so the answers come
out as clean halves and quarters instead of messy decimals:

```
mask  1 (      1) ->  1/2 of all squares
mask  2 (     10) ->  ALL of them        <- the 2s bit is never set
mask  5 (    101) ->  1/4
mask 12 (   1100) ->  1/2
mask 21 (  10101) ->  1/8
mask 255(11111111) -> 1/16
```

No approximation, no "about." Exact fractions, every time.

---

## The two results

`ntpar` answered this one bit at a time:

```
bit 0    exactly 1/2, forever
bit 1    exactly 0, forever
bit b    1/2 + 2^-(b//2+1)     -- creeps toward a half, never arrives
```

Two questions follow: whether you can just multiply those together for a
whole mask, and what happens for cubes and higher powers. Both have exact
answers, and the first one is a **no**.

### 1. The bits of squares are not independent

The obvious guess is that the chance of two bits both being clear is the product
of their separate chances. It is not, and the error is never messy:

| bits | truth | independent guess | ratio |
|---|---|---|---|
| {0,2} | 1/4 | 3/8 | **2/3** |
| {2,3} | 1/2 | 9/16 | **8/9** |
| {2,4} | 3/8 | 15/32 | **4/5** |
| {3,4} | 1/2 | 15/32 | **16/15** |
| {0,5} | 3/8 | 5/16 | **6/5** |
| {0,2,4} | 1/8 | 15/64 | **8/15** |

Some too high, some too low, all exactly rational. So the bits of squares
*talk to each other*, and by a definite amount.

**Bit 1 is the one exception.** Paired with anything, its ratio is exactly 1 —
because it is always clear, so it has nothing to correlate with.

### 2. But the exact answer is simpler than the correction

You do not need the correction at all. `k mod 2^n` completely determines
`k^2 mod 2^n`, and a mask `n` bits wide cannot see any further. So the whole
question is a finite count over `2^n` cases — no sampling, no limit, and the
answer is always a dyadic fraction.

```python
from fractions import Fraction

def joint_law(m, power=2):
    """Exact fraction of k**power with (k**power) & m == 0."""
    if m == 0:
        return Fraction(1)
    n = m.bit_length()
    M = 1 << n
    good = sum(1 for k in range(M) if (pow(k, power, M) & m) == 0)
    return Fraction(good, M)
```

Six lines, exact for every mask. **Verified against 200,000 squares across 64
masks: 0 mismatches.**

### 3. Squares are not special — and they are the only irregular case

The same machinery runs for any power, and the "forced" bits — the ones that
are clear for *every* value — follow an exact rule.

> **Theorem.** Write `p = 2^e · q` with `q` odd. Then `k^p mod 2^(e+2)` is 0 or
> 1 for every `k`, so bits `1 … e+1` are always clear — **provided `p ≥ e+2`.**
> The only exponent where that fails is `p = 2`, which forces bit 1 alone.
>
> **Odd powers force nothing**, because `x → x^q` is a bijection on the units
> mod `2^n` when `q` is odd, so the residues stay uniform.

| power | forced bits | why |
|---|---|---|
| 2 | 1 | `k²` is 0 or 1 mod 4 |
| 4 | 1,2,3 | `k⁴` is 0 or 1 mod 16 |
| 6 | 1,2 | `k⁶` is 0 or 1 mod 8 |
| 8 | 1,2,3,4 | `k⁸` is 0 or 1 mod 32 |
| 16 | 1,2,3,4,5 | `k¹⁶` is 0 or 1 mod 64 |
| any odd | none | uniform |

**Verified over p = 1..40, 0 mismatches.**

**The proof in a sentence.** For odd `k`, `k^q` is odd, an odd square is 1 mod
8, and squaring `e` times pushes that to `k^p ≡ 1 mod 2^(e+2)`. For even `k`,
`k^p ≡ 0 mod 2^p`, which covers `2^(e+2)` exactly when `p ≥ e+2`.

**And that is why squares are the odd one out.** `p = 2` gives `e = 1`, so the
rule wants `p ≥ 3` — and 2 is not 3. The residues mod 8 come out `{0, 1, 4}`
instead of `{0, 1}`, that stray 4 carries a 2s-place bit at position 2, and bit
2 escapes. Every other exponent in the family collapses cleanly. The most famous
case is the single broken one.

---

## So what does that say about every other power?

### Powers are pinned, not scattered

The rule is not a curiosity about small numbers. Every power collapses onto a
tiny set of possible binary endings, and the higher the power, the tinier it
gets:

| power | how many values it can be | which |
|---|---|---|
| every 2nd | **3** values mod 8 | 0, 1, 4 |
| every 4th | **2** values mod 16 | 0, 1 |
| every 8th | **2** values mod 32 | 0, 1 |
| every 16th | **2** values mod 64 | 0, 1 |
| every 32nd | **2** values mod 128 | 0, 1 |

Take any whole number at all — 7, or 1000, or a googol. Raise it to the 32nd
power. The answer is **0 or 1 mod 128**: two possibilities out of a hundred and
twenty-eight. Every 32nd power that will ever exist ends in one of two binary
tails.

### And XOR stops being exotic

A mask works *always* — `x XOR m == x + m` for every value, no exceptions —
exactly when every bit it uses is a forced-clear one. So they can be counted:

| power | masks where XOR is always addition | count |
|---|---|---|
| squares | 2 | **1** |
| 4th powers | 2, 4, 6, 8, 10, 12, 14 | 7 |
| 8th powers | 2, 4, 6, 8, 10, … | 15 |
| 16th powers | 2, 4, 6, 8, 10, … | 31 |
| 32nd powers | 2, 4, 6, 8, 10, … | **63** |

Always `2^(e+1) - 1`. For squares there is exactly **one** number you can XOR
and be guaranteed addition. For 32nd powers there are **sixty-three**. The
higher the power, the more of ordinary arithmetic collapses into plain
bit-flipping.

### The part that surprised me

The odd half of the exponent does nothing at all.

| these behave identically | forced bits |
|---|---|
| `k^4, k^12, k^20, k^36 …` | 1, 2, 3 |
| `k^6, k^10, k^14, k^22 …` | 1, 2 |
| `k^3, k^5, k^7, k^9, k^11 …` | none |

`k^12` behaves exactly like `k^4`. That extra factor of three is invisible.
Only `v2(p)` — how many times two divides the exponent — has any effect.
Cubing, quintupling, raising to the 999th power: none of it registers.

> **So the rigidity of a power has nothing to do with how big the exponent is.**
> `k^1024` and `k^1025` live in different worlds — one is pinned to two values
> mod 4096, the other is perfectly uniform. Meanwhile `k^1024` and `k^3072` are
> identical twins.

The whole infinite family of powers collapses into a ladder indexed by a single
number. And squares — the rung everybody studies — are the one step where the
ladder is broken.

---

## The whole chain

```
5 XOR 2 = 7                       noticed in the 1980s, staring at hex
   |
a XOR b == a + b  iff  a AND b == 0        XOR is addition minus the carry
   |
squares are 0 or 1 mod 4                   so the 2s bit is never set
   |
square_xor's rates are dyadic              measured
   |
ntpar's per-bit law                        1/2, 0, then 1/2 + 2^-(b//2+1)
   |
joint_law -- exact, independence is false  derived
   |
forced bits determined by v2(p)            every power, one rule
```

Each step is the previous one asked more carefully.

---

## Run it

```python
from sitoolz import ntpar_joint as nj

nj.joint_law(5)             # Fraction(1, 4) -- exact, no counting
nj.mask_frame(24)           # every mask: truth, guess, and the error
nj.independence_frame(7)    # pairwise dependence between bits
nj.add_sub_frame()          # reproduces numtools.square_xor exactly
nj.power_frame(16)          # forced bits per power, with residues
nj.verify_forced()          # {'powers': 40, 'mismatches': 0}
nj.verify()                 # {'masks': 64, 'mismatches': 0}
```

Stdlib plus `dfkit`, so the frames degrade instead of raising when pandas is
absent.

---

## Provenance

`sitoolz.ntpar` — the single-bit law, and the observation that started this.
`numtools.square_xor` — the measured dyadic rates, and the note that *"squares
are rigid in the low bits: every square is 0 or 1 mod 4, so bit 1 is NEVER
set."* That sentence is the whole thing in advance.
`sitoolz.ntpar_joint` — the exact joint law, the dependence table, and the
forced-bit theorem.

Both results are proved, with their verification counts alongside; nothing here
rests on a sample.
