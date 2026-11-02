# How the periodicity falls out

*Derive the right edge of rule 30's period spectrum from one running-sum identity, and turn Rowland's eight measured terms into a counting function.*

**Lars Rocha**

---

## The gap, stated exactly

Rowland's 2006 paper reaches the right edge of rule 30. Section 3 is titled *"The
right side of rule 30,"* and it gets there by reflecting into rule 86 *"for
simplicity of indices."* He has more than the object — he has the periodicity:

> **Lemma 2.** *Let `f` be a right bijective, `k`-color, range `[−d,0]` rule... if
> columns `m−d, …, m−1` are periodic with lengths `p_{m−d}, …, p_{m−1}`, then column
> `m` is periodic with period length **dividing** `l(k) · lcm(p_{m−d}, …, p_{m−1})`.*

And it is immediate, not eventual — he says so directly:

> *"A result of Jen guarantees the eventual periodicity of columns in any range
> `[−d,0]` cellular automaton with a rightful initial condition. **For right
> bijective rules, the periodicity is not eventual but immediate.**"*

So the periodicity on this edge is his. For `k = 2`, `l(2) = 2`, and the periods
are powers of two.

**What Lemma 2 does not give is a criterion for when the period doubles.** It
gives a divisor — an upper bound. Whether column `m` actually doubles or holds
steady is outside it.

On the *left* side that job is done by a second theorem, Proposition 2, keyed on
Wolfram's observation:

> *"Each period doubling occurs exactly when a diagonal eventually becomes a white
> stripe, and the diagonal to its left has an odd number of black cells in each
> repeating block."*

**And on the right edge that criterion has no instances.** Measured over 520
generations from a single cell: the left diagonals go white at `m = 2` and `m = 7`;
**no right diagonal is ever white**, and bit 0 is constantly 1. The white-stripe
condition is not false there — it never occurs.

So Section 3 is left with periodicity and no handle on the spectrum. It reports

> *"`ι_I(t)` begins 1, 3, 1, 4, 1, 3, 1, 6, 1, 3, 1, 4, 1, 3, 1, 7, … and **has the
> structure** `ι_I(t) = a(ord₂(t+1))` for the sequence `{a(n)}` beginning
> **1, 3, 4, 6, 7, 9, 15, 16, …**"*

— and stops. The section contains no theorem, lemma, or proposition; its only
labelled item is Conjecture 1. `a(n)` stays a list of measurements. Seventeen
years later OEIS A094603 still carries the statement as a conjecture, contributed
2023-02-17.

**That is the gap: periodicity without a spectrum.** What follows closes it.

---

## The move: one identity

Put the origin on the cone boundary and let it travel. Outside the cone every
cell is 0 forever, so nothing outside ever enters the update, and the edge closes
on itself:

```
c_b(t+1) = c_b(t) XOR ( c_{b-1}(t) OR c_{b-2}(t) )        c_b = 0 for b < 0
c_b(0)   = [b = 0]
```

Two facts about this frame, and everything else is a consequence.

**It is triangular.** Bit `b` reads only `b−1` and `b−2`. Information flows
strictly upward, never down.

**It is autonomous.** Truncate at any width and the answer below the cut is
exact. Widths 8, 16, 32, 64 all produce `1, 3, 1, 4, 1, 3, 1, 6, 1, 3, 1, 4, 1, 3,
1, 7` — identical, because there is no history off to the side to lose. A column
of a two-dimensional pattern has no such property; that is what "eventual" is
about.

Now the identity. Since `c_b(0) = 0` for every `b ≥ 1` and the update is an XOR,
**each bit is the running XOR-sum of its own driver.** Writing
`d_b := c_{b-1} OR c_{b-2}`,

```
c_b(t) = XOR_{h=0}^{t-1} d_b(h)                                          (1)
```

That is the whole move. A question about how a nonlinear automaton evolves has
become a question about **the parity of a weight**.

---

## The periodicity, and the criterion, from the same line

Induct on `b`. Suppose every bit below `b` is periodic with a power-of-two period
and that `p` is nondecreasing there. Then `d_b` is an OR of two signals whose
periods are `2^{p(b-1)}` and `2^{p(b-2)}`, and an lcm of powers of two is the
larger, so

```
Q := 2^{p(b-1)}   is a period of d_b.
```

Let `S` be the weight of `d_b` over one period `Q`, taken mod 2. By (1), advancing
`t` by `Q` adds exactly `S` to the running sum. So:

```
S = 0  ->  c_b has period Q         (holds)
S = 1  ->  c_b has period 2Q        (doubles)
```

**Read that again, because it is the point.** The same computation that proves the
period is a power of two also says *which* power. Periodicity and the doubling
criterion are not two theorems here. They are one line, because the running-sum
identity makes them the same question — is a weight even or odd?

Hence, for `b ≥ 2`:

> **C1.** `c_b` is periodic with period **exactly** `2^{p(b)}`, and
> `p(b) ∈ { p(b−1), p(b−1) + 1 }`. In particular `p` is **nondecreasing**.
>
> **C2 (the criterion).** `b` is a doubling bit **exactly when `d_b` has odd
> weight over one period.**

Not "divides." **Exactly**, plus the step structure. That monotonicity is what
Lemma 2 cannot see, and it is the load-bearing fact for everything after.

Measured, 24 bits, 32,768 generations:

```
  b   p(b-1)    Q   driver wt   parity   doubles?   match
  2        1     2           2        0    False     True
  3        1     2           1        1     True     True
  4        2     4           3        1     True     True
  5        3     8           6        0    False     True
  6        3     8           5        1     True     True
  7        4    16          11        1     True     True
  8        5    32          24        0    False     True
  9        5    32          21        1     True     True
 10        6    64          48        0    False     True
 11        6    64          48        0    False     True
```

Agreement over every bit: **complete**. And `p = 0,1,1,2,3,3,4,5,5,6,6,6,…`
nondecreasing, stepping by 0 or 1 — verified.

---

## What the criterion unlocks

With `p` in hand, two lemmas finish it, and both turn on parity again.

**In phase (C3).** *If `v₂(t+1) ≥ p(b)` then `c_b(t) = 1`.*

Write `t = k·2^P − 1` and split (1) into whole driver periods plus one leftover
term. **The whole-period part vanishes either way** — a doubling bit contributes an
*even number* of odd-weight periods; a non-doubling bit contributes even-weight
periods. So `c_b(t)` collapses to a single evaluation `d_b(Q−1)`, where
`v₂(Q) = p(b−1)`, and the induction hypothesis applies there and returns 1.

And here is why the nonlinearity never obstructs the argument: **the OR is only
ever needed to be at least `c_{b-1}`. It is bounded below, never inverted.** A
proof that had to invert the OR would fail; this one does not.

**Out of phase (C4).** *At the first bit `B` with `p(B) > v`, every `t` with
`v₂(t+1) = v` exactly gives `c_B(t) = 0`.*

Because a doubling bit satisfies `c_b(t + 2^{P-1}) = NOT c_b(t)` for **every** `t` —
that window is exactly one odd-weight period, so (1) flips. Step forward `2^v`,
apply C3 to get 1, invert.

---

## The theorem

`p` nondecreasing means `{ b : p(b) ≤ v }` is an **initial segment** `{0, …, B−1}`.
C3 fills every one of them with 1. C4 puts a 0 at `B`. So the leading run has
length exactly `B`:

> ```
> run(t) = #{ b : p(b) ≤ v₂(t+1) }
> ```

It depends on `t` only through `v₂(t+1)`. ∎

Which is Rowland's `ι_I(t) = a(ord₂(t+1))` — now with `a` identified:

> **`a(v) = #{ b : p(b) ≤ v }` — the cumulative period spectrum.**

The table stops being a list of measurements and becomes a counting function over
the spectrum. `1, 3, 4, 6, 7, 9, 15, 16` are his eight terms, and they agree
exactly; the difference is that they are now *derived*, and the derivation keeps
going — rung 30 gives `a(30) = 72`.

**And because it is a counting function, the sequences interlock.** `v` occurs
`a(v) − a(v−1)` times in A094606, so A094604 and A094606 determine each other.
A094604's 70 published terms pin down 137 of A094606 and 90 of A094605 — verified
against every published term, zero mismatches.

---

## Why the frame was the whole thing

A reflection moves *which side* has bijectivity. It never creates it.

```
rule 30  bijective in [left]     ->  right side NOT bijective
rule 86  bijective in [right]    ->  left  side NOT bijective
```

Rowland reflected into rule 86 to put the bijectivity where his indices wanted
it — and that is a real and useful move; it is how Lemma 2 gets to apply and how
the periodicity becomes immediate rather than eventual.

The co-moving frame does something a reflection cannot:

```
E'[j] = E[j] XOR (E[j-1] OR E[j-2])
    bijective in its OWN coordinate:  yes
    reads only strictly lower coords:  yes
```

**Both properties, at once, on the side that had neither.** And it is exactly
those two that produce identity (1) — autonomy makes the sum closed, triangularity
makes the induction well-founded. Without (1) there is no parity question to ask,
and the doubling criterion has to come from somewhere else: from inspecting
diagonals for whiteness, which on this edge never happens.

That is the difference between having the periodicity and having the spectrum.
Rowland has an upper bound from bijectivity and a criterion that lives on the
other side. The co-moving frame turns the criterion into the parity of a computable
weight — and the same induction that yields the periodicity yields the criterion,
because in this frame they are one question.

---

## Verification

Unconditional proof; every step checked before it was written
(`sitoolz.edgeproof`, `sitoolz.edgelaw`, `sitoolz.tri_xor`):

| step | checked | failures |
|---|---|---|
| C1 — `p` nondecreasing, steps in {0,1} | 30 bits | 0 |
| C2 — driver weight odd iff doubling | 29 bits | 0 |
| the bracket vanishes, `c_b(t) = d_b(k·2^P−1)` | 145 cases | 0 |
| `v₂(Q) = p(b−1)` and `c_{b-1}(Q−1) = 1` | 29 cases | 0 |
| C3 (in phase) | 181,072 cases | 0 |
| antipodal lemma | 12 doubling bits | 0 |
| C4 (out of phase) | 65,520 cases | 0 |
| law `run(g) = a(v₂(g+1))` | to generation 4096 | 0 |
| against real rule 30's rightmost run | 500 generations | 0 |
| table vs Rowland's 8 published terms | 8 terms | 0 |

Constant-time evaluation at any generation: `edge(g) = a(v₂(g+1))`, no simulation.
The triangular structure also bounds memory — truncation is exact below the cut —
so the edge of a forever-growing cone is trackable in a fixed 128-bit register at
any generation whatsoever; `sitoolz.edgereach` sustains ~900 million generations per
second on three 64-bit words.

## Attribution

The object, the law's statement, and the first eight terms of `a(n)` are in
Rowland, *Local nested structure in rule 30*, Complex Systems **16** (2006)
239–258, Section 3. The immediate periodicity of columns for right-bijective
rules is his Lemma 2. Wolfram's doubling observation is NKS p. 871, proved as his
Proposition 2 for the left side.

What is proved here: the doubling criterion on the right edge, the exact period
and its step structure, and the identification of `a` as the cumulative period
spectrum — hence a proof of the law A094603 has carried as a conjecture since
2023-02-17.
