# A proof of the v₂ law for rule 30's light-cone edge

*Prove the formula OEIS A094603 has carried as a conjecture since 2023: `a(n) = A094604(v₂(n))`, and read the edge run at any generation from its 2-adic valuation.*

**Lars Rocha**

---

## Statement

Let rule 30 run from a single black cell. Let `run(t)` be the length of the
maximal run of rightmost black cells in row `t`. Let `v₂` denote the 2-adic
valuation.

> **Theorem.** `run(t)` depends only on `v₂(t+1)`.
>
> Precisely, `run(t) = #{ b : p(b) ≤ v₂(t+1) }`, where `p(b)` is defined below.

Equivalently, in OEIS terms, `A094603(n) = A094604(A007814(n))` — the formula
carried as a conjecture in A094603 since 2023-02-17.

---

## Setup

In the frame co-moving with the light cone the edge is autonomous: cells
outside the cone are 0 forever, so they never enter the update. Writing
`c_b(t)` for bit `b` of the edge at generation `t`,

```
c_b(t+1) = c_b(t) XOR ( c_{b-1}(t) OR c_{b-2}(t) ),      c_b = 0 for b < 0
c_b(0)   = [b = 0]
```

and `run(t)` is the number of leading 1s of `(c_0(t), c_1(t), …)`.

Write `d_b := c_{b-1} OR c_{b-2}` for the **driver** of bit `b`. Since
`c_b(0) = 0` for `b ≥ 1` and the update is an XOR, bit `b` is the running
XOR-sum of its driver:

```
c_b(t) = XOR_{h=0}^{t-1} d_b(h)                                    (1)
```

**Lemma 0 (C1).** Every `c_b` is periodic in `t` with period exactly `2^{p(b)}`
for some `p(b) ≥ 0`; `p(0) = 0`, `p(1) = 1`; and for `b ≥ 2`,
`p(b) = p(b-1)` or `p(b-1) + 1`. In particular **`p` is nondecreasing**.

*Proof.* Induction on `b`. `c_0(t+1) = c_0(t)`, so `c_0 ≡ 1` and `p(0) = 0`.
`c_1(t+1) = c_1(t) XOR c_0(t) = c_1(t) XOR 1`, so `c_1(t) = t mod 2` and
`p(1) = 1`. For `b ≥ 2`, assume the claim below `b`; then `p` is nondecreasing
there, so `d_b` — an OR of signals of periods `2^{p(b-1)}` and `2^{p(b-2)}` —
has period `Q := 2^{p(b-1)}` (an lcm of powers of two is the larger). Let `S`
be the weight of `d_b` over `Q`, taken mod 2. By (1), advancing `t` by `Q` adds
`S` to the running sum, so `c_b` has period `Q` if `S = 0` and `2Q` if
`S = 1`. Either way a power of two, and `p(b) ∈ {p(b-1), p(b-1)+1}`. ∎

Call `b` a **doubling bit** when `p(b) = p(b-1) + 1`. Lemma 0's proof gives the
dichotomy we will use twice:

> **(D)** `d_b` has **odd** weight over `Q = 2^{p(b-1)}` exactly when `b` is a
> doubling bit, and **even** weight otherwise.

---

## The phase lemma

> **Lemma 1 (C3).** For every bit `b` and every `t` with `v₂(t+1) ≥ p(b)`,
> `c_b(t) = 1`.

*Proof.* Strong induction on `b`.

*Base.* `c_0 ≡ 1`, so the claim is immediate. For `b = 1`: `p(1) = 1`, and
`v₂(t+1) ≥ 1` forces `t` odd, where `c_1(t) = t mod 2 = 1`.

*Step.* Let `b ≥ 2`, put `P := p(b)`, and assume the lemma for all indices
below `b`. Let `t` satisfy `v₂(t+1) ≥ P`, so `t = k·2^P − 1` for some `k ≥ 1`.
Put `Q := 2^{p(b-1)}`, a period of `d_b`. By (1),

```
c_b(t) = XOR_{h=0}^{k·2^P − 2} d_b(h)
       = [ XOR_{h=0}^{k·2^P − 1} d_b(h) ] XOR d_b(k·2^P − 1).           (2)
```

**The bracket vanishes, in both cases.**

- If `b` is a doubling bit then `P = p(b-1) + 1`, so `Q = 2^{P-1}` and
  `k·2^P = 2k·Q` is an **even** number of full periods of `d_b`. By (D) each
  period has odd weight, and `2k` of them have weight `≡ 0 (mod 2)`.
- If `b` is not a doubling bit then `P = p(b-1)`, so `Q = 2^P` and `k·2^P`
  is `k` full periods. By (D) each has even weight, so the total is even.

Hence `c_b(t) = d_b(k·2^P − 1)`.

**Reduce the argument.** In both cases `Q` divides `2^P`, so
`k·2^P ≡ 0 (mod Q)` and therefore `k·2^P − 1 ≡ Q − 1 (mod Q)`. As `Q` is a
period of `d_b`,

```
c_b(t) = d_b(Q − 1).                                                   (3)
```

**Apply the hypothesis.** `v₂((Q−1)+1) = v₂(Q) = log₂ Q = p(b-1)`, so the
induction hypothesis at `b-1` applies at time `Q−1` and gives
`c_{b-1}(Q−1) = 1`. Therefore

```
d_b(Q − 1) = c_{b-1}(Q−1) OR c_{b-2}(Q−1) = 1,
```

and by (3), `c_b(t) = 1`. ∎

Note what the OR does here: it is needed only to be *at least* `c_{b-1}`. The
nonlinearity is never inverted, only bounded below — which is why the induction
survives it.

---

## The stopping lemma

> **Lemma 2 (antipodal).** If `b` is a doubling bit with `p(b) = P`, then
> `c_b(t + 2^{P-1}) = NOT c_b(t)` for **every** `t`.

*Proof.* By (1), `c_b(t + 2^{P-1}) XOR c_b(t)` is the XOR-sum of `d_b` over a
window of length `2^{P-1} = Q`, i.e. exactly one period of `d_b`. By (D) that
weight is odd, so the difference is 1 for every `t`. ∎

> **Lemma 3 (C4).** Fix `v ≥ 0` and let `B := min{ b : p(b) > v }`. Then for
> every `t` with `v₂(t+1) = v` **exactly**, `c_B(t) = 0`.

*Proof.* `p` steps by 0 or 1 (Lemma 0) and `p(B) > v ≥ p(B-1)`, so
`p(B) = v + 1` and `B` is a doubling bit. Write `t + 1 = m·2^v` with `m` odd.
Then `(t + 2^v) + 1 = (m+1)·2^v` with `m + 1` even, so
`v₂(t + 2^v + 1) ≥ v + 1 = p(B)`. Lemma 1 gives `c_B(t + 2^v) = 1`. Lemma 2
with `P = v + 1` gives `c_B(t + 2^v) = NOT c_B(t)`. Hence `c_B(t) = 0`. ∎

---

## The theorem

*Proof.* Fix `t` and let `v = v₂(t+1)`. Since `p` is nondecreasing,
`{ b : p(b) ≤ v }` is an initial segment `{0, 1, …, B−1}` where
`B = min{ b : p(b) > v }`. Lemma 1 gives `c_b(t) = 1` for every `b < B`, and
Lemma 3 gives `c_B(t) = 0`. So the leading run has length exactly `B`, that is

```
run(t) = #{ b : p(b) ≤ v₂(t+1) },
```

which depends on `t` only through `v₂(t+1)`. ∎

**Corollary.** `A094604[v] = #{ b : p(b) ≤ v }` — the table is the cumulative
period spectrum, not a table of measurements. Consequently `v` occurs
`A094604[v] − A094604[v−1]` times in A094606, so the two sequences determine
each other.

---

## Verification

The proof is unconditional, but every step was checked before it was written,
over 30 bits and 65,536 generations (`sitoolz.edgeproof`, and the
step-by-step battery in `sitoolz.parity_lab`):

| step | checked | failures |
|---|---|---|
| Lemma 0 — `p` nondecreasing, steps in {0,1} | 30 bits | 0 |
| (D) — driver weight odd on doubling bits, even otherwise | 29 bits | 0 |
| (2) — the bracket vanishes, `c_b(t) = d_b(k·2^P−1)` | 145 cases | 0 |
| (3) — `v₂(Q) = p(b-1)` and `c_{b-1}(Q−1) = 1` | 29 cases | 0 |
| Lemma 1 (C3) | 181,072 cases | 0 |
| Lemma 2 (antipodal) | 12 doubling bits | 0 |
| Lemma 3 (C4) | 65,520 cases | 0 |

---

## What was already known

**Wolfram**, *A New Kind of Science* (2002), p. 871, observes the period
doubling of the diagonals on the left side of rule 30. **Rowland**, *Local
nested structure in rule 30*, Complex Systems **16** (2006) 239–258, proves
that left-side doubling criterion (his Proposition 2) and, in Lemma 2, gives
immediate periodicity with a divisibility bound on column periods for
right-bijective rules; his Section 3 records the `ord₂(t+1)` structure of
the right edge as an observation with the first eight terms of `a(n)`.
**OEIS A094603–A094606** record the sequences, since 2004.

The `v₂` statement is carried in A094603 as a conjecture, contributed by
Alan Michael Gómez Calderón on 2023-02-17.

Proposition 2 is about the **left** side of rule 30; this note is about the
**right** light-cone edge, where Lemma 2 gives the periodicity and the
doubling criterion (D), the phase lemmas and the law are proved here from
the recurrence, self-contained and citing nothing. The comparison with
Rowland's route is in [EDGE_LAW.md](EDGE_LAW.md).
