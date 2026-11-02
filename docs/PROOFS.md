# Proofs — rule 30's light-cone edge, and triangular XOR-driven systems

*Every theorem behind the edge law, written out, with the proved kept separate from the measured.*

**Lars Rocha**

Sections 1–5 are proofs. Section 6 lists what is verified but not derived;
nothing there should be quoted as a theorem.

Notation throughout. The light-cone edge in the co-moving frame:

```
c_b(t+1) = c_b(t) XOR ( c_{b-1}(t) OR c_{b-2}(t) ),    c_b = 0 for b < 0
c_b(0)   = [b = 0]
```

`c_b(t)` is bit `b` at generation `t`; `d_b := c_{b-1} OR c_{b-2}` is the
**driver** of bit `b`; `run(t)` is the number of leading 1s; `v₂` is the 2-adic
valuation. The edge is autonomous — cells outside the cone are 0 forever and
never enter the update — so this is not an approximation of rule 30 but rule 30
restricted to a region nothing else can reach.

Since `c_b(0) = 0` for `b ≥ 1` and the update is an XOR, bit `b` is the running
XOR-sum of its driver:

```
c_b(t) = XOR_{h=0}^{t-1} d_b(h)                                          (1)
```

---

## 1. Periodicity

> **Lemma 0.** Every `c_b` is periodic in `t` with period exactly `2^p(b)`;
> `p(0) = 0`, `p(1) = 1`; and for `b ≥ 2`, `p(b) ∈ {p(b-1), p(b-1)+1}`. In
> particular **`p` is nondecreasing**.

*Proof.* Induction on `b`.

`c_0(t+1) = c_0(t) XOR (c_{-1} OR c_{-2}) = c_0(t)`, and `c_0(0) = 1`, so
`c_0 ≡ 1` and `p(0) = 0`. Next `c_1(t+1) = c_1(t) XOR (c_0(t) OR 0) =
c_1(t) XOR 1`, so `c_1(t) = t mod 2` and `p(1) = 1`.

For `b ≥ 2`, assume the claim below `b`; then `p` is nondecreasing there, so
`d_b` — an OR of signals of periods `2^p(b-1)` and `2^p(b-2)` — has period
`Q := 2^p(b-1)`, since an lcm of powers of two is the larger of them. Let `S`
be the weight of `d_b` over `Q`, taken mod 2. By (1), advancing `t` by `Q` adds
`S` to the running sum, so `c_b` has period `Q` when `S = 0` and `2Q` when
`S = 1`. Either way a power of two, and `p(b) ∈ {p(b-1), p(b-1)+1}`. ∎

Call `b` a **doubling bit** when `p(b) = p(b-1) + 1`. The proof gives:

> **(D)** `d_b` has **odd** weight over `Q = 2^p(b-1)` exactly when `b` is a
> doubling bit, and **even** weight otherwise.

The OR is the only nonlinearity in rule 30, and this argument routes around it:
the OR defines `d_b` and is never used again. What carries the induction is
that XOR-summation preserves power-of-two periodicity regardless of the summand.

---

## 2. The phase lemma

> **Lemma 1.** For every bit `b` and every `t` with `v₂(t+1) ≥ p(b)`,
> `c_b(t) = 1`.

*Proof.* Strong induction on `b`.

*Base.* `c_0 ≡ 1`, so the claim is immediate. For `b = 1`: `p(1) = 1`, and
`v₂(t+1) ≥ 1` forces `t` odd, where `c_1(t) = t mod 2 = 1`.

*Step.* Let `b ≥ 2`, put `P := p(b)`, and assume the lemma for all indices
below `b`. Let `t` satisfy `v₂(t+1) ≥ P`, so `t = k·2^P − 1` for some `k ≥ 1`.
Put `Q := 2^p(b-1)`, a period of `d_b`. By (1),

```
c_b(t) = XOR_{h=0}^{k·2^P − 2} d_b(h)
       = [ XOR_{h=0}^{k·2^P − 1} d_b(h) ] XOR d_b(k·2^P − 1).            (2)
```

**The bracket vanishes, in both cases.**

- If `b` is a doubling bit then `P = p(b-1) + 1`, so `Q = 2^{P-1}` and
  `k·2^P = 2k·Q` is an **even** number of full periods of `d_b`. By (D) each
  period has odd weight, and `2k` of them sum to `0 mod 2`.
- If `b` is not a doubling bit then `P = p(b-1)`, so `Q = 2^P` and `k·2^P` is
  `k` full periods. By (D) each has even weight, so the total is even.

Hence `c_b(t) = d_b(k·2^P − 1)`.

**Reduce the argument.** In both cases `Q` divides `2^P`, so
`k·2^P ≡ 0 (mod Q)` and therefore `k·2^P − 1 ≡ Q − 1 (mod Q)`. As `Q` is a
period of `d_b`,

```
c_b(t) = d_b(Q − 1).                                                     (3)
```

**Apply the hypothesis.** `v₂((Q−1)+1) = v₂(Q) = log₂ Q = p(b-1)`, so the
induction hypothesis at `b-1` applies at time `Q−1` and gives
`c_{b-1}(Q−1) = 1`. Therefore `d_b(Q−1) = c_{b-1}(Q−1) OR c_{b-2}(Q−1) = 1`,
and by (3), `c_b(t) = 1`. ∎

Note what the OR does here: it is needed only to be *at least* `c_{b-1}`. The
nonlinearity is never inverted, only bounded below — which is why the induction
survives it.

---

## 3. The antipodal lemma, and where the run stops

> **Lemma 2 (antipodal).** If `b` is a doubling bit with `p(b) = P`, then
> `c_b(t + 2^{P-1}) = NOT c_b(t)` for **every** `t`.

*Proof.* By (1), `c_b(t + 2^{P-1}) XOR c_b(t)` is the XOR-sum of `d_b` over a
window of length `2^{P-1} = Q`, i.e. exactly one period of `d_b`. By (D) that
weight is odd, so the difference is 1 for every `t`. ∎

> **Corollary 2.1.** A doubling bit has weight exactly `T/2` over one period
> `T = 2^{p(b)}`; hence its period-weight parity is 0 whenever `p(b) ≥ 2`.

*Proof.* By Lemma 2 the two halves of the period are complementary, so exactly
one of `{t, t + T/2}` carries a 1 for each `t < T/2`. The weight is `T/2`,
which is even for `p(b) ≥ 2`. ∎

> **Lemma 3.** Fix `v ≥ 0` and let `B := min{ b : p(b) > v }`. Then for every
> `t` with `v₂(t+1) = v` **exactly**, `c_B(t) = 0`.

*Proof.* `p` steps by 0 or 1 (Lemma 0) and `p(B) > v ≥ p(B-1)`, so
`p(B) = v + 1` and `B` is a doubling bit. Write `t + 1 = m·2^v` with `m` odd.
Then `(t + 2^v) + 1 = (m+1)·2^v` with `m + 1` even, so
`v₂(t + 2^v + 1) ≥ v + 1 = p(B)`. Lemma 1 gives `c_B(t + 2^v) = 1`. Lemma 2
with `P = v + 1` gives `c_B(t + 2^v) = NOT c_B(t)`. Hence `c_B(t) = 0`. ∎

---

## 4. The v₂ law

> **Theorem.** `run(t) = #{ b : p(b) ≤ v₂(t+1) }`. In particular `run(t)`
> depends on `t` only through `v₂(t+1)`.

*Proof.* Fix `t` and let `v = v₂(t+1)`. Since `p` is nondecreasing,
`{ b : p(b) ≤ v }` is an initial segment `{0, 1, …, B−1}` where
`B = min{ b : p(b) > v }`. Lemma 1 gives `c_b(t) = 1` for every `b < B`, and
Lemma 3 gives `c_B(t) = 0`. So the leading run has length exactly `B`. ∎

This is the statement carried in OEIS A094603 as a conjecture since
2023-02-17: `A094603(n) = A094604(A007814(n))`.

> **Corollary 4.1 (the cumulative spectrum).** `A094604[v] = #{b : p(b) ≤ v}`.
> Hence `v` occurs exactly `A094604[v] − A094604[v−1]` times in A094606, and
> the four sequences A094603–A094606 determine one another.

> **Corollary 4.2.** A094604 is **strictly** increasing.

*Proof.* `p` is nondecreasing, steps by 0 or 1 (Lemma 0), and is unbounded, so
every exponent `v` is attained by at least one bit. A repeated value in
A094604 would mean some period `2^v` is attained by none. ∎

---

## 5. Triangular XOR-driven systems

The class: `x_j(t+1) = x_j(t) XOR G(x_{j-1}, …, x_{j-d})`, `x_j = 0` for `j < 0`,
arbitrary Boolean `G`. Lemma 0's argument uses nothing about rule 30, so:

> **Theorem 5.1.** In any triangular XOR-driven system every coordinate is
> periodic in `t` — exactly, from `t = 0` — with period a power of two.

*Proof.* Verbatim Lemma 0: coordinate `j` is the running XOR-sum of its driver,
which is a Boolean function of coordinates below `j`, themselves periodic with
power-of-two periods by induction; the lcm of powers of two is the largest, and
a running XOR-sum of a period-`Q` signal has period `Q` or `2Q`. ∎

Break triangularity — let a driver read a coordinate at or above its own — and
this fails: periodicity becomes eventual at best.

### 5.1 The exact spectrum of a linear driver

Take `x' = x XOR Σ_{s∈S} (x << s)`, and write `P(z) = Σ_{s∈S} z^s` over GF(2),
so one generation is multiplication by `1 + P` and the single-cell orbit is
`x(t) = (1+P)^t`.

> **Theorem 5.2.** `p(j)` is the least `k` such that `[z^j](F_T · Q_k) = 0` for
> every finite subset `T`, where `Q_k := P(z^{2^k})` and
> `F_T := Π_{i∈T} (1 + Q_i)`.

*Proof.* By Frobenius over GF(2), `(1+P)^{2^k} = 1 + P(z^{2^k}) = 1 + Q_k`,
hence `x_j(t + 2^k) + x_j(t) = [z^j]((1+P)^t · Q_k)`. By Lucas, `(1+P)^t`
factors over the binary digits of `t` into a product of **distinct**
`(1 + Q_i)`, so the orbit is indexed by finite subsets `T`, and only indices
with `min(S)·2^i ≤ j` contribute to coefficient `j`. Coordinate `j` has period
dividing `2^k` exactly when all those differences vanish. ∎

This is algebraic: no generations are evolved. It is an exactness result, not a
speed result — affine drivers are tame, so their periods were never expensive
to simulate.

---

## 6. Verified but NOT proved

Nothing in this section is a theorem. Each is exhaustively checked over the
stated range and stated here so it is not mistaken for one.

**6.1 The tame set.** Over the 128 constant-free 3-input drivers
`G = L + N`, a driver appears to be tame — logarithmic period growth under
every seed — exactly when it is linear, or the indicator of a subcube that does
not contain the all-ones point. *128/128 agreement, 8 seeds, 40 coordinates,
periods to 2^13.* Route to a proof: Lemma 1 puts the low coordinates all-ones
at `t ≡ −1 mod 2^v`, so a subcube containing `(1,1,1)` fires exactly when the
run is solid and extends it, while one missing it never fires there.

**6.2 The closed form's range.** `p(j) = ceil(log2((j+1)/d))` appears exact
exactly when `min(S) = gcd(S)`. *Over all 127 shift sets `S ⊆ {1..7}`: fails on
53, precisely those with `min > gcd`.* Theorem 5.2 is exact on all of them —
6096 coordinates, 0 mismatches.

**6.3 The fresh rule.** Writing `σ_b` for the period-weight parity and `J_b`
for the joint AND-parity, `δ_b = σ_{b-1} ⊕ [p(b-2)=p(b-1)]·σ_{b-2} ⊕ J_b`
(24 coordinates, 0 mismatches). By Corollary 2.1, `σ_{b-1} = 0` when `b-1`
doubles, so the identity collapses to `δ_b = J_b` there — *8 coordinates,
0 failures.* The identity's derivation is routine; what is unproved is that the
measured `σ` and `J` values are the only content.

**6.4 The rendezvous, in general.** Lemma 1 is the rule 30 case. Across the 128
drivers the same property holds on the live coordinates of **76**. No predictor
found — `G(1,1,1) = 1` agrees on only 70 of 104. Open.

**6.5 The plateau sequence.** `#{b : p(b) = v}` — the primitive object; A094604
is its partial sums. Measured on 136 terms derived from A094604: longest gap
**8**; density **0.4191**, wandering (0.375 / 0.4219 / 0.4271 / 0.4191 at
32 / 64 / 96 / 136); **2-kernel classes 1, 2, 4, 8 through depth 3 with zero
collapse**, so no evidence of 2-automaticity; not periodic, not eventually
periodic, not a function of `b mod 2^k` for `k ≤ 7`, and uncorrelated with
Thue-Morse, Rudin-Shapiro and `v₂` parity. A rule for this sequence would make
every rung past 57 free. None found.

**6.6 The parity hierarchy.** Over a full period the shift `t → t+1` permutes
the window, so every monomial in the columns yields a linear relation among
monomial-parities. On a band closed downward no relation escapes, but the
system does not determine the higher parities from the singletons: at nine
columns, 511 monomials, GF(2) rank 462, **49 free**. Shift-invariance alone
collapses `2^N` unknowns to about `N²/2` — exponential to polynomial — and no
further. Lemma 1 (known values) and Lemma 2 (antipodal symmetry) supply
relation families that were **not** in that computation; re-running the rank
with them included is the untried experiment.

---

## 7. Verification record

The proofs above are unconditional. Each step was nonetheless checked before it
was written, over 30 bits and 65,536 generations.

| step | checked | failures |
|---|---|---|
| Lemma 0 — `p` nondecreasing, steps in {0,1} | 30 bits | 0 |
| (D) — driver weight odd on doubling bits, even otherwise | 29 bits | 0 |
| (2) — the bracket vanishes | 145 cases | 0 |
| (3) — `v₂(Q) = p(b-1)` and `c_{b-1}(Q−1) = 1` | 29 cases | 0 |
| Lemma 1 — the phase lemma | 181,072 cases | 0 |
| Lemma 2 — antipodal | 12 doubling bits | 0 |
| Corollary 2.1 — weight exactly `T/2` | 10 bits | 0 |
| Lemma 3 — the run terminates | 65,520 cases | 0 |
| Theorem 5.2 — exact spectrum | 6,096 coordinates | 0 |
| Corollary 4.1 — cross-derivation vs published terms | 70 terms | 0 |

Reproduce: `sitoolz.edgeproof`, `sitoolz.driver_lab`, `sitoolz.parity_lab`,
`sitoolz.oeis_crosscheck`.

---

## 8. Prior work

**Wolfram**, *A New Kind of Science* (2002), p. 871 — observes the period
doubling of the diagonals on the left side of rule 30.
**Rowland**, *Local nested structure in rule 30*, Complex Systems **16** (2006)
239–258 — proves Wolfram's left-side doubling criterion (Proposition 2), and
in Lemma 2 gives immediate periodicity with a divisibility bound on column
periods for right-bijective rules. For 2 colours, right-bijective and
triangular XOR-driven are the same class.
**Klimov & Shamir**, CHES 2002 — T-functions; the same class under another
name, studied for provable period rather than for the spectrum.
**OEIS** A094603–A094606, since 2004; the v₂ statement contributed as a
conjecture by Alan Michael Gómez Calderón, 2023-02-17.

Proposition 2 is a theorem about the **left** side of rule 30; everything
above concerns the **right** light-cone edge, where Lemma 2 supplies the
periodicity and the doubling criterion (D), the phase lemmas and the v₂ law
are proved here from the recurrence — self-contained, citing nothing. The
side-by-side comparison is in [EDGE_LAW.md](EDGE_LAW.md).
