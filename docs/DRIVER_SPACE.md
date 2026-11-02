# Driver space — where the complexity of rule 30's edge lives

*Map the Boolean drivers of triangular XOR-driven systems: which give a closed-form period spectrum, which grow like rule 30, and what the doubling bit is made of.*

**Lars Rocha**

---

## The setting

The edge-law theorem ([AN_INVITATION.md](AN_INVITATION.md),
[PROOFS.md](PROOFS.md) §5) fixes every period of a triangular XOR-driven
system

```
x_j(t+1) = x_j(t) XOR G(x_{j-1}, …, x_{j-d}),     x_j = 0 for j < 0
```

as a power of two, for any Boolean driver `G`. It says nothing about
**which coordinates double**. That profile — the spectrum `p(j)` — is the
open object, and driver space is the laboratory for it: the theorem holds
for every `G`, so each driver is one experiment.

Rule 30's light-cone edge is the driver `G = a OR b` with `a = x_{j-1}`,
`b = x_{j-2}`. Over GF(2),

```
OR(a,b)  =  a XOR b XOR ab   =   XOR  ⊕  AND
```

a linear part plus a quadratic part. Rule 90 is the linear part alone.

---

## What was measured

| experiment | scale | code |
|---|---|---|
| all 16 two-input drivers, spectrum and clean/opaque split | 40 coordinates, periods to `2^11` | `tri_xor.survey`, `tri_frame.survey_frame` |
| all 128 constant-free three-input drivers, tame or not | 8 seeds each (4 dense, 4 sparse), 40 coordinates, periods to `2^13`, classified on the envelope `max p(0..j)` | `driver_lab.mixture_deep` |
| every shift set `S ⊆ {1..7}`, linear driver, exact spectrum vs measurement | 127 sets, 6096 coordinates | `driver_lab.check_exact`, `txs.spectrum` |
| the three parity terms of the doubling bit, antipodal lemma, fresh rule | 24 coordinates, periods to `2^14` | `parity_lab.report` |
| the rendezvous (in-phase) property across drivers | 128 drivers | `txs.verify_rendezvous`, `txs.in_phase_profile` |
| sparse "shoulder" states per driver and depth | 8,192 generations, 32 coordinates | `tri_frame.shoulder_frame`, `tri_frame.event_frame`, `tri_frame.anf_frame` |
| the plateau sequence of A094604, kernel and correlation tests | 58 published terms, 136 doubling terms | `oeis_crosscheck` |
| rungs 31–43 of the edge table, by direct sweep | `2^43` generations | `edgelaw_deep` |

The bit-parallel evaluator in `driver_lab` is cross-validated against the
list-indexed one in `tri_xor` before use — 0 disagreements across the 2-input
drivers (`driver_lab.crossvalidate`).

---

## Two negative results

Both close a question rather than leave it open, and each kills a class of
hypotheses instead of one.

### The spectrum is a period ladder, not a mixture

Rule 30's driver splits as `XOR ⊕ AND`. Measured separately, at 64
coordinates:

| driver | envelope of p at j = 8, 16, 32, 63 | verdict |
|---|---|---|
| `XOR = a+b` | 4, 5, 6, 6 | logarithmic — +1 per doubling |
| `AND = ab` | 0, 0, 0, 0 | inert from a single-cell seed |
| `OR = a+b+ab` | 5, 8, 12, — | linear — outruns `2^13` by j=36 |

The linear part grows logarithmically on its own. The quadratic part does
nothing at all on its own from this seed. **Their sum grows linearly.** The
spectrum of `OR` is not assembled from the spectra of its parts; it is a
ladder climbed one parity decision at a time — bit `b` doubles exactly when
its driver has odd weight over one period (C2), and that parity depends on
the whole interaction. Nor is "tame linear part + any nonlinear part ⇒ linear
growth" a rule: over the 128 three-input drivers the exceptions are
characterised, and there are exactly 12 of them (next section).

### The spectrum is not 2-automatic

Plateau length at rung `v` = how many bits share `p(b) = v` = the first
difference of A094604. Over its 58 published terms:

```
v= 0..14 : 1  2  1  2  1  2  6  1  8  1  2  2  5  2  1
v=15..29 : 2  2  2  5  1  2  3  1  3  2  3  1  2  3  1
v=30..44 : 2  2  3  2  1  2  2  2  4  1  2  7  3  1  2
v=45..57 : 2  3  1  2  5  1  2  2  3  1  5  2  2

histogram:  1×16   2×26   3×8   4×1   5×4   6×1   7×1   8×1
sum = 137 = A094604[57]                    mean = 2.362
```

Forty-two of fifty-eight are 1 or 2 — an orderly `1,2,1,2,1,2` rhythm through
`v=5`. Then it breaks. At v=6 the plateau is 6; at v=8 it is 8 — the only two
places where plateau equals `v`, and then v=12→5, v=18→5, v=41→7, v=49→5,
v=55→5 kill any such rule. The minimum plateau is 1, so every rung is
attained — no period `2^v` is skipped, which is Corollary 4.2 of
[PROOFS.md](PROOFS.md).

| test | result |
|---|---|
| correlation with `v₂(v)` | means 2.28, 2.29, 2.43, 3.50, 2.00, 3.00 — no signal |
| autocorrelation, lags 1–14 | best 41% at lag 7; nothing periodic |
| plateau `= v`? | only at v=6 and v=8; coincidence |
| **2-kernel of the doubling sequence** | **does not collapse** |

The last one is the real result. Take the raw object rather than the
compressed one: `D(b) = 1` exactly when bit `b` doubles — which by C2 *is*
the parity of the driver weight. 136 terms, 57 ones (41.9%):

```
1011011010000011000000011010100001011010101000011010011001010011010011010100
101101010100011010000001001101010011010000110101001100001010
```

Its 2-kernel — the set of subsequences `D(2^k n + r)` — grows without
collapsing:

```
k ≤ 1:   3 distinct
k ≤ 2:   7 distinct
k ≤ 3:  15 distinct
```

Exactly `2^(k+1) − 1`. Every subsequence is new; nothing merges. A
2-automatic sequence has a finite 2-kernel. This one, on the evidence
available, does not. Further checks on the same object: not periodic, not
eventually periodic, not a function of `b mod 2^k` for `k ≤ 7`, and
uncorrelated with Thue–Morse, Rudin–Shapiro and `v₂` parity; density 0.4191,
wandering (0.375 / 0.4219 / 0.4271 / 0.4191 at 32 / 64 / 96 / 136).

> **There is no closed form for the spectrum as a finite automaton on the
> binary digits of `b`.** Any formula for `p(b)` must go through the driver
> parity.

So C2 is not a stepping stone toward something simpler — it is the answer.
The chain measurement → spectrum → closed form stops at "spectrum", and it
stops for a reason. 136 terms cannot *prove* non-automaticity — a finite
kernel could still close beyond the data — so this is strong evidence and a
falsifiable claim, recorded as measured, not derived.

---

## The driver families

### Two inputs: eleven clean, five opaque

Running the recurrence with each of the sixteen 2-input Boolean drivers, 40
coordinates, periods resolved to `2^11`:

| driver | coordinates resolved | spectrum |
|---|---|---|
| `XOR`, `XNOR`, `a`, `b`, `NOT a`, `NOT b`, `a AND NOT b`, `b AND NOT a`, `AND`, `0`, `1` | 40/40 | explicit |
| `OR`, `NOR`, `NAND`, `(NOT a) OR b`, `a OR (NOT b)` | 20–29/40 | outruns `2^11` |

In truth-table codes the clean set is {0,2,3,4,5,6,8,9,10,12,15} and the
opaque set {1,7,11,13,14} = {NOR, NAND, ~a|b, a|~b, OR}. Rule 30's edge is
OR (code 14), opaque; rule 90 is XOR (code 6), clean, with
`p(j) = ⌈log2(j+1)⌉`.

Under a single-cell seed the two-input split is characterised, 16/16, by

```
G gives a clean spectrum  ⟺  G is affine
                              OR  G is a subcube indicator whose subcube
                                  excludes the all-zeros point (0,0).
```

The two ends are one bit apart in places. `AND` freezes — nothing ever
propagates. `NAND`, which is `AND` XOR the constant 1, is opaque. Adding a
constant to the driver moves it from explicit to intractable: NOR is the
indicator of `(0,0)`, so `NOR(0,0) = 1` flips every quiet cell and the
system goes OR-like.

### Three inputs: the tame set is exactly two families

Take the 3-input drivers with no constant term, written over GF(2) in
`a = x_{j-1}`, `b = x_{j-2}`, `c = x_{j-3}`:

```
G  =  (linear part L ⊆ {a, b, c})  +  (nonlinear part N ⊆ {ab, ac, bc, abc})
```

128 drivers. Rule 30's light-cone edge is `L = {a,b}, N = {ab}`. Each was
run from four dense seeds and four sparse ones, 40 coordinates, periods
resolved to `2^13`, classified on the **envelope** `max p(0..j)` rather than
raw `p(j)`.

**20 drivers are tame under every seed. 108 reach linear growth under at
least one.** The 20 split cleanly:

> **The tame set.** A driver is tame exactly when it is **linear**, or it is
> the indicator of a **subcube that does not contain the all-ones point**.
>
> Verified: **128/128 drivers agree**, 8 seeds each.

Both clauses are needed and neither contains the other:

- `a`, `b`, `c` are subcube indicators that *do* contain `(1,1,1)`, and are
  tame anyway — they are linear, and the Frobenius argument does not care
  about the seed.
- `ab¬c`, `a¬bc`, `¬abc`, `a¬b`, `a¬c`, `a¬b¬c`, `b¬a`, `b¬c`, `b¬a¬c`,
  `c¬a`, `c¬b`, `c¬a¬b` are not linear, and are tame. Those 12 are the
  exceptions to "tame linear part + any nonlinear part ⇒ linear growth".

The near misses say the boundary is sharp rather than fuzzy:

| driver | subcube? | contains (1,1,1)? | verdict |
|---|---|---|---|
| `ab` | yes | yes | **linear growth** |
| `ab¬c` | yes | no | tame |
| `abc` | yes | yes | **linear growth** |
| `ab + ac` = `a(b+c)` | no | no | **linear growth** |

That last row matters: `G(1,1,1) = 0` on its own is *not* enough. The driver
has to be a single conjunction of literals, and that conjunction has to miss
the all-ones point.

**Why all-ones and not all-zeros.** For constant-free drivers
(`G(0,0,0) = 0`) the vacuum point is automatically excluded, so the
two-input "excludes all-zeros" framing and the "avoids all-ones" framing
coincide; they diverge exactly on NOR (which has a constant term) and on the
degenerate seed-death case AND. Under a single seed, AND-type drivers freeze
and look clean; under every seed they reach linear growth and are correctly
non-tame. The seed-robust characterisation is the all-ones one. A looser
reproduction — single seed, "clean = all bits resolve to ≤ 2^11" — counts
70 of 128 as clean; the envelope test under eight seeds gives 20, and the
difference is entirely those two relaxations.

**The mechanism, measured.** The exceptional drivers do not slow the growth
down — they *sparsify*. `x' = x XOR (b AND NOT c)` freezes every odd
coordinate and, on the even ones, reproduces the spectrum of the pure linear
driver `x' = x XOR (x<<2)` coordinate for coordinate. `c AND NOT a` freezes
everything off the multiples of 3 and runs the shift-3 linear spectrum on
what is left. The nonlinear term acts as a filter on which coordinates
participate, not as an accelerant.

**Why this is the interesting side of the question.** The all-ones point is
exactly where the phase lemma C3 puts the low coordinates: at
`g ≡ -1 mod 2^v` the run below the cut is solid. A driver whose subcube
contains `(1,1,1)` fires precisely at the moment the run is solid, and
extends it. One that misses it never fires there. That is a mechanism, and
it is the thing to turn into a proof.

### The rendezvous, across drivers

The phase lemma (C3: in phase means 1) is proved for rule 30's OR. Across
the 128 drivers the same property holds on the live coordinates of **76**.
`G(1,1,1) = 1` predicts it on only 70 of 104; no predictor found.
`txs.verify_rendezvous` and `txs.in_phase_profile` run the test per driver.

---

## Linear drivers: the spectrum, exactly, for every shift set

Take `x' = x XOR Σ_{s∈S} (x << s)`, a shift polynomial over GF(2). The map
is `1 + P` with `P(z) = Σ_{s∈S} z^s`, and by Frobenius

```
(1 + P)^(2^k)  =  1 + P(z^(2^k))
```

so if `P`'s lowest shift depth is `d`, the `2^k`-th iterate only reaches
coordinate `j` while `d·2^k <= j`. That predicts `p(j) = ceil(log2((j+1)/d))`.

**For `d = 1` this is exact.** Verified for shift sets
`{1}, {1,2}, {1,3}, {1,2,3}`, under both a single-cell seed and a dense
random seed, to 95 coordinates — zero mismatches in all eight runs.

Past `d = 1` the system splits into interleaved subsystems on residue
classes and `p(j)` goes non-monotone (`0,0,0,1,0,0,2,0,0,2,…` for
`x' = x XOR (x<<3)`). The failure is not a residue-class artifact needing a
per-class patch:

> **The closed form is exact exactly when `min(S) = gcd(S)`.**
>
> Verified over all 127 shift sets `S ⊆ {1..7}`: the form fails on 53 sets,
> and those 53 are precisely the sets with `min(S) > gcd(S)`.

When `min(S) > gcd(S)` the form is only an upper bound, because coordinates
below the lowest shift can cancel against each other.

The exact statement replaces it. The orbit from a single cell is
`x(t) = (1+P)^t`. Frobenius gives `(1+P)^(2^k) = 1 + P(z^(2^k)) =: 1 + Q_k`,
hence

```
x_j(t + 2^k)  +  x_j(t)   =   [z^j] ( (1+P)^t · Q_k )
```

and by Lucas, `(1+P)^t` factors over the binary digits of `t` into a product
of **distinct** `(1 + Q_i)`. So the whole orbit is indexed by finite subsets
`T` of the indices with `min(S)·2^i ≤ j` — about `log2(j)` of them — and

> **Exact rule.** `p(j)` is the least `k` such that `[z^j](F_T · Q_k) = 0`
> for every subset `T`, where `F_T = Π_{i∈T} (1 + Q_i)`.
>
> Verified: **127 shift sets, 6096 coordinates, 0 mismatches** against
> measurement. `driver_lab.spectrum_linear_exact`, `txs.spectrum`.

This is algebraic, not evolved: no generations are run at all. It is an
exactness win rather than a speed win — affine drivers are tame, so their
periods are small and simulating them was never expensive — and it is the
structural reason the old form worked when it worked. `txs.synthesize` runs
it backwards: name a spectrum, get the shift sets that produce it.

So: **affine drivers are tame — logarithmic period growth, closed form at
depth 1, exact rule at every depth.** This is provable, not merely measured,
and seed-independent, because the Frobenius argument does not mention the
initial condition.

---

## What the doubling bit is made of

Coordinate `j` is the running XOR-sum of its driver `d_j = c_{j-1} OR c_{j-2}`,
so the period doubles exactly when `d_j` has odd weight over one period. With
`w(A OR B) = w(A) + w(B) − w(A AND B)` that is three parity terms. Write
`σ_j` for the weight parity of coordinate `j` over one of **its own** periods
and `J_j` for the joint term:

```
δ_j  =  σ_{j-1}  ⊕  [p_{j-2} = p_{j-1}]·σ_{j-2}  ⊕  J_j            (D)
```

Verified: **24 coordinates, 0 mismatches.** Two of those three terms are then
free:

> **Antipodal lemma.** If coordinate `j` doubles, then
> `c_j(t + T_j/2) = NOT c_j(t)` for **every** `t`.
>
> *Proof.* `c_j(t + T_j/2) ⊕ c_j(t)` is the XOR-sum of `d_j` over a window of
> length `T_j/2`, which is exactly one period of `d_j`. That sum is 1 for
> every `t` — odd weight over one period is what made `j` double. ∎
>
> **Corollary.** A doubling coordinate has weight exactly `T_j/2`, so
> `σ_j = 0` whenever `p_j ≥ 2`.

Verified: 10 doubling coordinates, 0 failures, both statements.

> **Fresh rule.** If coordinate `j-1` doubled and `p_{j-1} ≥ 2`, then
> `σ_{j-1} = 0` and `p_{j-2} < p_{j-1}`, so (D) collapses to
>
> ```
> δ_j  =  J_j
> ```
>
> The doubling of a coordinate that follows a doubling is *entirely* the
> joint AND-parity. Nothing else survives.

Verified: 8 coordinates, 0 failures.

The wall was never in the OR. It is one term —
`J_j = parity of Σ_t c_{j-1}(t)·c_{j-2}(t)` — and on the coordinates that
follow a doubling it is not merely the hard part, it is the whole quantity.
Anything that computes `J` cheaply computes `δ` cheaply on a third of the
coordinates immediately. The remaining coordinates sit on plateaus and carry
a `σ` from a non-doubling coordinate as well, and those `σ` values are the
next thing to characterise. [CENTER_LAW.md](CENTER_LAW.md) carries this
further: the Sierpiński-coefficient form of `δ`, the post-doubling law, and
the disjoint-convolution form of `J`.

The obvious attack — expand `Σ_t c_j d_{j-1}` by shift-invariance over a full
period — produces triples at the first step, so closure is not free. On a
band of nine columns, 511 monomials, the shift-invariance relations have
GF(2) rank 462, leaving **49 free** ([PROOFS.md](PROOFS.md) §6.6).

---

## Sparse shoulders: a constant term interacting with depth

Seven of the sixteen two-input drivers produce sparse "shoulder" states
(`2^a + 2^b`); nine never do. No predicate on the algebraic normal form
separates them — every conjunction to three literals and every linear
combination was searched. One row pair looked like the explanation:

```
NOT a  =  1 + a   ->  128 sparse states
NOT b  =  1 + b   ->    0 sparse states
```

Same function, same degree, same weight, opposite behaviour — which reads as
depth-coupling, the driver attaching to `a = x_(j-1)` at depth 1 versus
`b = x_(j-2)` at depth 2. The control shows it is narrower than that. Run
the *bare* shift at each depth alongside the *negated* one:

```
BARE      G = x_(j-d)              solid    2-bit   sparse?
  d = 1                             1536     1280    yes
  d = 2                              512     2048    yes
  d = 3                              512     2048    yes

NEGATED   G = 1 + x_(j-d)          solid    2-bit   sparse?
  d = 1                              128      128    yes
  d = 2                              256        0    NO
  d = 3                              512        0    NO
```

**Bare shifts produce sparse states at every depth.** Depth alone does nothing.
The negated ones survive only at depth 1.

> It is not depth. It is the **constant term interacting with depth**, and it
> survives only at depth 1.

There is a mechanism visible in that. `G = 1 + x_(j-d)` means *fire when
`x_(j-d)` is zero*. At depth 1 the cell being read is the immediate neighbour,
so "fire on zero" tracks the leading edge and can still leave a sparse state
behind it. At depth ≥ 2 the driver fires on a cell two or more back — cells the
edge has already passed — so it keeps flipping settled territory and nothing
stays sparse.

That also accounts for the outlier the whole question started from:

```
a  (rule 90):   1536 solid   1280 two-bit   1280 shoulder pairs
OR (rule 30):      6 solid      0 two-bit      0 pairs
```

`OR = a + b + ab` fires whenever *either* cell below is set, so past the leading
edge it is firing almost always. Six solid runs in 8,192 generations and never a
sparse state. Its three one-bit neighbours — `XOR`, `b`, `a` — each drop one of
those terms and immediately regain sparsity.

The mechanism is a reading, not a proof. What would settle it is writing down
explicitly what a two-bit state requires under `x' = x XOR G(...)` — all but
two coordinates zero simultaneously — and showing the condition fails for
negated shifts at depth ≥ 2. That step is small.

---

## Rungs 31–43, from the recurrence alone

A compiled sweeper written from the recurrence — 256-bit window, four
`uint64` words, `E' = E ^ ((E << 1) | (E << 2))`, one sequential pass reading
the run at every `g = 2^v − 1` it passes, ~670 million generations per second
— reproduced rungs 0–30 exactly as `edgelaw.TABLE` has them, which is what
licenses the rest:

| rung | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 39 | 40 | 41 | 42 | 43 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| edge run | 74 | 77 | 79 | 80 | 82 | 84 | 86 | 90 | 91 | 93 | 100 | 103 | 104 |
| plateau | 2 | 3 | 2 | 1 | 2 | 2 | 2 | 4 | 1 | 2 | **7** | 3 | 1 |

Checked against OEIS A094604 afterwards, all thirteen match exactly, and
A094604 runs on to rung 57. `edgelaw_deep` carries these; `edgelaw` keeps
the original table. The law covers every generation to
**8,796,093,022,207**. Rung 41's plateau of 7 is the largest since rung 8;
rungs 42 and 43 came back to 3 and 1, so the 7 is a spike rather than a
change of regime. Rung 43 took 4h 12m; each rung doubles.

The value here is the independent agreement: a 2.2-trillion-generation sweep
written from the recurrence, with no reference to the published sequence,
landing on all thirteen of its terms.

---

## Measurement rules

Three things that make the measurements above trustworthy, each learned the
expensive way:

- **Shift direction.** In the bit-parallel evaluator, `x << s` reads the
  coordinate *below* `j`; `x >> s` reads above it and runs the non-triangular
  system. Every evaluator here is cross-validated against the list-indexed
  implementation in `tri_xor` before use; all sixteen 2-input spectra agree
  exactly.
- **Classify on the envelope.** `p(j)` is not monotone for most drivers.
  Comparing raw endpoints calls `x' = x XOR (x<<3)` wild when it is cleanly
  logarithmic. Use `max p(0..j)`, never raw `p`.
- **Seed from dense rows as well as a single cell.** A driver with only
  degree-≥2 monomials fires only when two of the lower coordinates are set at
  once, and a single seeded cell never produces that, so the system freezes
  and looks tame. Re-run from a dense seed, 11 of 15 purely nonlinear
  drivers go wild. The affine ⇒ tame direction is the only one
  that is seed-independent.

---

## Still open

1. What property of a driver produces linear rather than logarithmic period
   growth? Nonlinearity is necessary (affine ⇒ tame, proved); the tame set
   says exactly which nonlinear drivers escape (measured, 128/128); a proof
   via the all-ones mechanism is the next step.
2. The rendezvous property across drivers: 76 of 128 hold it, and no
   predictor is known.
3. For rule 30's `OR`: are the gaps in the doubling indicator bounded, does
   its density converge from 0.4191, and does the 2-kernel ever close? The
   evidence to 136 terms says no finite automaton; a proof would need the
   driver parity.
4. The joint parity `J_j`: compute it without summing the window.
5. The sparse-shoulder mechanism, as a proof rather than a reading.

---

## Reproducing

```python
from sitoolz import tri_xor as tx
tx.spectrum(tx.XOR)    # logarithmic
tx.spectrum(tx.AND)    # inert from a single cell
tx.spectrum(tx.OR)     # rule 30's edge — linear
tx.show_survey()       # all sixteen 2-input drivers

from sitoolz import driver_lab as dl
dl.crossvalidate()                       # bit-parallel vs list-indexed, 0 disagreements
dl.mixture_deep()                        # the tame set, 128 drivers x 8 seeds
dl.check_exact(smax=7)                   # exact rule vs measurement
dl.spectrum_linear_exact((3,5,8), 5000)  # a spectrum with nothing evolved

from sitoolz import txs
txs.spectrum((1,3))                      # exact spectrum of a shift set
txs.synthesize([0,1,1,2,3,3])            # which shift sets produce a given spectrum
txs.verify_rendezvous(3, 1)              # the in-phase property for one driver

from sitoolz import parity_lab
parity_lab.report()                      # identity (D), antipodal lemma, fresh rule

from sitoolz import tri_frame as tf
tf.shoulder_frame()                      # solid runs and two-bit states, all 16
tf.event_frame(tx.A)                     # when the shoulder events happen
tf.anf_frame()                           # algebraic normal form vs shoulder capability

from sitoolz import edgelaw_deep         # the table to rung 43
from sitoolz import oeis_crosscheck      # A094604's published terms and the cross-derivations
oeis_crosscheck.check()
```
