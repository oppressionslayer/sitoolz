# The Center Law — Rule 30's center column, located and metered

*Locate rule 30's center column inside a bank of exact power-of-two clocks, index any near-edge cell at generation 10^100 in milliseconds, and price the center problem as one exponential staircase.*

**Lars Rocha** · companion to [EDGE_LAW.md](EDGE_LAW.md). Everything below is
verified by code in `sitoolz.centerrandom` unless marked open.

## The law

    center(t)  =  bit t of E(t),      E' = E ^ ((E << 1) | (E << 2)),   E(0) = 1

The edge recurrence, run at full width, **is Rule 30 in the frame co-moving with
the right light-cone edge**: bit `j` of `E(t)` is the cell at distance `j` behind
the edge. The center column is the **anti-diagonal** through this system —
sample coordinate `t` at time `t`. Verified for 2,000 generations against the
direct simulation (which is itself verified against OEIS A051023).

Consequences, all verified:

1. **Every diagonal of Rule 30 parallel to the right edge is exactly periodic**
   from generation 0, with period `2^p(j)` (the `tri_xor` theorem applied to the
   full cone; 25 coordinates checked, 0 failures). The "chaotic" triangle is a
   bank of exact power-of-two clocks. The famous randomness is a property of
   the *cut*, not the material: cut horizontally → pure clocks; cut on the
   anti-diagonal → the center column.

2. **Why the center looks random — quantified.** The center reads clock `t`
   exactly once, at time `t`, but clock `t` has period `2^p(t)` with
   `p(t) ~ t/2.4`. The sample lands inside the clock's *first* period —
   under-sampled by ×4 at t=7, ×20,971 at t=50, ×43,980,465,111 at t=100.
   Not chaos: **under-sampling of exact clocks.**

3. **Three center bits are pinned by the phase lemmas** (C3, in-phase → 1):
   the sample hits the in-phase slot exactly when `2^p(t) | t+1`, which happens
   at t = 0, 1, 3 (and then never again below the horizon). All three verified.

## The wedge index (random access into Rule 30)

Triangularity makes the low-`k` window of E **self-contained**, hence periodic
with period `2^p(k)`. So time reduces like an index:

    E(g) mod 2^k  =  E(g mod 2^p(k)) mod 2^k        -> reduce, walk the remainder

`cell(g, j)` / `strip(g, k)` give **any near-edge cell at ANY generation** —
10^100 as easily as 100 (walks are masked to k bits; milliseconds). Verified by
700 random probes against direct simulation, all matching.

Fun consequence: any `g` divisible by `2^20` has its 50-cell wedge **identical
to generation 0**. Generations 10^20 through 10^100 all have déjà vu.

## The left wedge

The left light cone is not triangular (coordinates read themselves), so left
columns are only **eventually periodic** (Rowland 2006) — but that still gives
a lookup after the transient. Measured to offset 199: periods are tiny (≤ 8),
transients grow ~1.3·j. `left_cell` / `left_strip` verified 600/600 against an
independent simulation; thin-evidence entries refuse to answer rather than
guess (offset 65 is provisional).

Cost profiles are mirror images:

|                  | RIGHT wedge          | LEFT wedge          |
|------------------|----------------------|---------------------|
| periods          | 2^48 (exponential)   | ≤ 8 (tiny)          |
| transients       | ~0                   | ~1.3·j (linear)     |
| width            | 113                  | 199                 |
| cost to widen    | 2^v per rung         | linear (simulate)   |

## The metering (the main result)

The center bit at generation `t` sits at offset `t` behind the right edge, so

    center(t) = cell(t, t)      — valid exactly for t <= 113.

`cell(113,113)` works; `cell(114,114)` raises KeyError. Therefore:

    center bits known  =  right-wedge depth  =  rungs of A094604

* 114 center bits are indexed (60 ones, 54 zeros).
* Each new rung `v` of A094604 buys center bits up to `TABLE[v] − 1`.
* Each rung costs `2^v` **sequential** steps (rung 49 ≈ a week of optimized C;
  rung 60 ≈ centuries). The index reduction *consumes* periods; it cannot
  discover them.
* The center escapes the LEFT wedge at the same rate (it is offset t from both
  edges) — it swims up the dark middle, equidistant from both lit shores.

**Wolfram's center problem, restated:** the wedge-widening problem. One
exponential staircase, clearly priced, instead of a mystery.

## Spec: what a new structure must do to beat the wall

Any structure that cracks (or dents) the center must do one of these — this is
the precise target list:

1. **Discover periods without elapsing them.** Produce `p(j)` for `j ≥ 114`
   in under `2^p(j)` sequential steps. Equivalent to extending A094604 cheaply.
   (The doubling indicator δ(j) is the parity of the OR-driver's weight over one
   period — any algebraic handle on that parity beats the staircase.)

2. **Decide Q2: is δ 2-automatic?** If yes, the whole period ladder becomes a
   finite automaton and rungs are free. 114 terms cannot decide it; a proof
   could. This is the highest-leverage single question.

3. **Jump the anti-diagonal.** Evaluate bit `t` of `E(t)` in `o(t)` work
   without knowing p(t) — a shortcut past the walk itself. (Same shape as the
   VDF / repeated-squaring wall from the Lucas–Lehmer work: transparent chain,
   distance is the only guard.)

4. **Find a third frame.** The right frame is triangular (exact clocks); the
   left frame is self-reading (eventual clocks). A frame in which the CENTER's
   coordinate reads only below itself would make the center a T-function and
   the theorem fires. Nobody has exhibited one; nobody has proven one cannot
   exist. This is the door a genuinely new object could open.

## What the hunt established

Results of the computational hunt against the spec targets. Positives first;
closed routes after, so nobody re-walks them.

### The parity handle (verified 48/48)

δ(j) = parity of the OR-driver's weight over one window 2^p(j-1) — the
[AN_INVITATION.md](AN_INVITATION.md) formulation, confirmed for every
resolvable coordinate j ≤ 50. The whole staircase reduces to predicting one
weight-parity per rung.

### par_and: the refined target

Decomposing δ(j) = par_u(j-1) ⊕ par_v(j-2) ⊕ par_and over the window shows
par_u is sparse (11/48 ones) — the content concentrates in **the parity of
coincidences of two adjacent clocks over one window**. That correlation
parity is the single quantity to attack: predict it without summing the
window, and the rungs are free. (It is itself unstructured at accessible
scales: not constant, no small period.)

### The Sierpiński algebra of the staircase (discovered and verified)

In subset-zeta (XOR-Pascal / binomial mod 2) coordinates — Lars Rocha's
Sierpiński transform, c_i = XOR of S(t) over t ⊆ i — the ladder becomes exact
algebra (each law verified on random words AND on the real system):

  - prefix-XOR (the ladder integrator)  =  index shift up by 1
    (hockey-stick identity mod 2: Σ_{s<t} C(s,i) = C(t,i+1))
  - pointwise AND (the OR-driver's nonlinearity)  =  OR-convolution of
    coefficients (Lucas: C(t,i)·C(t,j) ≡ C(t, i|j) mod 2)
  - weight parity over window 2^b  =  the TOP coefficient c_{2^b−1}
  - period 2^k  ⟺  coefficients supported below index 2^k

Hence **δ(j) = the top Sierpiński coefficient of the driver word — verified
48/48 on all resolvable rungs.** The staircase is diagonalized by the same
XOR-Pascal object as the square-index work. The hard core, restated exactly:
*predict the top coefficient of an OR-convolution without materializing it.*

### The post-doubling law (theorem)

Derived from the antipodal lemma of
[THE_SPECTRUM_BREAKTHROUGH.md](THE_SPECTRUM_BREAKTHROUGH.md), verified
perfectly: if bit `b` doubles, then

  (i)  its period-word has even weight: `w_b = 0`   (all doubled bits pass), and
  (ii) **the next rung is free:  δ(b+1) = w_{b−1}** — the weight-parity of the
       sub-bit's word, a half-size quantity.  **Verified 19/19.**

Proof shape: antipodality `c_b(t+Q/2) = ¬c_b(t)` telescopes the OR-weight —
the pure terms vanish (each half contributes Q/2 ≡ 0), and the product term
collapses to Σ c_{b−1} over its own full period. No window sum. This predicts
every post-doubling rung (~42% of all rungs) from strictly smaller data — the
first brute-force-free rung predictions. Held-bit subcases probed for further
collapses (grouped by δ(b−1), half-window candidates): all negative — the
held-pair coincidence parity remains the irreducible core.

### The full-degree reformulation

The subset/Möbius transform IS the algebraic normal form over GF(2), so the
"top Sierpiński coefficient" statement becomes: **bit b doubles ⟺ its driver
c_{b-1} ∨ c_{b-2}, as a Boolean function of the time-bits, has FULL algebraic
degree** (its ANF contains the product of all p(b-1) input bits). Makes
density ≈ ½ expected; reframes doubling as a degree condition. Clock words
themselves run near-maximal ANF degree (deficit p−deg ∈ {0,1,2}).

### parAND as a disjoint convolution (verified 17/17)

The irreducible coincidence parity equals `⊕` over *disjoint* index pairs
(i,j) of C̃ᵢ^{b-1}·C̃ⱼ^{b-2} (complemented Sierpiński spectra) — equivalently
the covering-pairs parity i∪j = full. Cleanest closed form of the core, but
still costs the full window (spectra are dense).

### Independent confirmation of rungs 29–40

`edgereach.reach(40)` independently confirms rungs 29–40 against A094604
(MATCH, 737M gen/s).

### Routes closed

* **The word-kernel route.** The one-period words S_j are NOT uniformly
  automatic: their dyadic-decimation kernel grows with length (3 → 8218 for
  j = 2..40, tracking ~10–20% of word length). No finite machine lives in raw
  time-decimations of the orbits.

* **Small-DFA 2-automaticity of δ (leaning closed).** Myhill–Nerode probe on
  the 113 known δ terms (both MSB- and LSB-first reading of binary(j)) needs
  ≥ ~27 state-classes with thousands of hard conflicts. If δ is 2-automatic,
  its automaton is not small; the easy win is off the table. (Large machines
  are not excluded — 113 terms cannot do that.)

* **The sparsity route.** The words' Sierpiński coefficient vectors are DENSE
  (20–50% ones; ~393k ones at j=50), so iterating the exact algebra
  symbolically costs the same exponential. No free lunch from sparseness.

* **Direct rung-sequence prediction.** A094604 grows near-linearly
  (~2.4–2.6/rung) but the residual drifts (held-out errors reach ±10 by rung
  48); gaps are irregular; square-index lenses (fold family, XOR of
  consecutive rungs) show no membership law. δ needs exact rungs, so
  approximate trend-fits buy nothing.

* **2-regularity of the rungs.** Kernel-rank test on the gap sequence is
  data-starved (rank fills the row count at every reachable depth); no small
  linear representation visible in 48 terms. Inconclusive, leaning no — needs
  more rungs to decide.

* **δ(b) has no cheap predicate.** popcount, v2, mod-3, digit tests: best is
  popcount(b+1) at 65%, no clean law. A094604 is not stride-self-similar.

* **Diagonals do NOT nest into few germs.** A loose prefix-match suggested
  every clock decimates to one of 4 germs; the strict canonical test refutes
  it — distinct germs GROW with period (2,4,13,31,38,…), and the apparent
  collapse at large periods is data starvation, not structure. Single
  time-decimation D[c_j] stays a clock only for j ≤ 4. (Mild residue: at
  period 4 the germ hits 4 of 5 necklace classes, missing [0111]; not enough
  to build on.)

* **No spatial RG.** Block-spin (E[2j]⊕E[2j+1]) and even-cell subsystems are
  not deterministic local CAs at radius ≤ 5. The self-similarity, if any, is
  not a spatial coarse-graining.

Net: the real gains are algebraic reformulations (full-ANF-degree,
disjoint-convolution), not a nesting shortcut. The irreducible core remains
the disjoint-convolution parity of two adjacent clocks' complemented spectra.

## Modules

* `sitoolz.centerrandom` — center law, clocks, wedge index, left wedge,
  streams, pinned bits, `caveat_frame`.
* `sitoolz.edgelaw` — table extended rungs 30 → 48 via A094604 (verified on
  overlap); `edge(g)` valid to ~5.6×10^14 generations.
* `sitoolz.edgereach` — the compiled sweep.

## Still open

* The center law *locates* the center; it does not predict it past t = 113.
* The left-wedge table is empirical (≥6 periods of evidence per entry); a late
  period-doubling beyond the measured window could revise an entry.
* The right wedge is theorem-grade (T-function periodicity + A094604).
* Held-out statistics on the center remain at exactly 0.5 — there is no local
  statistical structure to exploit; whatever cracks this will be algebraic.
