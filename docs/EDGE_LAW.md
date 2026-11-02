# The edge law is a period spectrum

*Read rule 30's right light-cone edge straight off a lookup table, and see why it only sees the 2-adic valuation.*

**Lars Rocha**

---

## The observation being explained

`edgelaw` records that the run of live cells at the leading edge of rule
30's light cone is a function of `v2(g+1)` alone — the odd part of `g+1`
does not matter. Generations 255, 767, 1279 and 1791 all give 24, because
`g+1` is `2^8`, `3*2^8`, `5*2^8` and `7*2^8`. Verified single-valued to
generation 16,777,215, and stored as a measured lookup table.

The module states the law; this note gives the mechanism. The conclusion
is that the table is not a table of measurements at all — it is a
**cumulative period spectrum**, and the `v2` dependence is a consequence
of periodicity rather than a separate fact.

The edge in question is the **right** edge of rule 30's light cone. Derive
both edges from rule 30 directly, in their co-moving frames, and they come
out differently:

```
LEFT   D'[j] = D[j-2] ^ (D[j-1] | D[j])    reads its own coordinate; not triangular; only eventually periodic
RIGHT  R'[j] = R[j] ^ (R[j-1] | R[j-2])    edgelaw's recurrence
```

Checked over 2999 generations of real rule 30: the right-edge run matches
`edgelaw.edge(g)` with 0 mismatches, and the left-edge run mismatches 2998
times. OEIS agrees — A094603 says "rightmost", A094605 says "from the
right". Everything below is about the right edge, and the diagonals it
speaks of are the diagonals parallel to that edge.

---

## Setup

The edge evolves in the frame co-moving with the light cone:

```
E' = E ^ ((E << 1) | (E << 2))          E[0] = 1
```

equivalently `E'[j] = E[j] XOR (E[j-1] OR E[j-2])`.

Two properties matter, and both are structural rather than observed:

**The edge is autonomous.** Cells outside the cone are 0 forever, so they
never enter the update. The edge is a closed system; it is not an
approximation to the full automaton, it is exactly the automaton
restricted to a region that nothing else can reach.

**Information flows strictly low-to-high.** Bit `j` depends only on bits
`j-1` and `j-2`. So truncating at any width `W` is exact for every bit
below `W`, at any generation. This is why the whole thing is computable
in a 128-bit window forever, and why `edgelaw.measure` is O(1) per
generation.

---

## C1 — every bit is exactly periodic, with a power-of-two period

Not eventually periodic, not approximately. Exactly, from generation 0.

```
b     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
p(b)  0  1  1  2  3  3  4  5  5  6  6  6  6  6  6  7  8  8  8  8

b    20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38
p(b)  8  8  8  8  9 10 10 11 11 12 12 12 12 12 13 13 14 15 15
```

where bit `b` has period exactly `2^p(b)` in `g`. Resolved for 39 bits,
to period `2^15`, over 524,288 generations.

### Proof

Induction on the bit index.

*Base.* Bit 0: `E'[0] = E[0] XOR (E[-1] OR E[-2]) = E[0]`, so bit 0 is
constant 1 — period `2^0`. Bit 1: `E'[1] = E[1] XOR E[0] = E[1] XOR 1`,
so it flips every generation — period `2^1`.

*Step.* Fix `b >= 2` and suppose bits `b-1` and `b-2` are periodic with
periods `2^{p(b-1)}` and `2^{p(b-2)}`. Write the **driver**

```
d(g) = bit_{b-1}(g) OR bit_{b-2}(g)
```

which is periodic with period `P = 2^max(p(b-1), p(b-2))`, since an OR of
two signals is periodic with period the lcm of theirs, and an lcm of
powers of two is the larger of them.

The update is `bit_b(g+1) = bit_b(g) XOR d(g)`, so

```
bit_b(g) = bit_b(0) XOR (d(0) XOR d(1) XOR ... XOR d(g-1))
```

— bit `b` is the **running XOR-sum** of its driver. Let `S` be the weight
of `d` over one period, mod 2. Then:

- `S = 0` — the sum returns to where it started after `P` steps, so
  `bit_b` has period `P`.
- `S = 1` — the sum is inverted after `P` steps and restored after `2P`,
  so `bit_b` has period `2P`.

Either way the period is a power of two. ∎

The OR is the only nonlinearity in rule 30, and this argument routes
around it: the OR is used to *define* `d`, and then never again. What
carries the induction is that XOR-summation preserves power-of-two
periodicity regardless of what the summand is.

---

## C2 — where the period doubles

The proof above also says exactly *when* the period doubles: precisely
when the driver has odd weight over one of its own periods. That is the
`S = 1` branch, and it is the only empirical content left in the whole
problem.

Verified: 32 consecutive bits, 0 disagreements between "driver weight is
odd" and "period doubled here".  (This claim is checked in the standard
battery only; `deep()` skips it, since the driver-weight test needs the
unpacked columns.)

So the sequence of doubling positions — bits 1, 3, 4, 6, 7, 9, 15, 16,
24, 25, 27, 29, 34, 36, ... — is the sequence of bits whose driver has
odd weight. Characterising it is the open problem; see "Still open".

---

## C3 — bits in phase at `g = -1 mod 2^v`

A bit with period `2^p` is at the same point in its cycle at every
generation `g = -1 (mod 2^p)`. The observed fact is that this point is
always the **1**:

> bit `b` is 1 at every `g = -1 (mod 2^{p(b)})`.

Verified for all 39 resolved bits, 0 failures.

This is the hinge. It says the low bits are not merely periodic but
*synchronised*, all pointing the same way at the same moments — and the
moments are exactly the generations where `g+1` is divisible by a large
power of two.

---

## C4 — and where it stops

The complement: at level `v`, the first bit whose period exceeds `2^v` is
**0** there, so the run terminates exactly at that bit.

Verified for 15 levels, 0 failures.

---

## C5 — the table is the spectrum

Putting C3 and C4 together:

```
r(v)  =  #{ b : period(b) divides 2^v }
```

Checked against every entry of `edgelaw.TABLE` in range:

```
 v    TABLE   spectrum count
 0        1        1   ok
 1        3        3   ok
 2        4        4   ok
 3        6        6   ok
 4        7        7   ok
 5        9        9   ok
 6       15       15   ok
 7       16       16   ok
 8       24       24   ok
 9       25       25   ok
10       27       27   ok
11       29       29   ok
12       34       34   ok
13       36       36   ok
14       37       37   ok
```

15 levels, 0 mismatches.

---

## The plateaus — how many bits share a wheel

The period spectrum is flat in stretches: several consecutive bits can
carry the same period. Those plateaus are exactly the deltas of the
table,

```
run(v)  =  #{ b : p(b) = v }  =  TABLE[v] - TABLE[v-1]
```

so a long plateau *is* a big jump in the edge law. When eight bits share
period `2^8`, they all come into phase together and level 8 gains eight
cells at once.

```
 v : run   bits
 0 :   1   0..0    #
 1 :   2   1..2    ##
 2 :   1   3..3    #
 3 :   2   4..5    ##
 4 :   1   6..6    #
 5 :   2   7..8    ##
 6 :   6   9..14   ######
 7 :   1   15..15  #
 8 :   8   16..23  ########
 9 :   1   24..24  #
10 :   2   25..26  ##
11 :   2   27..28  ##
12 :   5   29..33  #####
13 :   2   34..35  ##
14 :   1   36..36  #
15 :   2   37..38  ##
16 :   2   39..40  ##
17 :   2   41..42  ##
18 :   5   43..47  #####
19 :   1   48..48  #
20 :   2   49..50  ##
21 :   3   51..53  ###
22 :   1   54..54  #
23 :   3   55..57  ###
24 :   2   58..59  ##
25 :   3   60..62  ###
26 :   1   63..63  #
27 :   2   64..65  ##
28 :   3   66..68  ###
29 :   1   69..69  #
30 :   2   70..71  ##
```

Distribution over the 31 known levels:

| run length | 1 | 2 | 3 | 5 | 6 | 8 |
|---|---|---|---|---|---|---|
| count | 10 | 13 | 4 | 2 | 1 | 1 |

Three observations, each still unexplained:

- **The maximum is 8, at v=8** *(within the first 31 rungs)*.
- **The monsters are early.** Every run above 3 happens at `v <= 18`.
  From v=20 on it is only 1s, 2s and 3s, and the growth of `r(v)`
  settles to roughly 2.4 cells per level.
- **Long plateaus are flanked by singletons.** The 8 at v=8 has `run(7)=1`
  and `run(9)=1` on either side; the 6 at v=6 has `run(7)=1` after it.
  In the odd-weight language of C2: a long stretch of even-weight drivers
  is bracketed by odd-weight ones.

`edgeframe.bits()` makes the last point mechanical: `driver_odd` is True
at exactly one bit per plateau — the first — and False for every bit
after it. So a plateau of length L is one odd-weight driver followed by
`L-1` even-weight ones, and the plateau-length question is "how long do
the even runs get".

### The gaps were nothing

Within the first 31 rungs the plateau lengths hit 1, 2, 3, 5, 6, 8 and
**skipped 4 and 7**, which looked like structure. It is not. A094604
continues to v=46, and both gaps fill:

```
plateaus, v = 0..46
1 2 1 2 1 2 6 1 8 1 2 2 5 2 1 2 2 2 5 1 2 3 1 3 2 3 1 2 3 1 2 2 3 2 1 2 2 2 4 1 2 7 3 1 2 2 3
                                                                            ^           ^
                                                                        4 at v=38   7 at v=41

length      1    2    3    4    5    6    7    8
count      13   21    7    1    2    1    1    1
```

Small-sample noise, and 31 draws was not enough to see it. Mean plateau
over v>=20 is 2.30.

`run(6) = 6` and `run(8) = 8` likewise: no third instance in 47 rungs, and
none is possible much higher, since plateau lengths stay in single digits
while v grows without bound. Coincidence.

Whether the plateau lengths are bounded at all is still open, and is the
same question as C2 from the other side.

---

## Why `v2` — and why the odd part is invisible

Now the original question answers itself.

`v2(g+1) = v` means `g+1 = 2^v * odd`, which means

```
g  =  -1  (mod 2^v)
```

Every bit whose period divides `2^v` is therefore in phase, and by C3 in
phase means showing 1. They form a solid run. The first bit with a longer
period has not completed its cycle, is out of phase, and by C4 shows 0.
The run stops there.

The odd part of `g+1` is irrelevant because it lives entirely **above the
modulus that any of those bits can see**. A bit with period `2^v` cannot
distinguish `2^v` from `3*2^v` from `5*2^v` — they are the same residue.
Generations 255, 767, 1279 and 1791 do not merely happen to give 24; the
bits involved are physically incapable of telling them apart.

**The edge is not computing anything. It is an odometer, and the run
length is how many wheels have just rolled over.**

That is also why `v2` and not some other statistic of `g`: `v2(g+1)` is
the carry-chain length when a binary counter increments, and the edge run
*is* a carry chain.

---

## Prior art — this is a known object

Checked against OEIS after the fact. The whole structure is catalogued,
and has been since 2004.

| ours | OEIS | their description |
|---|---|---|
| `edgelaw.TABLE[v]`, the edge run at rung v | **A094604** | largest run of rightmost black cells up to row `2^n` |
| the edge run at generation g | **A094603** | maximal run of rightmost black cells in row n |
| the period `2^p(b)` of bit b | **A094605** | *the period of the n-th diagonal, from the right, of rule 30* |
| the spectrum `p(b)` | **A094606** | base-2 logarithm of A094605 |

Agreement is exact: all 31 entries of `TABLE` against A094604, and all 39
resolved terms of `p(b)` against A094606.

So "every bit of the edge is periodic with a power-of-two period" is the
statement that rule 30's right-side **diagonals** are periodic — the same
object seen along a different axis. A bit at fixed offset from the moving
cone edge *is* a diagonal of the static triangle. The co-moving frame is
what makes the diagonals look like columns.

The lineage:

- **Wolfram**, *A New Kind of Science* (2002), p. 871 — observes the
  period doubling of the diagonals on the **left** side of rule 30.
- **Eric Rowland**, *Local nested structure in rule 30*, Complex Systems
  **16** (2006) 239–258 — proves Wolfram's left-side doubling criterion
  (his Proposition 2), inside a general theory of local nested structure
  for right-bijective rules. The mechanism is that at row `2^n` a region of
  the initial condition reappears and the automaton locally begins again,
  because rule 30 is reversible in time given that one half of each row is
  white. His **Lemma 2** covers the right edge's periodicity: for a
  right-bijective k-color range `[-d,0]` rule, column `m` is periodic with
  period dividing `l(k)·lcm(lower periods)`, `l(k) = lcm(1..k)`, and the
  periodicity is immediate rather than eventual. For `k=2` that is `2·lcm`,
  i.e. powers of two.
- **Eric Rowland**, *The Hidden Structure of Rule 30* (video, 2026), from
  07:42.

**One class, two names.** For `k=2`, `f(c,x)` bijective in `x` means
`f(c,x) = x XOR g(c)`, so **right-bijective and triangular XOR-driven are
the same class.** C1 states exact-from-zero periodicity for that class
and derives it from the recurrence; Rowland's Lemma 2 reaches the same
periodicity by a different route. The two derivations are laid side by
side in "Two routes to the same periodicity".

**What this note adds.** Lemma 2 gives a divisibility bound; it does not
say *when* the period doubles. Wolfram's and Rowland's doubling criterion
(Proposition 2) is a theorem about the other edge. C2 supplies the
doubling criterion for the right edge, and C3/C4 supply the phase lemmas
that turn it into the closed-form `v2` law. The `v2(g+1)` statement is
also sharper than what OEIS records: A094604's comment says only that
*new* values appear at positions `2^m` — a statement about records.
`edgelaw` asserts the stronger thing: the run at **every** generation is
`TABLE[v2(g+1)]`, with the odd part of `g+1` irrelevant. A094603 carries
that form as a conjecture, and Rowland's Section 3 records it as an
observed structure with eight terms of the table; the proof is in
[PROOF_V2_LAW.md](PROOF_V2_LAW.md).

**One nice confirmation.** A094605 carries the note *"One 64 term is
missing in Wolfram 2002 (p. 871). (There are six.)"* — the plateau at
v=6 has length six, and NKS printed five. Our independent measurement
gives six. We agree with the correction, not with the book.

---

## Two routes to the same periodicity

Both proofs are laid side by side here, with what each one buys. They are
not the same argument, and they do not prove quite the same statement.

### The two setups

**Rowland** works in the full two-dimensional pattern. He takes the
left-justified form of a range-`[-d,0]` rule — `f(c₋₂, c₋₁, x)` for rule
30 — and studies its *columns*, which are the diagonals of the usual
triangle. For right-bijective rules, Lemma 2 gives their periodicity
immediately, with a divisibility bound on the period. Proposition 2, his
doubling *criterion*, concerns the left side of rule 30, where Wolfram's
white-stripe observation lives.

**Ours** works in the one-dimensional frame that moves with the right
light-cone edge: `E'[j] = E[j] XOR (E[j-1] OR E[j-2])`, started from a
single cell. Bits of `E` are periodic **exactly, from generation zero**,
and periodicity and the doubling criterion come out of the same four-line
induction.

### The two criteria

Wolfram's observation, NKS p. 871, which Rowland's Proposition 2 proves
for the left side:

> Each period doubling occurs exactly when a diagonal becomes a white stripe, and
> the diagonal to its left has an odd number of black cells in each repeating block.

Ours (C2), for the right edge:

> The period doubles at bit `b` exactly when the driver `bit(b-1) OR bit(b-2)` has odd
> weight over one period.

The two have the same shape: where a white-stripe condition holds, the OR
collapses to a single column, and "odd weight of the driver" is literally
"odd number of black cells in the block." On the right edge the
white-stripe condition never occurs — bit 0 is constantly 1, and no bit of
the run is ever white — yet the periods double, and the odd-weight test
predicts every one of them (32 bits checked, 0 disagreements). The two
criteria are statements about different edges; C2 is the one that applies
here.

The enumerations differ accordingly. Rowland prints his eventual period
lengths for the left side, from the leftmost nonwhite diagonal, as

```
1, 1, 1, 2, 1, 2, 2, 1, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 4, ...
```

which is **not** A094605, the right-side sequence we match exactly:

```
1, 2, 2, 4, 8, 8, 16, 32, 32, 64, 64, 64, 64, 64, 64, 128, 256, ...
```

Two edges, two sequences. A094605 cites Rowland and NKS p. 871 as its
references.

### What each proof needs

| | Rowland 2006 | this note |
|---|---|---|
| machinery | right bijectivity, Jen's Thm 4, long-term-effect lemma | the recurrence, and parity |
| periodicity | immediate (Lemma 2), eventual in general (Jen) | exact, from t=0 |
| external results | yes (Jen) | none |
| generality | any right-bijective k-color range-`[-d,0]` rule, any rightful initial row | rule 30, single-cell start, cone edge |
| also proves | local nested structure; the "begin again" phenomenon at row 2ⁿ; the left-side doubling criterion | the right-side doubling criterion; the `v2` law |
| length | a section of a 20-page paper | four lines |

Rowland's is a general theorem about a class of automata, embedded in a
theory of where local nested structure comes from, and it explains *why* —
the pattern reconstructs a copy of its own initial condition at row 2ⁿ
because the rule is reversible when one half of each row is white.

This note's is a local computation in a frame where the system closes on
itself. It assumes nothing, cites nothing, gets periodicity from
generation zero, fits in a docstring next to the code it justifies, and
ends in a constant-time law.

The two buy different things.

### What the two routes do computationally

Rowland's paper is not a computational one and ours is nothing else.

**From the co-moving frame we get a memory bound.** Bit `j` reads only
`j-1` and `j-2`, so information flows strictly upward and truncation at
any width is exact below the cut. The edge of a light cone that grows
forever is trackable in a fixed 128-bit register at any generation
whatsoever. `edgereach` does ~900 million generations per second on three
`uint64` words.

**From the spectrum we get constant-time lookup.** `edge(g)` is `v₂(g+1)`
and a table index — no simulation at any size. That is the whole point of
`edgelaw`, and it is what the spectral reading justifies: the table stops
being measurements and becomes the cumulative period spectrum.

**From Rowland's route you get a different kind of shortcut** — the local
reducibility itself. If the pattern re-seeds at row 2ⁿ, you can restart
the local evaluation there instead of running through. His Conjecture 3
goes further and proposes the associated integer sequences as a *detector*
of reducibility in automata that show no global structure. That is a
research program; ours is a lookup table.

### What transfers

The transferable parts are the techniques, and they differ:

- **From Rowland:** bijectivity is what makes information survive, and
  information surviving is what makes structure. If you want to know
  whether a system re-seeds itself, ask what it is reversible with respect
  to.
- **From here:** find the sub-system that cannot see the chaos; check
  which way information flows before sizing your buffer; and quarantine
  the nonlinearity into a single parity bit instead of fighting it.

The second list came out of writing the code rather than out of the
mathematics.

---

## What this changes

**The law is no longer empirical.** The `v2`-only dependence is a theorem
given C1 and C3, and C1 is proved above. The verification to 16,777,215
was never testing whether the law holds — it was sampling the period
spectrum.

**`edge()` need not raise past `VMAX`.** The implementation refuses to
extrapolate beyond the measured table, which was the right call under the
empirical reading. Under the spectral reading, extending the table means
resolving more of `p(b)`, which needs `2^{v+1}` generations in a
fixed-width window rather than a full simulation — cheap, and the
narrow-window observation above makes it cheaper still.

**The open question got smaller and better-shaped.** It is no longer
"prove the table". It is:

> Characterise `p(b)` — equivalently, for which bits does the driver
> `bit_{b-1} OR bit_{b-2}` have odd weight over one period?

That is a question about one bit at a time, with a finite test, instead of
a question about the whole light cone.

---

## Reproducing

```python
from sitoolz import edgeproof
edgeproof.run()     # C1-C5, ~1 min
edgeproof.deep()    # C1, C3-C5 at depth, ~3 min
```

Stdlib only. `run()` defaults to 65,536 generations, 48 bits, periods
trusted to `2^12`. `deep()` goes to 524,288 generations and `2^15` by
storing columns as byte strings, so each period test is one `bytes`
comparison instead of a Python loop. Both take their constants as
arguments; raising them is the extension path.

Or read it down as frames:

```python
from sitoolz import edgeframe
edgeframe.bits()        # per bit: period, plateau, in_phase, driver_odd
edgeframe.levels()      # per rung: run, plateau, spectrum count, agreement
edgeframe.ladder(1,40)  # per generation: v2, lookup vs measurement
edgeframe.plateaus(3)   # the flat stretches, longest first
edgeframe.show()        # the same picture without pandas
```

## Modules

- `sitoolz.edgelaw` — the law as lookup
- `sitoolz.edgeproof` — the verification battery for this note
- `sitoolz.edgeframe` — the same material as frames, for reading down the rows

## Status of each claim

| | claim | status |
|---|---|---|
| C1 | every bit exactly periodic, power-of-two period | **proved** (induction), verified 39 bits to period 2^15 |
| C2 | period doubles iff driver has odd weight | proved as part of C1; the odd-weight positions remain uncharacterised |
| C3 | in-phase bits show 1 at `g = -1 mod 2^p` | verified 39 bits, 0 failures; proved in [PROOF_V2_LAW.md](PROOF_V2_LAW.md) |
| C4 | first out-of-phase bit shows 0 | verified 15 levels, 0 failures; proved in [PROOF_V2_LAW.md](PROOF_V2_LAW.md) |
| C5 | `TABLE[v] = #{b : period divides 2^v}` | verified 15 levels, 0 mismatches; follows from C3 + C4 |

## Still open

- **The doubling positions.** Which bits have an odd-weight driver — the
  sequence 1, 3, 4, 6, 7, 9, 15, 16, 24, ... — is measured, not derived.
  Characterising it is the whole remaining content of the edge law.
- **Plateau lengths.** Whether the even-weight runs between doublings are
  bounded. This is C2 seen from the other side.
- **The driver families.** Which drivers in the triangular XOR-driven
  class give a computable spectrum and which do not. See
  [DRIVER_SPACE.md](DRIVER_SPACE.md).
