# An Invitation — the doubling has nothing to do with rule 30

*Take the period-doubling theorem out of rule 30, state it for every triangular XOR-driven system, and work on the five questions it leaves open.*

**Lars Rocha**

---

## The one-paragraph version

Rule 30's right light-cone edge is a binary odometer: the run of live cells
at the edge is a function of `v2(g+1)` and nothing else, because each cell
is a clock with a power-of-two period, and `v2(g+1)` counts how many of
them have just come into phase. Wolfram observed period doubling on the
left side of rule 30 in 2002; Rowland proved that criterion in 2006, and
his Lemma 2 gives immediate power-of-two periodicity on the right side.

What this note is about is that **the mechanism does not use rule 30**. It
does not use the OR. It does not use cellular automata. It is a property of
triangular XOR-driven systems over GF(2), it holds for arbitrary Boolean
drivers, and the only thing rule 30 contributes is *which* positions double
— a single 0/1 sequence about which almost nothing is known.

That sequence is the ask.

---

## What is proved

Let `x_0, x_1, x_2, ...` be bits evolving in lock-step by

```
x_j(t+1)  =  x_j(t)  XOR  G_j( x_{j-1}(t), x_{j-2}(t), ..., x_{j-d}(t) )
```

with `x_j = 0` for `j < 0`, arbitrary Boolean functions `G_j`, arbitrary depth
`d`, arbitrary initial condition. Call this a **triangular XOR-driven
system**: every coordinate is updated by XOR-ing in some function of
coordinates strictly below it. For two colours this is exactly Rowland's
right-bijective class — `f(c,x)` bijective in `x` means `f(c,x) = x XOR g(c)`.

> **Theorem.** Every coordinate of a triangular XOR-driven system is exactly
> periodic in `t`, with period a power of two.

*Proof.* Induct on `j`. Coordinate 0 has no arguments, so `G_0` is a constant
and `x_0` has period 1 or 2. For `j > 0`, write `d_j(t) = G_j(...)`, the
**driver**. By hypothesis the coordinates it reads are periodic with
power-of-two periods, so `d_j` is periodic with period `P`, the largest of
them — an lcm of powers of two is the largest of them. The update telescopes:

```
x_j(t)  =  x_j(0)  XOR  d_j(0) XOR d_j(1) XOR ... XOR d_j(t-1)
```

so `x_j` is the running XOR-sum of `d_j`. Let `S` be the parity of `d_j`'s
weight over one period. If `S` is even, the sum returns home after `P` steps
and `x_j` has period dividing `P`. If `S` is odd, the sum is inverted after
`P` and restored after `2P`, so the period is exactly `2P`. Either way, a
power of two. ∎

Three things to notice.

**The nonlinearity is quarantined.** `G_j` is used exactly once, to name a
signal, and never again. Everything downstream is XOR-summation, which is
linear and does not care what it is summing. This is why the proof does not
have to know anything about `G_j` — not its algebraic degree, not its
symmetry, nothing.

**Periodicity is exact, not eventual.** From `t = 0`. There is no transient.

**Triangularity is load-bearing.** Let any `G_j` read `x_j` itself or a
coordinate above it and the property dies. In a small experiment — six random
systems where `G_j` reads one coordinate above — five failed, several with no
power-of-two-periodic coordinates at all.

### Rule 30's edge is one instance

In the frame co-moving with the right light-cone edge, rule 30 evolves by

```
E'[j] = E[j] XOR ( E[j-1] OR E[j-2] )
```

which is the theorem's shape with `G_j = OR` of the two below, for every `j`.
So the periods of rule 30's right-side diagonals are a corollary of a
statement that mentions neither rule 30 nor the OR.

### And the phase, which closes the law

Two further claims, both proved by the same induction:

- **In phase means 1.** `x_j(t) = 1` at every `t ≡ -1 (mod 2^{p(j)})`, where
  `2^{p(j)}` is `x_j`'s period. *Sketch:* evaluate the telescoped sum at
  `t = 2^{p(j)} - 1`. It collapses to `d_j(P-1)`, and `P-1 ≡ -1` modulo every
  lower period, so both arguments are 1 by induction, so the driver is 1.
- **Out of phase means 0.** For the first `j` with `p(j) = v+1`, `x_j = 0` at
  every `t ≡ 2^v - 1 (mod 2^{v+1})`. Same evaluation, one term shorter, and
  the odd-weight branch flips the answer.

Together with the theorem these give the edge law outright:

```
run at generation g  =  #{ j : period(x_j) divides 2^{v2(g+1)} }
```

— not verified-to-16-million, but derived. Every step checked numerically at
39 coordinates with zero failures before being believed. The full proof is
in [PROOF_V2_LAW.md](PROOF_V2_LAW.md).

---

## What is open

The theorem says the periods are powers of two. It says nothing about **which
positions double**, and that is the entire remaining content. Define the
**doubling indicator**

```
δ(j) = 1  if  p(j) = p(j-1) + 1,  else 0
```

equivalently: `δ(j) = 1` iff the driver `G_j(...)` has odd weight over one
period. For rule 30's edge, from `j = 1`:

```
1 0 1 1 0 1 1 0 1 0 0 0 0 0 1 1 0 0 0 0 0 0 0 1 1 0 1 0 1 0 0 0 0 1 0 1 1 0 1 0
1 0 1 0 0 0 0 1 1 0 1 0 0 1 1 0 0 1 0 1 0 0 1 1 0 1 0 0 1 1 0 1 0 1 0 0 1 0 1 1
```

(110 terms, derived from OEIS A094604.) Everything below is open.

**Q1 — Are the gaps bounded?** Runs of consecutive zeros in `δ` are the
plateaus of the period spectrum. Observed run lengths over 110 positions:

| zeros in a row | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| count | 13 | 21 | 6 | 1 | 2 | 1 | 1 | 1 |

Longest is 7, at `j = 16..23`. Is there an upper bound at all? Equivalently:
can arbitrarily many consecutive drivers have even weight?

**Q2 — Is `δ` 2-automatic?** The spectrum `p(j)` grows linearly and so cannot
be automatic, but `δ` is bounded, and every structure in sight is 2-adic. If
`δ` is 2-automatic, or `p` is 2-regular in the Allouche–Shallit sense, the
whole edge law becomes a finite object and the table stops being a table.
110 terms is far too few to decide this, and it is the question we would most
like answered. The 2-kernel evidence to 136 terms points against it; see
[DRIVER_SPACE.md](DRIVER_SPACE.md).

**Q3 — Density.** `δ` has density 0.418 over the first 110 positions. A fair
coin would give 0.5. Does the density converge, and to what?

**Q4 — Which drivers give which spectra?** The theorem holds for *all*
`G_j`, so driver-space is a laboratory. Running the same recurrence with each
of the sixteen 2-input Boolean drivers, 40 coordinates, periods resolved to
`2^11`, splits them cleanly in two:

| driver | coordinates resolved | spectrum |
|---|---|---|
| `XOR`, `XNOR`, `a`, `b`, `NOT a`, `NOT b`, `a AND NOT b`, `b AND NOT a`, `AND`, `0`, `1` | 40/40 | explicit |
| `OR`, `NOR`, `NAND`, `(NOT a) OR b`, `a OR (NOT b)` | 20–29/40 | outruns `2^11` |

At the transparent end the answer is closed form. For `G = XOR`:

```
p(j) = ceil(log2(j+1))        for j >= 1        (verified, 64 coordinates)
```

— plateau lengths exactly `1, 1, 2, 4, 8, 16, ...`. Nothing to ask.

At the other end sits `G = OR`, which is rule 30's edge, where plateau
lengths are `1, 2, 1, 2, 1, 2, 6, 1, 8, ...` and the density is 0.418.

The two ends are one bit apart in places. `AND` freezes — nothing ever
propagates. `NAND`, which is `AND` XOR the constant 1, is one of the opaque
five. **Adding a constant to the driver moves it from explicit to
intractable.**

So: **where in driver-space does the spectrum stop being computable, and what
is the invariant that separates the two columns above?** It is not linearity
(`AND` is nonlinear and transparent), not weight (`NOR` and `AND` both have
weight 1 and land on opposite sides), and not degeneracy alone. This is the
real question, and it is not a question about cellular automata. The
3-input answer — the tame set — is in [DRIVER_SPACE.md](DRIVER_SPACE.md).

**Q5 — Where is the theorem already written down?** Triangular XOR-driven
systems are **T-functions** (Klimov–Shamir, 2002), equivalently 1-Lipschitz
maps on the 2-adic integers, and that literature studies the same class
under another name. See "Prior art, part two" below. What we want is the
exact citation, and confirmation of the part that appears not to be covered:
that literature asks whether a period is **maximal**, and we are asking for
the **spectrum**.

---

## Prior art, part two: these are T-functions

After the rule 30 lineage was traced, the general class turned out to have a
literature too — in cryptography, not in cellular automata, which is why
neither Wolfram nor Rowland would have cited it.

**Klimov and Shamir, 2002.** A *T-function* ("T" for triangular) is a mapping
on n-bit words in which bit `i` of the output depends only on bits `0..i` of
the input. That is our condition, word for word. Addition, multiplication,
AND, OR and XOR are all T-functions, because carry propagates only upward —
the same low-to-high flow that lets rule 30's edge live in a fixed register.

They introduced them to build stream ciphers with a *provable* period: the
triangular structure is what lets you prove a single cycle of length `2^n`
rather than hope for one. Multiword constructions ran over 5 Gbit/s in
software, roughly an order of magnitude faster than RC4 at the time, and the
idea was later extended to block ciphers and hash functions. The specific
ciphers did not survive — distinguishing attacks and linear weaknesses were
found — but the structural claim is untouched.

**The mathematical home is 2-adic dynamics.** A T-function is precisely a
**1-Lipschitz map on the 2-adic integers**: depending only on lower bits *is*
being non-expanding in the 2-adic metric. Anashin and others have criteria for
when such maps are measure-preserving and ergodic on `Z_2`, with van der Put
and Mahler expansions as the working tools.

**And the oldest instance is a counter.** A ripple counter is
`x_j ^= AND(x_0, ..., x_{j-1})`, coordinate `j` has period `2^{j+1}`, and each
stage halves the frequency. Every ripple counter ever fabricated is a
triangular XOR-driven system with the AND driver. The odometer reading was not
a metaphor.

### What this leaves open

The important thing is what that literature asks. **T-function work asks
whether the period is maximal** — single cycle of length `2^n`, yes or no,
because that is what a keystream needs. **We are asking for the spectrum**:
the per-coordinate period profile `p(j)`, and which positions double. That is
a finer object than a yes/no, and the spectrum is not a question that
literature poses.

So Q5 splits:

- The *theorem* is almost certainly a corollary in the T-function or 2-adic
  literature. We would like the exact citation.
- Q1–Q4, which are about the **spectrum** rather than about maximality, may
  survive intact even after Q5 is settled.

And a consequence worth stating: rule 30's light-cone edge and a
Klimov–Shamir stream cipher are the same kind of object, differing only in the
driver. Theirs is chosen for maximal period. Rule 30's is `OR`.

---

## Why this might interest you

**If you work on automatic sequences.** Q2 is your question, and the object
is unusually clean: a 0/1 sequence generated by a parity condition on a
running XOR-sum, with a natural 2-adic scale built in. Rowland's own work
sits at exactly this intersection.

**If you work on GF(2) linear algebra.** A triangular XOR-driven system is a
unipotent map on a filtered GF(2) space. The theorem says every unipotent
orbit has 2-power period, which is unsurprising; the content is that *the
filtration degree of each coordinate* is what the doubling indicator measures.
That is a Schur-complement-shaped question, and we would like to know whether
the spectrum can be read off a triangular decomposition rather than simulated.

**If you work on carry propagation or shift registers.** `v2(n)` is the carry
chain length when a binary counter increments, and the edge run *is* a carry
chain. The general theorem says any triangular XOR-driven system is a
generalized counter with a nonstandard carry rule. The plateau lengths are how
far a carry propagates. Bounded carries and Q1 are the same question.

**If you work on cellular automata.** The result de-emphasizes the rule. Rule
30's right-side regularity is not a fact about rule 30; it is a fact about
triangular (right-bijective) updates, and rule 30 supplies only the driver.
Which of the elementary rules give a decidable spectrum is unexplored here.

**If you work on stream ciphers or T-functions.** You have the theorem
already, in the form "these constructions have provable period". What we are
asking is one level finer: not whether the cycle is maximal, but what the
per-coordinate period profile looks like when it *is not* — which is the
regime a cryptographer discards and we are working in. `OR` is the driver
rule 30 hands us, and its spectrum has density 0.393 with no closed form
yet found.

**If you work on 2-adic dynamics.** These are 1-Lipschitz maps on `Z_2`.
Measure-preservation and ergodicity are well studied; the spectrum `p(j)`
is not. Q2 asks whether it is 2-regular, which would make it a finite object
in exactly your language.

---

## Ledger

| claim | status |
|---|---|
| rule 30's left-side diagonals have doubling periods | **known** — Wolfram, NKS p. 871 (2002) |
| that left-side doubling criterion, proved | **known** — Rowland, Complex Systems 16 (2006), Proposition 2 |
| immediate power-of-two periodicity on the right edge | **known** — Rowland, Lemma 2 (divisibility bound) |
| the run-length sequences | **known** — OEIS A094603–A094606, since 2004 |
| exact-from-zero periodicity in the co-moving frame, from the recurrence | ours |
| the right-edge doubling criterion | ours |
| the general triangular XOR-driven theorem | ours; the same class is studied in the T-function literature — see *Prior art, part two* |
| the phase lemmas closing the `v2` law | ours, proved |
| everything in **What is open** | open |

We arrived at the rule-30 facts by measurement, independently, and traced the
prior art after proving the mechanism. That is recorded in
[THE_ODOMETER.md](THE_ODOMETER.md) in full, including the false leads that
died.

---

## The ask

Any of Q1–Q5. Q2 and Q5 are the ones a specialist could settle in an
afternoon.

If Q4 is as open as it looks, it is a program rather than a problem: a map
from Boolean drivers to period spectra, with rule 90 at the transparent end
and rule 30 somewhere out in the dark. [DRIVER_SPACE.md](DRIVER_SPACE.md)
is the first survey of that map.

---

## The code

Everything in this note reproduces from **one module**, `sitoolz.tri_xor`
— stdlib only.

```python
from sitoolz import tri_xor as tx

tx.check()          # the theorem and both phase lemmas, verified
tx.show_survey()    # all sixteen 2-input drivers
tx.spectrum(tx.OR)  # rule 30's edge: period exponents p(j)
tx.delta(tx.OR)     # the doubling indicator -- the open object
tx.edge_run(1791)   # 24, by valuation and table index, no simulation
```

```
$ python3 -m sitoolz.tri_xor
C1  34 coordinates resolved, all periods powers of two
C3  in-phase bits show 1        -- 0 failures
C4  out-of-phase bits show 0    -- 0 failures over 12 levels
law run(g) == table[v2(g+1)]    -- 0 mismatches to generation 4096
```

The rest is convenience and speed, not necessity:

```python
from sitoolz import edgelaw     # the law as lookup; edge(10**12) is instant
from sitoolz import edgeframe   # spectrum as pandas frames: bits(), levels(), plateaus()
from sitoolz import edgeproof   # the claims, checked at depth: run(), deep()
from sitoolz import edgereach   # compiled sweep, ~790M generations/sec (numba)
```

| where | what it is |
|---|---|
| `sitoolz.tri_xor` | **standalone** — the whole note, one stdlib module |
| [THE_ODOMETER.md](THE_ODOMETER.md) | the narrative, with the derivation and the dead ends |
| [EDGE_LAW.md](EDGE_LAW.md) | the technical note: claims, proofs, verification counts, prior art |
| [PROOF_V2_LAW.md](PROOF_V2_LAW.md) | the full proof of the `v2` law |
| [DRIVER_SPACE.md](DRIVER_SPACE.md) | the driver families and what is known about the spectrum |

## Sources

- OEIS [A094603](https://oeis.org/A094603) · [A094604](https://oeis.org/A094604) · [A094605](https://oeis.org/A094605) · [A094606](https://oeis.org/A094606)
- Eric Rowland, [*Local nested structure in rule 30*](https://ericrowland.github.io/papers/Local_nested_structure_in_rule_30.pdf), Complex Systems **16** (2006) 239–258
- Stephen Wolfram, [*A New Kind of Science*, p. 871](http://www.wolframscience.com/nksonline/page-871)
- A. Klimov and A. Shamir, *A New Class of Invertible Mappings* (CHES 2002); [*New Applications of T-Functions in Block Ciphers and Hash Functions*](https://link.springer.com/chapter/10.1007/11502760_2) (FSE 2005)
- [*Distinguishing Attacks on T-Functions*](https://link.springer.com/chapter/10.1007/11554868_2) (Mycrypt 2005)
- [T-function](https://en.wikipedia.org/wiki/T-function) · [T-functions revisited: new criteria for bijectivity/transitivity](https://arxiv.org/pdf/1111.3093)
- V. Anashin et al., [*Toward the ergodicity of p-adic 1-Lipschitz functions represented by the van der Put series*](https://arxiv.org/pdf/1210.5001)
