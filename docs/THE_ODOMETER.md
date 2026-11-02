# The Odometer — rule 30's edge, from measurement to law

*Turn a measured lookup table into a theorem, then read the edge of rule 30's light cone at generation 10¹² without simulating a step.*

**Lars Rocha**

---

## Where it starts

`edgelaw` holds a lookup table and a claim. The claim is that the run of
live cells at the leading edge of rule 30's light cone depends only on the
2-adic valuation of `g+1` — how many times two divides it — and on nothing
else. The odd part is irrelevant. Generations 255, 767, 1279 and 1791 all
give a run of 24, because `g+1` is `2^8`, `3·2^8`, `5·2^8` and `7·2^8`.

The table was measured, not derived. Thirty-one rows, verified
single-valued across millions of generations, and no reason given.

Rule 30 is the automaton people reach for when they want an example of
something that plainly has no shortcut. Its centre column has been used as
a random number generator. So a lookup table that answers a question about
it in constant time is either wrong or interesting.

The edge in question is the **right** edge of the light cone. Derived from
rule 30 directly, the two edges have different co-moving recurrences:

```
LEFT   D'[j] = D[j-2] ^ (D[j-1] | D[j])    reads its own coordinate; not triangular; only eventually periodic
RIGHT  R'[j] = R[j] ^ (R[j-1] | R[j-2])    edgelaw's recurrence
```

Over 2999 generations of real rule 30 the right-edge run matches
`edgelaw.edge(g)` with 0 mismatches; the left-edge run mismatches 2998
times. OEIS A094603 says "rightmost" and A094605 says "from the right",
and that is the edge this story is about.

---

## The frame

The edge is easiest to see in the reference frame that moves with it:

```
E' = E ^ ((E << 1) | (E << 2))          E starts at 1
```

which is `E'[j] = E[j] XOR (E[j-1] OR E[j-2])`. Two properties fall out
immediately, and both are structural rather than lucky.

**The edge is a closed system.** Cells outside the light cone are 0 and stay
0, so they never enter the update. This is not an approximation to rule 30
restricted to a region — it *is* rule 30, on a region nothing else can
reach.

**Information only flows upward.** Bit `j` reads bits `j-1` and `j-2` and
nothing above. So truncating the state at any width is exact for every bit
below the cut, at any generation whatsoever. The light cone nominally grows
one cell per step forever; you can follow its edge in a 128-bit register
until the heat death and never lose a bit that matters.

That second property is why any of this is checkable. It's also worth
sitting with for a second — an object famous for being computationally
irreducible has a boundary you can track in constant memory.

---

## Every bit is a clock

Here is the thing that turns the table from data into structure.

**Every bit of the edge is exactly periodic in the generation, and the
period is always a power of two.** Not eventually periodic. Exactly, from
generation zero.

```
bit    0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
period 1  2  2  4  8  8 16 32 32 64 64 64 64 64 64 128 256 ................. 512
```

The proof is four lines, and the pleasing part is what it does with the OR.

Rule 30's entire reputation rests on that OR. It is the only nonlinearity in
the rule; without it you'd have rule 90, Sierpinski, and a closed form.
Watch where it goes.

Induct on the bit index. Bit 0 never changes — its neighbours below don't
exist, so `E'[0] = E[0]`, constant 1, period 1. Bit 1 flips every step,
period 2. Now take any bit `b` and call

```
d(g)  =  bit_{b-1}(g)  OR  bit_{b-2}(g)
```

its **driver**. By hypothesis the two bits below are periodic with
power-of-two periods, so `d` is periodic with period `P`, the larger of the
two. And the update says

```
bit_b(g)  =  bit_b(0)  XOR  d(0) XOR d(1) XOR ... XOR d(g-1)
```

— bit `b` is nothing but the **running XOR-sum of its driver**. Let `S` be
the parity of `d`'s weight over one period. If `S` is even the sum returns
home after `P` steps and bit `b` has period `P`. If `S` is odd the sum comes
back inverted after `P` and true after `2P`, so the period doubles. Either
way it stays a power of two. Done.

The OR appears once, to define `d`, and never again. Everything after it is
XOR-summation, which is linear and preserves power-of-two periodicity no
matter what it is fed. The whole chaotic reputation of rule 30 gets
compressed, at the edge, into **one bit of information per bit position**:
is the driver's weight odd?

---

## Why the valuation

Now the original question answers itself, and the answer is a machine you
already own.

A bit with period `2^p` is at the same point in its cycle at every
generation `g ≡ -1 (mod 2^p)`. Measurement says that point is always the
**1**. So at any generation where

```
v2(g+1) = v      i.e.      g ≡ -1  (mod 2^v)
```

every bit whose period divides `2^v` is *in phase*, and in phase means
showing 1. They form a solid run. The first bit with a longer period hasn't
come around yet, shows 0, and the run stops there.

Which gives the law in one line:

```
run(v)  =  #{ bits whose period divides 2^v }
```

The table isn't a table of measurements. It's the **cumulative period
spectrum** of the edge, read off rung by rung. Checked against every entry
in range: exact.

And the odd part of `g+1` is invisible for a reason that isn't subtle at
all. A bit with period `2^v` **cannot distinguish** `2^v` from `3·2^v` from
`5·2^v` — same residue, same phase, same cell. Generations 255, 767, 1279
and 1791 don't happen to agree. The bits doing the work are physically
incapable of telling them apart.

> **The edge is not computing anything. It is an odometer, and the run
> length is how many wheels have just rolled over.**

That is also why the answer is `v2` rather than any other statistic of `g`.
The carry chain when a binary counter increments has length exactly
`v2(n)`. The edge run *is* a carry chain. Rule 30's boundary is counting,
and has been the whole time.

---

## The plateaus

Several consecutive bits can share a period. Those flat stretches are
exactly the jumps in the table — a plateau of length `L` means `L` bits come
into phase together and the run gains `L` cells at that rung.

```
 v : plateau
 6 :   6      ######
 7 :   1      #
 8 :   8      ########      <- eight bits, all period 2^8
 9 :   1      #
12 :   5      #####
18 :   5      #####
```

Reading it per-bit sharpens it further: the driver is odd-weight at exactly
one bit per plateau — the first — and even for every bit after. So a
plateau of length `L` is one odd driver trailed by `L-1` even ones, and the
only open question in the whole construction is **how long the even runs
get.**

---

## The false lead

Over the first 31 rungs the plateau lengths came out

```
length      1    2    3    4    5    6    7    8
count      10   13    4    0    2    1    0    1
```

Ten ones, thirteen twos, no fours, no sevens — but a five, a six and an
eight. Zero fours in 31 draws while an eight sits there is the shape of
something. The mean was drifting toward 2, which is what a fair coin would
give, while the distribution refused to be one.

It was nothing. The sequence continues to rung 46 and both gaps fill: a
**4 at v=38**, a **7 at v=41**.

```
1 2 1 2 1 2 6 1 8 1 2 2 5 2 1 2 2 2 5 1 2 3 1 3 2 3 1 2 3 1 2 2 3 2 1 2 2 2 4 1 2 7 3 1 2 2 3
                                                                          ^           ^
```

Both were then re-measured from scratch rather than read off the catalogue.
A compiled sweep to generation `2^41` — 2.2 trillion steps, 46 minutes at
790M generations/sec — reproduces rungs 0–30 exactly and independently
produces the 4 at v=38 and the 7 at v=41:

```
 v    run   plateau   status
30     72         2   ok
31     74         2   NEW
...
38     90         4   NEW     <- the "impossible" four
39     91         1   NEW
40     93         2   NEW
41    100         7   NEW     <- the "impossible" seven
```

The retraction stands on our own arithmetic, not on trust.

Thirty-one samples was not enough to see a 4, and the hole was noise wearing
a pattern's clothes. Same for the other coincidence worth mentioning:
`run(6) = 6` and `run(8) = 8`, plateau length equal to its own rung, twice,
at the two biggest plateaus — and never again in 47 rungs. It can't recur
much higher either, since plateau lengths stay in single digits while `v`
grows without bound.

Both were small-number effects. Worth writing down that they were checked
and died, because that is most of what checking is for.

---

## Then we looked it up

The sequence of run lengths is **OEIS A094604**. It has been there since
2004.

| ours | OEIS | theirs |
|---|---|---|
| the table, `run(v)` | A094604 | largest run of rightmost black cells up to row `2^n` |
| the run at generation `g` | A094603 | maximal run of rightmost black cells in row `n` |
| the period of bit `b` | **A094605** | **period of the n-th diagonal, from the right, of rule 30** |
| the spectrum, `log2` of it | A094606 | base-2 logarithm of A094605 |

Agreement is exact — all 31 table entries against A094604, all 39 measured
periods against A094606.

The translation is the part worth keeping. **A bit at fixed offset from the
moving cone edge is a diagonal of the static triangle.** Sitting in the
co-moving frame turns diagonals into columns, which is why they look like
clocks from here and like diagonals from the printed page. "Every bit of the
edge is periodic with a power-of-two period" and "the right-side diagonals
of rule 30 have doubling periods" are one sentence in two languages.

The lineage runs:

- **Wolfram**, *A New Kind of Science* (2002), p. 871 — observes the
  doubling on the left side of rule 30.
- **Eric Rowland**, *Local nested structure in rule 30*, Complex Systems
  **16** (2006) 239–258 — proves the left-side criterion (Proposition 2),
  and in Lemma 2 gives immediate periodicity with a power-of-two bound for
  the right-bijective class, which for two colours is exactly the
  triangular XOR-driven class the edge frame lives in. His mechanism is a
  different one: at row `2^n` a region of the initial condition reappears
  and the automaton locally begins again, because rule 30 is reversible in
  time provided one half of each row is white.
- **Eric Rowland**, *The Hidden Structure of Rule 30* (video, 2026), from
  07:42.

So the induction above proves the periodicity from the recurrence, arrived
at independently here, and adds the right-edge doubling criterion and the
phase lemmas that Lemma 2's divisibility bound does not contain. The
prior-art comparison is in [EDGE_LAW.md](EDGE_LAW.md).

---

## What survives

**The `v2` statement is sharper than the catalogue.** A094604's comment says
only that *new* values appear at positions `2^m` — a claim about records.
`edgelaw` claims the stronger thing: the run at **every** generation is
`table[v2(g+1)]`, odd part irrelevant. A094603 carries that as a
conjecture; [PROOF_V2_LAW.md](PROOF_V2_LAW.md) proves it.

**The proof is short and lives next to the code.** Rowland's argument goes
through time-reversibility. The XOR-running-sum argument is four lines, needs
nothing but the recurrence, and can sit in a docstring.

**The odometer reading is a way to hold it.** Carry chains and 2-adic
valuations are the same object, and the edge of rule 30's light cone is a
counter incrementing. That's a sentence you can keep.

**And one small pleasure.** A094605 carries the note *"One 64 term is
missing in Wolfram 2002 (p. 871). (There are six.)"* The plateau at rung 6
has length six; the book printed five. Our independent measurement says six.
We agree with the correction, not with the book — which is a nice place to
land after finding out you were second.

---

## Running it

```python
from sitoolz import edgelaw         # the law, as lookup
edgelaw.edge(10**12)                # instant at any size

from sitoolz import edgeframe       # the same material, as frames
edgeframe.bits()                    # per bit: period, plateau, phase, driver parity
edgeframe.levels()                  # per rung: run, plateau, spectrum, agreement
edgeframe.ladder(1, 40)             # per generation: v2, lookup vs measurement
edgeframe.plateaus(3)               # the flat stretches, longest first

from sitoolz import edgeproof       # the claims, checked
edgeproof.run()                     # ~1 min
edgeproof.deep()                    # ~3 min, 39 bits out to period 2^15

from sitoolz import edgereach       # push to new rungs, compiled
edgereach.compare(38)               # ~900M generations/sec under numba
```

Stdlib only, except `edgeframe` (pandas, soft) and `edgereach` (numba, with
a slow fallback).

---

## Modules

| module | what it is |
|---|---|
| `sitoolz.edgelaw` | the law as lookup |
| [EDGE_LAW.md](EDGE_LAW.md) | the technical note: claims C1–C5, proofs, verification counts, prior art |
| `sitoolz.edgeproof` | the verification battery for those claims |
| `sitoolz.edgeframe` | the spectrum as pandas frames |
| `sitoolz.edgereach` | compiled sweep to higher rungs |

---

## Sources

- OEIS [A094603](https://oeis.org/A094603) · [A094604](https://oeis.org/A094604) · [A094605](https://oeis.org/A094605) · [A094606](https://oeis.org/A094606)
- Eric Rowland, [*Local nested structure in rule 30*](https://ericrowland.github.io/papers/Local_nested_structure_in_rule_30.pdf), Complex Systems 16 (2006) 239–258
- Stephen Wolfram, [*A New Kind of Science*, p. 871](http://www.wolframscience.com/nksonline/page-871)
