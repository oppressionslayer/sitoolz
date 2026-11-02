# The column-53209 case

*Settle the decidable half of Rowland's closing paragraph by exhaustive scan — column 53209 is the unique first occurrence — and map the search floor for the other half.*

**Lars Rocha**

---

## The sentence

*Local nested structure in rule 30* (Complex Systems **16** (2006) 239–258)
ends, just before the acknowledgment, by identifying where a counterexample to
the "only one left side of rule 30" conjecture would first arise — and leaving
it conditional:

> *"Upon examining the long term behavior for several initial conditions, one is
> tempted to conjecture that among nontrivial rightful rows R with a white left
> tail, the eventual period of column m in {f^t R} is independent of R — that
> there is really only one 'left side' of rule 30... However, this conjecture is
> **likely** false.*
>
> *...It turns out that this case first occurs at column m = 53209: Column 53208
> is eventually white, and the period of column 53207 is (…), resulting in two
> possible periods that **might occur** as the eventual period of column 53209
> and providing a counterexample to the conjecture **(if in fact they do occur
> for some initial conditions)**. Moreover, each of these periods leads to two
> other possible periods (at columns 58288 and 72577 respectively), and one
> **surmises** that in fact there are infinitely many possible left sides of
> rule 30."*

Four hedges in one paragraph. **Two separable questions live in it.**

1. **Does the case arise at all, and is 53209 really the first column where it
   does?** Decidable from the single-cell initial condition. Settled below.
2. **Are the two candidate periods realised by some rightful initial condition?**
   A search over initial conditions. **Still open**; the search floor is below.

---

## Setup, and the index check that gates everything

Left-justified rule 30, range `[−2,0]`:

```
new_i = c_(i-2) XOR (c_(i-1) OR c_i)
```

started from `I`, a single black cell.

His column numbering is not ours, and getting that wrong would make every
number below meaningless. So the index map is **validated against his own
published eventual-period list before the scan is allowed to run**:

```
his  : 1,1,1,2,1,2,2,1,4,1,4,4,4,4,4,4,4,4,4,4,4,4,4,2,4
ours : 1,1,1,2,1,2,2,1,4,1,4,4,4,4,4,4,4,4,4,4,4,4,4,2,4
match: True                                    (25 of 25)
```

His column `M` is our bit index `M − 1`. The scan refuses to proceed if that
check fails.

---

## Question 1 — settled

His three conditions, stated in his numbering:

- **(a)** column `m−1` is eventually white
- **(b)** `period(m)` is **not** greater than `period(m−2)`
- **(c)** the period of column `m−2` contains an **even** number of black cells

Measured, 131,072 generations, the last 8,192 recorded:

```
column 53207:  eventually white = False   period = 16
column 53208:  eventually white = True    period = 1     <- (a) holds
column 53209:  eventually white = False   period = 16

period of column 53207, 16 cells:
   [1,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0]     6 black -> EVEN  <- (c) holds
period(53209) = 16, period(53207) = 16, greater? False    <- (b) holds
```

Then the full scan — **every column from 3 to 53,260**, no sampling:

```
hits: [(53209, 16, 16, 6)]
first such column = 53209    Rowland states 53209    -> CONFIRMED
total hits in range: 1
```

> **His number is right, and it is the unique first occurrence.** One hit across
> the whole range, from a single black cell.

Runtime: **0.9 seconds.**

This **confirms the paper's analysis**. He worked out where the case must
first appear, by hand, in 2006, and got it exactly right.

---

## Question 2 — the search floor

His open half is whether either candidate period is actually *realised*. That is
a search over rightful initial rows, not a deeper simulation.

**Scale.** He verified his Conjecture 2 for central parts up to 10 cells, so
that scale is known to be tractable. The search covers the same.

**Shift-equivalence matters.** He states `a_{f^r R}(n) = a_R(n)` — a shifted row
is the *same* left side. An early run of this search reported "2 distinct
periods" and both were artifacts: one was a shift of `I`, whose column 53209
aligns with column 53208 of the unshifted run, which is white. Excluding shifts
and white columns is not a convenience; without it the search reports
counterexamples that are not counterexamples.

**Result:**

```
512 rightful rows tried (central parts to 10 cells, shifts excluded)
distinct eventual periods of column 53209: [16]
```

> **One period, every time.** Across every genuinely distinct initial condition
> at that scale, the "only one left side" conjecture **holds**.

**What this establishes.** It is not a refutation of his surmise, and it is not
a proof of the conjecture. He wrote *"one surmises that in fact there are
infinitely many possible left sides"* — this says that **if he is right, the
witnesses live outside central parts of ten cells.** Before this, *"if in fact
they do occur"* had no evidence in either direction. Now the nearest
neighbourhood is ruled out, and anyone continuing has a floor to start from
rather than a blank page.

---

## Still open

He names two further columns: *"each of these periods leads to two other
possible periods (at columns 58288 and 72577 respectively)."* Neither has been
examined here. Both are reachable with the same scan — the evolution is a few
seconds — and either could carry a witness that column 53209 does not.

Widening the central part is the other direction. Cost grows as `2^(w−1)` in the
number of rows and each row is one full evolution, so width 14 is roughly
8,000 rows at about 0.4 s each — an hour, not a project.

---

## Method

Question 1 is a single evolution of left-justified rule 30 from one black
cell for 131,072 generations, recording the last 8,192 rows, then a scan of
every column from 3 to 53,260 for conditions (a)–(c), gated by the 25-term
index check above. Question 2 repeats the evolution for every rightful row
with a central part of at most 10 cells, shifts excluded, and collects the
eventual period of column 53209.

Rowland, *Local nested structure in rule 30*, Complex Systems **16** (2006)
239–258. The conjecture, the column, the two candidate periods, and the
surmise are all his.
