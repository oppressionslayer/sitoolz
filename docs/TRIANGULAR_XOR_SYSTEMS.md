# The triangular register

*Run an unbounded XOR-driven system in a fixed register with no error term, and build clock
banks, rendezvous schedules and spectrum-specified drivers out of it. Rule 30 is one line of
one section.*

**Lars Rocha**

---

## The class

```
x_j(t+1) = x_j(t) XOR G(x_{j-1}, x_{j-2}, …, x_{j-d}),      x_j = 0 for j < 0
```

Arbitrary Boolean `G`, any depth, any initial condition. Three structural facts
follow, and none of them mentions any particular automaton.

**Autonomy.** Bit `j` reads only bits below it. So any prefix of the state is a
closed system — nothing above can reach down into it. A prefix is not an
approximation of the whole; it *is* the whole, restricted to a region that the
rest cannot touch.

**Exact width.** Therefore truncating at width `W` is exact for every bit below
`W`, at every generation, forever. An object that grows without bound runs in a
fixed register with no error term. This is why `E' = E ⊕ ((E≪1) | (E≪2))` can
be pushed to 2.2 trillion generations in four `uint64` words at 670 million
generations a second, and why the answer at the end is exact rather than
converged.

**Exact period.** Every coordinate is periodic in `t` — from `t = 0`, not
eventually — with period a power of two. Break triangularity by letting a
driver read a coordinate above itself and it dies immediately: periodicity
becomes eventual at best.

Those three together are the engineering content. What follows is what you can
build with them, and each item ships as a function in `sitoolz.txs`.

---

## 1. A clock bank off one shift-xor

A driver plus a width gives you a bank of signals whose periods are known
before you run anything:

```
x' = x XOR (x<<1) XOR (x<<3),  20 coordinates

      2^0   [0]
      2^1   [1]
      2^2   [2, 3]
      2^3   [4, 5, 6, 7]
      2^4   [8, 9, 10, 11, 12, 13, 14, 15]
      2^5   [16, 17, 18, 19]
```

The rows are the interesting part. Coordinates sharing a period are **distinct
waveforms of the same length**, produced by the same two shifts and one XOR, in
the same register. A binary counter gives you one signal per period; this gives
you a family per period, and the family sizes are the plateau structure of the
spectrum.

```python
txs.clock_bank((1, 3), 20)
```

---

## 2. The rendezvous

Where it holds, every coordinate whose period divides `2^v` reads **1** at
every `t ≡ −1 mod 2^v`. The whole bank re-synchronises, at times known in
advance, for no extra hardware — and the number of signals in phase at time `t`
is a function of `v₂(t+1)` alone.

That last sentence is rule 30's edge law, read as a scheduling primitive rather
than as a fact about a cellular automaton. The run of live cells at the edge of
the light cone *is* the count of signals currently in phase.

```
a OR b   (rule 30's edge)   holds    in phase at rung v: 1, 3, 4, 6, 7, 9, 15, 16, 24
a        (pure shift)       holds    in phase at rung v: 1, 2, 4, 8, 16, 24, 24, …
a + b + c                   fails    12 of 24 live coordinates
```

**Measured, over the 128 constant-free 3-input drivers: it holds on the live
coordinates of 76 of them.** Rule 30's `OR` is one. No simple test predicts
which — `G(1,1,1) = 1` agrees on only 70 of 104 — so `txs.verify_rendezvous`
checks a driver rather than assuming. Finding the invariant that separates the
76 from the 52 is an open question, and a good one.

---

## 3. Synthesis: choose a driver from a period spec

Because the spectrum of a linear driver is computable algebraically — least `k`
with `[z^j](F_T·Q_k) = 0` over every subset `T`, no generations evolved — you
can *search driver space against a specification* instead of simulating
candidates.

```python
txs.synthesize({20: 5, 8: 4})     # coordinate 20 of period 2^5, 8 of period 2^4
```

The answer is never a single driver. It is a class:

> **The spectrum is a coarse invariant.** 31 shift sets over shifts 1..5 give
> **13** distinct spectra; 127 sets over 1..7 give **41**. Every set whose
> minimum equals its gcd collapses onto the gcd alone — all 64 sets containing
> 1 have the identical spectrum. Past `min > gcd` the classes split further,
> so `(gcd, min)` does not determine it either.

Which is a free optimisation. Every member of the class has the same period
profile, so take the one with the fewest terms:

```
drivers meeting the spec: 16 of 31
cheapest member:    (1,)              1 XOR per generation
most expensive:     (1,2,3,4,5)       5 XOR per generation
```

Same bank, one fifth of the work per step.

---

## 4. Structured sparsity, by construction

The tame nonlinear drivers do not accelerate anything — they **filter**. A
subcube driver that misses the all-ones point freezes a sublattice of
coordinates at zero and runs the linear part's spectrum on the survivors:

```
b AND NOT c     stride 2     live: 0, 2, 4, 6, 8, 10, …
c AND NOT a     stride 3     live: 0, 3, 6, 9, 12, 15, …
```

So sparsity becomes a property of the update rule rather than a mask applied
afterwards. `txs.live_lattice` reports the stride.

---

## 5. The predictability dial

The tame set — linear drivers, plus subcube indicators missing the all-ones
point — is exactly the set whose period profile is closed-form. Everything else
has a spectrum that nobody can currently write down, rule 30's `OR` included.

```
a               closed form      (linear)
a AND NOT b     closed form      (subcube missing (1,1,1))
a OR b          no closed form   (not a subcube)
ab + ac         no closed form   (not a subcube)
```

That is a design knob rather than a limitation. Pick tame when you need the
periods known at design time. Pick outside it when you want a spectrum that is
not predictable from the rule — while still keeping the guarantee that every
period is a power of two and exact from `t = 0`, which is the part that never
depends on which driver you chose.

`txs.predictable(L, N)` answers it and says why.

---

## Where this class already lives

Triangular XOR-driven systems are **T-functions** — "T" for triangular — in the
sense of Klimov and Shamir (CHES 2002): mappings on n-bit words where output
bit `i` depends only on input bits `0..i`. Equivalently they are 1-Lipschitz
maps on the 2-adic integers. Addition, multiplication, AND, OR and XOR are all
T-functions, because carry propagates only upward — the same low-to-high flow
that puts an unbounded light-cone edge in a fixed register.

They were introduced to build stream ciphers with a *provable* period: the
triangular structure is what lets you prove a single cycle of length `2^n`
rather than hope for one. Multiword constructions ran over 5 Gbit/s in
software. The specific ciphers did not survive — distinguishing attacks and
linear weaknesses were found — while the structural claim was untouched.

That literature asks whether a period is **maximal**. The work here asks for
the **spectrum**: which coordinate has which period, which coordinates share
one, when they align. Both questions are about the same class of maps.
References in `AN_INVITATION.md`.

---

## The same property, in the other project

Triangularity is why a system resolves in dependency order without search. Bit
`j` depends only on what is already settled below it, so there is an order in
which every value falls out with nothing guessed — and the moment a driver
reads *above* itself, that order stops existing and you are searching.

That is the same property the L2 GF(2) rung leans on in the larsdoku solver,
where the null-space work over GF(2) is what lets a puzzle advance by algebra
instead of by branching. Both projects turn on the same structural fact from
opposite ends: here it says a light-cone edge fits in a register exactly; there
it says a constraint system resolves without a backtracker.

---

## Still open

1. **What separates the 76 from the 52?** Which drivers have the rendezvous.
   `G(1,1,1)` is not it.
2. **Prove the tame set.** The characterisation is 128/128 measured; the
   induction is the next step. The route: a subcube containing `(1,1,1)` fires exactly when the
   run below it is solid, which is exactly the rendezvous moment.
3. **Does the tame set survive at 4 inputs?** If it does, it is a statement
   about subcubes, not about three coordinates.
4. **Spectra for nonlinear drivers.** The exact rule is linear-only. The
   parity work in `MIXTURE_AND_SPECTRUM.md` §3 localises the obstruction to a
   single joint AND-parity term.

---

## Running it

```bash
python3 -m sitoolz.txs          # all five uses, demonstrated and verified
```

```python
from sitoolz import txs
txs.clock_bank((1,3), 20)          # periods and the families sharing them
txs.verify_rendezvous(3, 1)        # does this driver re-synchronise?
txs.synthesize({20: 5, 8: 4})      # drivers meeting a spec, nothing evolved
txs.cheapest(_)                    # the member of the class that costs least
txs.live_lattice(2, 4)             # which coordinates stay alive
txs.predictable(3, 1)              # closed-form spectrum, or not
txs.Bank.linear((1,3), 32).verify_period(20, 5)
```

Stdlib only, plus `sitoolz.driver_lab`.
