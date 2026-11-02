"""txs -- triangular XOR-driven systems, as a design tool.

The class is

    x_j(t+1) = x_j(t) XOR G(x_{j-1}, x_{j-2}, x_{j-3}),      x_j = 0 for j < 0

and the three structural facts that make it useful are not about rule 30:

    AUTONOMY      bit j reads only bits below it, so any prefix of the state
                  is a closed system.  Nothing above can reach down into it.
    EXACT WIDTH   therefore truncating at width W is exact for every bit below
                  W, at every generation, forever.  An unbounded object runs
                  in a fixed register with no approximation.
    EXACT PERIOD  every coordinate is periodic from t = 0 -- not eventually --
                  with period a power of two.

What you get from that is a bank of signals whose periods are known at design
time, off one shift-xor update, in one register.  This file is the tooling for
picking the driver.

    from sitoolz import txs

    txs.clock_bank((1,), 16)            # what periods a driver gives you
    txs.tap((1,3), 5)                   # cheapest coordinate of period 2^5
    txs.synthesize({20: 5, 8: 3})       # drivers meeting a spec
    txs.predictable(3, 1)               # is the spectrum closed-form?
    txs.live_lattice(2, 4)              # which coordinates a driver keeps alive
    txs.Bank.linear((1,2), 24).run(64)  # actually run one

Depends on driver_lab.py (same directory) for the exact spectrum rule and the
subcube test.  Stdlib otherwise.
"""
from sitoolz import driver_lab as dl


# ------------------------------------------------------------------ spectrum
def spectrum(S, jmax=63):
    """Period exponent p(j) for the linear driver with shift set S.

    Exact, and nothing is evolved -- see MIXTURE_AND_SPECTRUM.md section 2.
    """
    return dl.spectrum_linear_exact(tuple(S), jmax)


def clock_bank(S, width=24):
    """Every coordinate's period, and the families that share one.

    Coordinates sharing a period are the useful part: they are distinct
    waveforms of the same length, off the same register, for free.
    """
    p = spectrum(S, width - 1)
    fam = {}
    for j, k in enumerate(p):
        if k is None:
            continue
        fam.setdefault(k, []).append(j)
    return {"p": p, "families": dict(sorted(fam.items()))}


def tap(S, k, jmax=255):
    """Lowest coordinate whose period is exactly 2^k, or None."""
    p = spectrum(S, jmax)
    for j, v in enumerate(p):
        if v == k:
            return j
    return None


def synthesize(spec, smax=5, jmax=127, limit=None):
    """Every shift set S ⊆ {1..smax} whose spectrum satisfies `spec`.

    `spec` maps coordinate -> required period exponent.  Because the spectrum
    is computed algebraically, this is a search over drivers with no system
    ever run.
    """
    out = []
    for mask in range(1, 1 << smax):
        S = tuple(s for s in range(1, smax + 1) if mask >> (s - 1) & 1)
        p = spectrum(S, jmax)
        if all(j <= jmax and p[j] == k for j, k in spec.items()):
            out.append(S)
            if limit and len(out) >= limit:
                break
    return out


def spectrum_classes(smax=5, jmax=63):
    """Group every shift set over shifts 1..smax by its spectrum.

    The spectrum is a coarse invariant: 31 sets over shifts 1..5 give 13
    distinct spectra, 127 sets over 1..7 give 41.  Every set whose minimum
    equals its gcd collapses onto the gcd -- all 64 sets containing 1 have the
    identical spectrum.  Past min > gcd the classes split further, so (gcd,
    min) alone does not determine it.

    The practical consequence: a period spec never picks out one driver, it
    picks out a class, and you are free to take the cheapest member.
    """
    cls = {}
    for mask in range(1, 1 << smax):
        S = tuple(s for s in range(1, smax + 1) if mask >> (s - 1) & 1)
        cls.setdefault(tuple(spectrum(S, jmax)), []).append(S)
    return list(cls.values())


def cost(S):
    """XORs per generation: one per shift."""
    return len(S)


def cheapest(sets):
    """The member of a class with the fewest shifts -- same spectrum, less work."""
    return min(sets, key=lambda S: (cost(S), S)) if sets else None


# ------------------------------------------------------- driver classification
def predictable(L, N):
    """Is this driver's spectrum closed-form -- tame -- or not?

    The tame set is the linear drivers plus the subcube indicators that miss
    the all-ones point (MIXTURE_AND_SPECTRUM.md section 1, 128/128).  Tame is
    the setting to pick when you want the period profile known at design time;
    the others are the setting to pick when you do not.
    """
    tt = dl.truth(L, N)
    sub = dl.is_subcube(tt)
    ones = dl.fires_on_ones(tt)
    tame = (tt == 0) or (N == 0) or (sub and not ones)
    why = ("zero driver" if tt == 0 else
           "linear" if N == 0 else
           "subcube missing (1,1,1)" if tame else
           "subcube containing (1,1,1)" if sub else
           "not a subcube")
    return {"tame": tame, "reason": why, "subcube": sub, "fires_on_ones": bool(ones)}


def live_lattice(L, N, seed=1, width=24, steps=1 << 12):
    """Which coordinates a driver keeps alive, and which it freezes at 0.

    The tame nonlinear drivers sparsify: they freeze a sublattice and run the
    linear part's spectrum on what is left.  That is a way to build structured
    sparsity into a generator by construction rather than by masking.
    """
    R = dl.rows(L, N, width, steps, seed)
    live = [j for j in range(width) if any((r >> j) & 1 for r in R)]
    step = None
    if len(live) > 2:
        gaps = {live[i + 1] - live[i] for i in range(len(live) - 1)}
        if len(gaps) == 1:
            step = gaps.pop()
    return {"live": live, "frozen": [j for j in range(width) if j not in live],
            "stride": step}


# ------------------------------------------------------------- the rendezvous
def verify_rendezvous(L, N, seed=1, width=24, steps=1 << 13, kmax=12):
    """Do all coordinates of period dividing 2^v read 1 at t = -1 mod 2^v?

    Where it holds, the bank re-synchronises at known times for free, and the
    number of signals in phase at time t is a function of v2(t+1) alone --
    which is exactly rule 30's edge law, read as a scheduling primitive.

    Measured over the 128 constant-free 3-input drivers: it holds on the live
    coordinates of 76 of them, rule 30's OR among them.  No simple test
    predicts which -- G(1,1,1) agrees on only 70 of 104.  So: check, per
    driver, rather than assume.
    """
    R = dl.rows(L, N, width, steps, seed)
    p = dl.spectrum(R, width, kmax)
    cols = [bytes((r >> j) & 1 for r in R) for j in range(width)]
    live = [j for j in range(width) if any(cols[j])]
    tested = bad = 0
    for j in live:
        if p[j] is None:
            continue
        st = 1 << p[j]
        tested += 1
        if any(cols[j][t] != 1 for t in range(st - 1, steps, st)):
            bad += 1
    return {"live": len(live), "tested": tested, "failures": bad,
            "holds": tested > 0 and bad == 0}


def in_phase_profile(L, N, seed=1, width=24, steps=1 << 13, kmax=12):
    """How many coordinates are in phase at each rung v -- the edge law shape."""
    R = dl.rows(L, N, width, steps, seed)
    p = dl.spectrum(R, width, kmax)
    return [sum(1 for v in p if v is not None and v <= k) for k in range(kmax + 1)]


# ----------------------------------------------------------------- running one
class Bank:
    """A running triangular XOR system in one register.

    Build it explicitly -- `Bank.linear((1,3), 32)` for a shift set, or
    `Bank.boolean(3, 1, 32)` for the 3-input form.  Guessing which one a tuple
    means is exactly the kind of ambiguity that produces a silently different
    system, so the constructors are separate.

    Width is exact, not a truncation error: nothing above `width` can affect
    anything below it, ever.
    """

    def __init__(self, kind, driver, width=32, seed=1):
        if kind not in ("linear", "boolean"):
            raise ValueError("kind must be 'linear' or 'boolean'")
        self.kind, self.driver = kind, driver
        self.width = width
        self.mask = (1 << width) - 1
        self.x = seed & self.mask
        self.seed = seed
        self.t = 0

    @classmethod
    def linear(cls, S, width=32, seed=1):
        """x' = x XOR sum_{s in S} (x << s)."""
        return cls("linear", tuple(S), width, seed)

    @classmethod
    def boolean(cls, L, N, width=32, seed=1):
        """x' = x XOR G(a,b,c) with G given by its (L, N) parts."""
        return cls("boolean", (L, N), width, seed)

    def copy(self):
        return Bank(self.kind, self.driver, self.width, self.seed)

    def step(self):
        x, m = self.x, self.mask
        if self.kind == "linear":
            g = 0
            for s in self.driver:
                g ^= (x << s) & m
        else:
            L, N = self.driver
            a, b, c = (x << 1) & m, (x << 2) & m, (x << 3) & m
            g = 0
            if L & 1: g ^= a
            if L & 2: g ^= b
            if L & 4: g ^= c
            if N & 1: g ^= a & b
            if N & 2: g ^= a & c
            if N & 4: g ^= b & c
            if N & 8: g ^= a & b & c
        self.x = (x ^ g) & m
        self.t += 1
        return self.x

    def run(self, steps):
        return [self.step() for _ in range(steps)]

    def tapped(self, j, steps):
        """The signal on coordinate j, as a list of bits, from t = 0."""
        out = []
        for _ in range(steps):
            out.append((self.x >> j) & 1)
            self.step()
        return out

    def verify_period(self, j, k, cycles=4):
        """Check coordinate j has period 2^k exactly, and that 2^k is minimal."""
        n = (1 << k) * cycles
        sig = self.copy().tapped(j, n)
        P = 1 << k
        exact = all(sig[i] == sig[i + P] for i in range(n - P))
        smaller = any(all(sig[i] == sig[i + (1 << q)] for i in range(n - (1 << q)))
                      for q in range(k))
        return {"period_holds": exact, "minimal": exact and not smaller}


# --------------------------------------------------------------------- demo
def demo():
    print("=" * 70)
    print("1. a clock bank -- one shift-xor update, many known periods")
    print("=" * 70)
    S = (1, 3)
    cb = clock_bank(S, 20)
    print(f"driver x' = x XOR (x<<1) XOR (x<<3),  20 coordinates\n")
    print(f"{'period':>10}   coordinates sharing it")
    for k, js in cb["families"].items():
        print(f"{'2^'+str(k):>10}   {js}")
    print("\ncoordinates in the same row are distinct waveforms of the same")
    print("length, off the same register, for no extra work.")

    print("\n" + "=" * 70)
    print("2. the rendezvous -- when the bank re-synchronises")
    print("=" * 70)
    for name, (L, N) in (("a OR b   (rule 30's edge)", (3, 1)),
                         ("a        (pure shift)", (1, 0)),
                         ("a + b + c", (7, 0))):
        r = verify_rendezvous(L, N)
        prof = in_phase_profile(L, N)[:9]
        print(f"{name:<26} holds={str(r['holds']):>5}  "
              f"({r['tested']} live coordinates, {r['failures']} failures)")
        print(f"{'':<26} in phase at rung v: {prof}")
    print("\nwhere it holds, every signal of period dividing 2^v reads 1 at")
    print("t = -1 mod 2^v.  The count in phase is a function of v2(t+1) alone.")

    print("\n" + "=" * 70)
    print("3. synthesis -- pick a driver from a spec, with nothing run")
    print("=" * 70)
    spec = {20: 5, 8: 4}
    hits = synthesize(spec, smax=5)
    print("want: coordinate 20 of period 2^5 AND coordinate 8 of period 2^4")
    print("search space: 31 shift sets over shifts 1..5, none of them run\n")
    print(f"drivers meeting it: {len(hits)} of 31")
    best = cheapest(hits)
    print(f"cheapest member:    {best}  ({cost(best)} XOR per generation)")
    print(f"most expensive:     {max(hits, key=cost)}  "
          f"({cost(max(hits, key=cost))} XOR per generation)")
    b = Bank.linear(best, 32)
    print(f"verify {best}: j=20 -> {b.verify_period(20, 5)}")

    classes = spectrum_classes(5)
    print(f"\nthe spec picks out a class, never a single driver, because the")
    print(f"spectrum is a coarse invariant: 31 shift sets -> "
          f"{len(classes)} distinct spectra.")
    print(f"over shifts 1..7 it is 127 sets -> {len(spectrum_classes(7))} spectra.")
    print("so: meet the spec, then take the member with the fewest terms.")

    print("\n" + "=" * 70)
    print("4. structured sparsity -- freeze a sublattice by construction")
    print("=" * 70)
    for L, N, nm in ((2, 4, "b AND NOT c"), (4, 2, "c AND NOT a"),
                     (1, 1, "a AND NOT b")):
        ll = live_lattice(L, N)
        pr = predictable(L, N)
        print(f"{nm:>14}  stride={ll['stride']}  live={ll['live'][:8]}...  "
              f"{pr['reason']}")
    print("\nthe nonlinear term is a filter on which coordinates participate,")
    print("not an accelerant -- the survivors run the linear part's spectrum.")

    print("\n" + "=" * 70)
    print("5. the predictability dial")
    print("=" * 70)
    for L, N, nm in ((1, 0, "a"), (1, 1, "a AND NOT b"), (3, 1, "a OR b"),
                     (0, 3, "ab + ac")):
        pr = predictable(L, N)
        print(f"{nm:>14}  {'closed form' if pr['tame'] else 'no closed form':>15}"
              f"   ({pr['reason']})")
    print("\ntame when you want the period profile known at design time.")
    print("not tame when you want a spectrum nobody can write down.")


if __name__ == "__main__":
    demo()


# ---------------------------------------------------------------- dataframes
def clock_frame(S=(1, 3), width=24):
    """The clock bank as a table: coordinate, period, and its family size."""
    from sitoolz.dfkit import frame
    cb = clock_bank(S, width)
    sizes = {k: len(js) for k, js in cb["families"].items()}
    rows = []
    for j, k in enumerate(cb["p"]):
        if k is None:
            continue
        rows.append({"j": j, "p": k, "period": 1 << k,
                     "family_size": sizes[k],
                     "family": ",".join(map(str, cb["families"][k]))})
    return frame(rows, index="j")


def rendezvous_frame(width=24, steps=1 << 13, kmax=12, seed=1):
    """Which drivers re-synchronise -- one row per driver.

    Columns: name, L, N, live, tested, failures, holds, ones (G(1,1,1)),
    subcube, tame.  All 128 drivers, about a second.

    Note the two properties are independent: `a+b` is tame and does not
    rendezvous; plenty of wild drivers do.
    """
    from sitoolz.dfkit import frame
    rows = []
    for L in range(8):
        for N in range(16):
            r = verify_rendezvous(L, N, seed, width, steps, kmax)
            tt = dl.truth(L, N)
            pr = predictable(L, N)
            rows.append({"name": dl.name(L, N), "L": L, "N": N,
                         "live": r["live"], "tested": r["tested"],
                         "failures": r["failures"], "holds": r["holds"],
                         "ones": bool(dl.fires_on_ones(tt)),
                         "subcube": dl.is_subcube(tt), "tame": pr["tame"]})
    return frame(rows)


def class_frame(smax=5, jmax=63):
    """The spectrum equivalence classes -- what the design space collapses to.

    One row per distinct spectrum: how many shift sets share it, the cheapest
    member, its cost, and the gcd/min the class sits at.
    """
    from sitoolz.dfkit import frame
    rows = []
    for i, sets in enumerate(spectrum_classes(smax, jmax)):
        best = cheapest(sets)
        g = 0
        for s in best:
            a, b = g, s
            while b:
                a, b = b, a % b
            g = a
        rows.append({"class": i, "n_sets": len(sets), "gcd": g,
                     "min": min(best), "cheapest": "(" + ",".join(map(str, best)) + ")",
                     "cost": cost(best),
                     "dearest_cost": max(len(S) for S in sets),
                     "members": " ".join("(" + ",".join(map(str, S)) + ")"
                                         for S in sets[:6])
                                + (" …" if len(sets) > 6 else "")})
    return frame(sorted(rows, key=lambda r: -r["n_sets"]))
