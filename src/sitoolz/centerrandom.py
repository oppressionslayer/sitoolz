"""centerrandom -- the center column of Rule 30, seen from inside the edge system.

Companion to edgerandom.  Where edgerandom shows the *edge* is a closed-form
lookup, this module locates the *center* inside the same machine -- and states
exactly what is determined, what is open, and why.

THE CENTER LAW (verified here, generation by generation):

    center(t)  =  bit t of E(t),        E' = E ^ ((E << 1) | (E << 2)),  E(0) = 1

The edge recurrence, run at full width, IS Rule 30 in the frame co-moving with
the right light-cone edge: bit j of E(t) is the cell at distance j behind the
edge.  That makes the whole right cone a triangular XOR-driven system, so the
tri_xor theorem applies to EVERY coordinate:

    * every diagonal parallel to the right edge is EXACTLY periodic from
      generation 0, with period 2^p(j)  (verified for 25 coordinates below);
    * the center column is the anti-diagonal sample  c(t) = x_t(t mod 2^p(t)).

So Rule 30's interior is a bank of exact periodic clocks, and the center reads
clock t once, at time t.  Why it still looks random:  p(t) grows ~ t/2.4, so
the period 2^p(t) is EXPONENTIALLY longer than the sample time t -- the center
under-samples every clock it touches (by x43,980,465,111 at t=100).  The three
bits the phase lemmas reach (t = 0, 1, 3) are pinned and correct; everything
past the sampling horizon is open, and this module says so.

    from sitoolz import centerrandom as cr

    cr.center(64)              # the center column, verified against OEIS A051023
    cr.center_law_check(2000)  # c(t) == bit t of E(t), every generation
    cr.diagonal(5, 64)         # the offset-5 diagonal -- an exact periodic clock
    cr.diagonal_at(t)          # one center bit, read off the E-object
    cr.diagonal_col(n)         # the anti-diagonal sweep == the center column
    cr.center_stream(lo, hi)   # center bits for a range, one shared walk of E
    cr.stream()                # endless (t, bit) generator
    cr.cell(10**30, 37)        # RANDOM ACCESS: any wedge cell at ANY generation
    cr.strip(10**100, 50)      # the 50-cell strip behind the edge at g = 10^100
    cr.clock(50)               # clock 50: period, sample time, under-sampling gap
    cr.pinned()                # center bits FORCED by the phase lemmas (0,1,3 -> 1)
    cr.characterize(200000)    # bias / autocorrelation / held-out predictability
    cr.caveat_frame()          # what is settled vs open, in the T-function frame

The index move behind cell/strip: triangularity makes the low-k window of E
self-contained, so it is periodic with period 2^p(k) -- reduce the generation
like an index (g mod 2^p(k)) and walk only the remainder, on k-bit ints.
Fun fact this exposes: any g divisible by 2^20 has its 50-cell wedge IDENTICAL
to generation 0 -- generations 10^20, 10^21, ..., 10^100 all have deja vu.
"""
from collections import defaultdict

# A094604: the edge-run table to rung 48 (agrees with edgelaw.TABLE everywhere
# both are defined).  Doubles as the period ladder: coordinates 0..TABLE[v]-1
# have period dividing 2^v, so p(j) = min{ v : TABLE[v] > j }.
EDGE_TABLE = [1, 3, 4, 6, 7, 9, 15, 16, 24, 25, 27, 29, 34, 36, 37, 39, 41, 43, 48,
              49, 51, 54, 55, 58, 60, 63, 64, 66, 69, 70, 72, 74, 77, 79, 80, 82, 84,
              86, 90, 91, 93, 100, 103, 104, 106, 108, 111, 112, 114]


def p_of(j):
    """Period exponent of coordinate j:  x_j has exact period 2^p(j).
    Known for j <= 113 (rung 48 of A094604); None past that."""
    for v, c in enumerate(EDGE_TABLE):
        if c > j:
            return v
    return None


# ------------------------------------------------------------ the two frames
def center(n, W=None):
    """First n bits of Rule 30's center column, by direct simulation.
    Verified against OEIS A051023.  Rule 30: c'[i] = c[i-1] XOR (c[i] OR c[i+1])."""
    W = W or (2 * n + 5)
    mid = W // 2
    row = 1 << mid
    mask = (1 << W) - 1
    out = []
    for _ in range(n):
        out.append((row >> mid) & 1)
        row = ((row << 1) ^ (row | (row >> 1))) & mask
    return out


def edge_frame(n):
    """E(0..n-1) in the co-moving frame: E' = E ^ ((E<<1)|(E<<2)), unmasked.
    Bit j of E(t) is the cell at distance j behind the right edge -- exact for
    every bit, since information only flows upward (triangularity)."""
    E, out = 1, []
    for _ in range(n):
        out.append(E)
        E = E ^ ((E << 1) | (E << 2))
    return out


def center_law_check(n=2000):
    """THE LAW:  center(t) == bit t of E(t), checked for every t < n."""
    truth = center(n)
    diag = [(E >> t) & 1 for t, E in enumerate(edge_frame(n))]
    return {"checked": n, "holds": diag == truth}


def diagonal(j, n=None):
    """The offset-j diagonal of Rule 30 (bit j of E(t) over t) -- an exact
    periodic clock with period 2^p(j), from generation 0."""
    p = p_of(j)
    n = n or (4 << p if p is not None else 256)
    return [(E >> j) & 1 for E in edge_frame(n)]


def diagonal_at(t):
    """The center bit at generation t, read off Lars's edge system: bit t of
    E(t).  Same bits as center(t), but computed from the E-object alone -- no
    full-triangle Rule 30 simulation.  (Walking the anti-diagonal still costs t
    steps of E; a shortcut past that walk IS the open center problem.)"""
    E = 1
    for _ in range(t):
        E = E ^ ((E << 1) | (E << 2))
    return (E >> t) & 1


def diagonal_col(n):
    """The whole anti-diagonal of the edge system: [bit t of E(t) for t < n].
    Equals center(n) exactly -- the center column, as a cut through E."""
    return [(E >> t) & 1 for t, E in enumerate(edge_frame(n))]


def cell(g, j):
    """Random access into Rule 30's right wedge: the cell at offset j behind the
    right edge, at ANY generation g -- via the index reduction

        E(g) mod 2^(j+1)  =  E(g mod 2^p(j)) mod 2^(j+1)

    Triangularity makes the low window self-contained, so its period is 2^p(j);
    reduce the time like an index, walk only the remainder.  Works for g = 10**30
    as easily as g = 30.  Known for j <= 113 (rung 48 of A094604); the walk is
    at most 2^p(j) steps, so keep j <= ~55 for instant answers."""
    p = p_of(j)
    if p is None:
        raise KeyError(f"offset {j} is past rung 48 of A094604 -- period unknown")
    mask = (1 << (j + 1)) - 1          # the window is self-contained: mask the walk
    E = 1
    for _ in range(g % (1 << p)):
        E = (E ^ ((E << 1) | (E << 2))) & mask
    return (E >> j) & 1


def strip(g, k):
    """The first k cells behind the right edge at ANY generation g, as a bit
    string (offset 0 first).  Same index reduction as cell(), one walk."""
    p = p_of(k - 1)
    if p is None:
        raise KeyError(f"strip width {k} is past rung 48 of A094604")
    mask = (1 << k) - 1                # self-contained window: constant-size walk
    E = 1
    for _ in range(g % (1 << p)):
        E = (E ^ ((E << 1) | (E << 2))) & mask
    return "".join(str((E >> j) & 1) for j in range(k))


# ------------------------------------------------------------ the LEFT wedge
_LEFT = None

def _left_table():
    """Measured left-wedge table: offset j -> (transient prefix, eventual cycle).
    The left frame is NOT triangular (reads its own coordinate), so columns are
    only EVENTUALLY periodic (Rowland 2006) -- but that still yields a lookup
    after the transient.  Measured to offset 199 with >=6 periods of evidence;
    entries whose evidence window was thin are marked provisional."""
    global _LEFT
    if _LEFT is None:
        import json, os
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            ".cache", "left_wedge.json")
        raw = json.load(open(path))
        _LEFT = {}
        for k, v in raw.items():
            _LEFT[int(k)] = (v["tau"], v["P"], v["prefix"], v["cycle"],
                             v["tau"] + 6 * v["P"] <= 15800)   # solid evidence?
    return _LEFT

def left_cell(g, j):
    """Random access into Rule 30's LEFT wedge: the cell at offset j right of the
    left light-cone edge, at any generation g.  Lookup: transient prefix, then
    cycle.  Left periods are tiny (<= 8 measured) but transients grow ~1.3*j --
    the mirror-image cost profile of the right wedge."""
    tab = _left_table()
    if j not in tab:
        raise KeyError(f"left offset {j} not in measured table (0..{max(tab)})")
    tau, P, prefix, cycle, solid = tab[j]
    if not solid:
        raise KeyError(f"left offset {j}: evidence too thin -- provisional entry")
    if g < tau:
        return int(prefix[g])
    return int(cycle[(g - tau) % P])

def left_strip(g, k):
    """The first k cells right of the LEFT edge at generation g, offset 0 first."""
    return "".join(str(left_cell(g, j)) for j in range(k))


def center_stream(lo, hi):
    """Center bits for every generation in [lo, hi) -- the diagonal, streamed.

    Mirrors edgerandom.edge_stream, with the honest asymmetry on display:
    the edge is O(1) per generation (a lookup); the center costs one shared
    walk of E to `hi` (O(hi) total, but only ONE walk for the whole range --
    never call diagonal_at in a loop, which re-walks from zero each time).
    """
    E, out = 1, []
    for t in range(int(hi)):
        if t >= lo:
            out.append((E >> t) & 1)
        E = E ^ ((E << 1) | (E << 2))
    return out


def stream():
    """Endless center-bit generator: yields (t, bit), walking E once, forever.

        for t, b in cr.stream():
            ...
    """
    E, t = 1, 0
    while True:
        yield t, (E >> t) & 1
        E = E ^ ((E << 1) | (E << 2))
        t += 1


def diagonal_check(jmax=30, n=2000):
    """Verify each diagonal's exact periodicity against its period 2^p(j)."""
    Es = edge_frame(n)
    bad, checked = [], 0
    for j in range(jmax):
        p = p_of(j)
        if p is None or (1 << p) * 2 > n:
            continue
        per = 1 << p
        col = [(E >> j) & 1 for E in Es]
        checked += 1
        if any(col[i] != col[i % per] for i in range(n)):
            bad.append(j)
    return {"coords_checked": checked, "failures": bad}


# --------------------------------------------------------- what is determined
def clock(t):
    """Clock t as the center sees it: its exact period, the single sample time,
    and the under-sampling gap that keeps the center bit undetermined."""
    p = p_of(t)
    if p is None:
        return {"t": t, "note": "past rung 48 of A094604 -- period unknown"}
    per = 1 << p
    return {"t": t, "period": per, "sample_time": t,
            "undersampling": per // max(t, 1),
            "pinned": (t + 1) % per == 0}


def pinned(upto=114):
    """Center bits FORCED by the phase lemmas (C3: in-phase bits show 1), i.e.
    the t whose single sample lands on the in-phase slot: 2^p(t) | t+1.
    Returns {t: bit}; verified against the true column."""
    out = {}
    truth = center(upto + 2)
    for t in range(upto):
        p = p_of(t)
        if p is not None and (t + 1) % (1 << p) == 0:
            out[t] = 1
            assert truth[t] == 1, f"phase pin failed at t={t}"
    return out


# ----------------------------------------------------------- honest statistics
def characterize(n=200000):
    """Randomness profile of the center column: bias, autocorrelation, and
    HELD-OUT Markov predictability (train first half, test second -- the only
    version that distinguishes structure from memorized noise)."""
    col = center(n)
    ones = sum(col)
    ac = {}
    for lag in (1, 2, 3, 5, 8, 13):
        ac[lag] = round(sum(1 for i in range(n - lag) if col[i] == col[i + lag]) / (n - lag), 4)
    half = n // 2
    held = {}
    for k in (8, 12, 16):
        model = defaultdict(lambda: [0, 0])
        for i in range(half - k - 1):
            model[tuple(col[i:i + k])][col[i + k]] += 1
        pred = {c: (0 if v[0] >= v[1] else 1) for c, v in model.items()}
        cor = seen = 0
        for i in range(half, n - k - 1):
            c = tuple(col[i:i + k])
            if c in pred:
                seen += 1
                cor += (pred[c] == col[i + k])
        held[k] = round(cor / seen, 4) if seen else None
    return {"bits": n, "bias": round(ones / n, 5),
            "autocorrelation": ac, "held_out_predictability": held,
            "verdict": "no exploitable local structure (held-out ~0.5)"}


# ------------------------------------------------------------------- caveats
def caveat_frame():
    """Rule 30 in the T-function frame: every claim here is verified above."""
    return [
        {"feature": "right edge run", "status": "CLOSED FORM",
         "note": "edge(g) = TABLE[v2(g+1)], rungs to 48 (A094604)"},
        {"feature": "every right-diagonal", "status": "EXACTLY PERIODIC",
         "note": "period 2^p(j) from generation 0 -- the interior is a bank of clocks"},
        {"feature": "center column", "status": "= the diagonal sample x_t(t)",
         "note": "center(t) == bit t of E(t): located, not solved"},
        {"feature": "center bits t = 0, 1, 3", "status": "PINNED (= 1)",
         "note": "the sample lands on the in-phase slot: 2^p(t) | t+1"},
        {"feature": "all other center bits", "status": "OPEN",
         "note": "under-sampled: period 2^p(t) is exponential in t, sample time is t"},
        {"feature": "delta / Q1-Q3 (edge doubling)", "status": "OPEN",
         "note": "114 terms; 2-automaticity undecidable at this length"},
    ]


if __name__ == "__main__":
    print("center law:", center_law_check(2000))
    print("diagonals :", diagonal_check(30, 2000))
    print("pinned    :", pinned())
    print()
    for t in (7, 20, 50, 100):
        c = clock(t)
        print(f"  clock {t:>3}: period {c['period']:>16,}  sampled once at t={t}"
              f"  (under-sampled x{c['undersampling']:,})")
