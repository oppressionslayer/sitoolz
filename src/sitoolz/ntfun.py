"""ntfun -- poking at rule 30 for fun.

Frames for the structure that turns up when you look at rule 30 rows as
numbers rather than pictures.

    from sitoolz import ntfun

    ntfun.spikes(2040, 2060)      # zero runs, and what happens at 2^k
    ntfun.where(2048)             # where the long run actually sits
    ntfun.parity(1, 20)           # mod 4, the k that satisfy v^k == v+k
    ntfun.periods(16, 2000)       # low bits as functions of generation
    ntfun.gaps(1, 500)            # gap distribution between set bits
    ntfun.squares(1, 3000)        # values that are perfect squares
    ntfun.centre_index(4096)      # centre column vs the square machines
"""
from math import isqrt
from collections import Counter

from sitoolz import ntmc

try:
    import pandas as pd
except ImportError:
    pd = None


try:
    import gmpy2 as _gmpy2
except ImportError:
    _gmpy2 = None

#: gmpy2 present?  When it is, the tower functions use GMP for isqrt and
#: multiplication -- 23x on a 2,000,002-bit isqrt, 21x on the squaring, so
#: 14x on `deep_pack` and 5.5x on `deep_unpack`.  Results are identical and
#: everything still hands back plain python ints; the conversion at the
#: boundary costs 0.2 ms.  Pass fast=False to force CPython, fast=True to
#: insist on gmpy2 and fail loudly if it is missing.
HAVE_GMPY2 = _gmpy2 is not None


def _need_pd():
    if pd is None:
        raise ImportError("pandas required")


def _fast(fast):
    """Resolve the fast= argument to (mpz-or-int, isqrt).  None means auto."""
    if fast is None:
        fast = HAVE_GMPY2
    if fast:
        if _gmpy2 is None:
            raise ImportError("fast=True needs gmpy2; pip install gmpy2")
        return _gmpy2.mpz, _gmpy2.isqrt
    return int, isqrt


def _val(g, w=None):
    return ntmc.flipstrip(ntmc.row(g, w))


def _zero_runs(v):
    """(length, start_index_from_top) of every zero run in the value."""
    s = bin(v)[2:]
    out, i = [], 0
    while i < len(s):
        if s[i] == "0":
            j = i
            while j < len(s) and s[j] == "0":
                j += 1
            out.append((j - i, i))
            i = j
        else:
            i += 1
    return out


# =====================================================================
# the 2^k spike
# =====================================================================
def spikes(lo, hi, w=None):
    """Longest zero run per generation.  Powers of two are flagged."""
    _need_pd()
    w = w or ntmc.width_for(hi)
    rows = []
    for g in range(lo, hi + 1):
        v = _val(g, w)
        rr = _zero_runs(v)
        best = max(rr, default=(0, 0))
        rows.append({
            "gen": g,
            "longest_zero_run": best[0],
            "at_index": best[1],
            "span": v.bit_length(),
            "pow2": (g & (g - 1)) == 0 and g > 0,
            "popcount": bin(v).count("1"),
        })
    return pd.DataFrame(rows).set_index("gen")


def where(g, w=None, context=8):
    """Locate the longest zero run inside the row and report it several ways."""
    _need_pd()
    w = w or ntmc.width_for(g)
    v = _val(g, w)
    s = bin(v)[2:]
    n = len(s)
    L, i = max(_zero_runs(v), default=(0, 0))
    mid = i + L / 2
    return pd.DataFrame([{
        "gen": g,
        "run_len": L,
        "span": n,
        "from_top": i,
        "from_bottom": n - (i + L),
        "centre_of_run": mid,
        "row_centre": (n - 1) / 2,
        "offset_from_centre": mid - (n - 1) / 2,
        "frac_from_top": round(i / n, 4),
        "context": s[max(0, i - context): i + L + context],
    }]).set_index("gen")


def decay(g, steps=12, w=None):
    """The 2^k spike eroding: run length per generation after g."""
    _need_pd()
    w = w or ntmc.width_for(g + steps)
    rows = []
    prev = None
    for k in range(steps + 1):
        v = _val(g + k, w)
        L = max((a for a, _ in _zero_runs(v)), default=0)
        rows.append({"gen": g + k, "run": L,
                     "delta": None if prev is None else L - prev})
        prev = L
    return pd.DataFrame(rows).set_index("gen")


# =====================================================================
# arithmetic structure
# =====================================================================
def parity(lo, hi, w=None, kmax=16):
    """mod 4, the k with v^k == v+k, and whether the value can be a square."""
    _need_pd()
    w = w or ntmc.width_for(hi)
    rows = []
    for g in range(lo, hi + 1):
        v = _val(g, w)
        ks = [k for k in range(1, kmax) if (v ^ k) == (v + k)]
        r = isqrt(v)
        rows.append({
            "gen": g,
            "value": v,
            "mod4": v % 4,
            "bit1": (v >> 1) & 1,
            "bit2": (v >> 2) & 1,
            "xor_add_k": ks,
            "square_possible": v % 4 in (0, 1),
            "is_square": r * r == v,
        })
    return pd.DataFrame(rows).set_index("gen")


def periods(nbits=16, gens=2000, w=None, pmax=512):
    """Smallest p (power of 2) for which bit b is a function of gen mod p."""
    _need_pd()
    w = w or ntmc.width_for(gens)
    vals = [(g, _val(g, w)) for g in range(1, gens + 1)]
    rows = []
    for b in range(nbits):
        found, pat = None, None
        p = 1
        while p <= pmax:
            cls = [{(v >> b) & 1 for g, v in vals if g % p == r} for r in range(p)]
            if all(len(c) <= 1 for c in cls):
                found = p
                pat = "".join(str(sorted(c)[0]) if c else "-" for c in cls)
                break
            p *= 2
        rows.append({"bit": b, "period": found,
                     "pattern": pat if pat and len(pat) <= 40 else
                                (pat[:37] + "..." if pat else None)})
    return pd.DataFrame(rows).set_index("bit")


def squares(lo, hi, w=None):
    """Values that are perfect squares, and how many were ruled out for free."""
    _need_pd()
    w = w or ntmc.width_for(hi)
    hits, ruled = [], 0
    for g in range(lo, hi + 1):
        v = _val(g, w)
        if v % 4 == 3:
            ruled += 1
            continue
        r = isqrt(v)
        if r * r == v:
            hits.append({"gen": g, "value": v, "root": r})
    df = pd.DataFrame(hits)
    df.attrs["ruled_out_by_mod4"] = ruled
    df.attrs["checked"] = hi - lo + 1
    return df


# =====================================================================
# gaps
# =====================================================================
def gaps(lo, hi, w=None):
    """Distribution of gaps between consecutive set bits, vs geometric p=1/2."""
    _need_pd()
    w = w or ntmc.width_for(hi)
    h = Counter()
    for g in range(lo, hi + 1):
        v = _val(g, w)
        b = [i for i in range(v.bit_length() - 1, -1, -1) if (v >> i) & 1]
        h.update(b[i] - b[i + 1] for i in range(len(b) - 1))
    tot = sum(h.values())
    rows = []
    for k in sorted(h):
        rows.append({
            "gap": k,
            "count": h[k],
            "frac": h[k] / tot,
            "geometric": 2.0 ** -k,
            "ratio_to_prev": h[k] / h[k - 1] if k - 1 in h else None,
        })
    return pd.DataFrame(rows).set_index("gap")


# =====================================================================
# the edge law
# =====================================================================
def _v2(n):
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def edge_evolve(gens, W=256):
    """Evolve ONLY the light-cone edge.  O(1) per generation, not O(row).

    Let E[j] be the cell j below the leading edge.  The cone grows one cell
    per generation and nothing above it is live, so

        E'[j] = E[j] XOR ( E[j-1] OR E[j-2] )

    which as whole-word operations is just

        E' = E ^ ((E << 1) | (E << 2))

    The edge never consults the interior, so a narrow window suffices.
    Verified against full rows for every generation to 600.  100,000
    generations take ~0.01 s, against ~3 s for full evolution.
    """
    m = (1 << W) - 1
    E = 1                                    # generation 0: the seed
    out = [E]
    for _ in range(gens):
        E = (E ^ ((E << 1) | (E << 2))) & m
        out.append(E)
    return out


def edge_run_at(gens, W=256):
    """Length of the run of ones at the edge, per generation, via edge_evolve."""
    out = []
    for E in edge_evolve(gens, W):
        k = 0
        while (E >> k) & 1:
            k += 1
        out.append(k)
    return out


def edge(lo=1, hi=4095, w=None):
    """Per generation: the run of ones at the light-cone edge, and v2(g+1).

    The edge run is a pure function of v2(g+1) -- the odd part of g+1
    does not matter.  Verified exact for every generation to 4095.
    """
    _need_pd()
    w = w or ntmc.width_for(hi)
    rows = []
    for g in range(lo, hi + 1):
        s = bin(_val(g, w))[2:]
        rows.append({
            "gen": g,
            "v2_of_g_plus_1": _v2(g + 1),
            "edge_ones": len(s) - len(s.rstrip("1")),
            "longest_ones": max(len(x) for x in s.split("0")),
            "longest_zeros": max((len(x) for x in s.split("1")), default=0),
        })
    return pd.DataFrame(rows).set_index("gen")


def edge_law(lo=1, hi=4095, w=None):
    """The law itself: v2(g+1) -> edge run, with a check that it is single-valued."""
    _need_pd()
    df = edge(lo, hi, w)
    out = []
    for v, grp in df.groupby("v2_of_g_plus_1"):
        vals = sorted(grp.edge_ones.unique())
        out.append({
            "v2": v,
            "edge_ones": vals[0] if len(vals) == 1 else vals,
            "single_valued": len(vals) == 1,
            "n_gens": len(grp),
            "example_gens": sorted(grp.index)[:4],
        })
    return pd.DataFrame(out).set_index("v2")


# =====================================================================
# the edge index -- look it up instead of running it
# =====================================================================
EDGE_TABLE = {
    0: 1,   1: 3,   2: 4,   3: 6,   4: 7,   5: 9,   6: 15,  7: 16,
    8: 24,  9: 25,  10: 27, 11: 29, 12: 34, 13: 36, 14: 37, 15: 39,
    16: 41, 17: 43, 18: 48, 19: 49, 20: 51, 21: 54, 22: 55, 23: 58,
    24: 60,
}


class EdgeIndex:
    """The light-cone edge run, by lookup rather than simulation.

    The run of ones at rule 30's leading edge depends only on the 2-adic
    valuation of g+1 -- the odd part of g+1 is irrelevant.  Verified exact
    for every generation to 16,777,215.

        ix = ntfun.EdgeIndex()
        ix.run(16777215)      # 60, with no simulation at all
        ix.rank(1791)         # 8   -- 1792 = 7 * 2^8
        ix.unrank(8)          # 255 -- the first generation at that level
        ix.frame(1, 40)       # rows, values, and the lookup side by side
    """

    def __init__(self, table=None):
        self.table = dict(table or EDGE_TABLE)
        self.vmax = max(self.table)

    # ---- rank / unrank -------------------------------------------------
    def rank(self, g):
        """v2(g+1) -- which level of the law this generation sits on."""
        n = g + 1
        return (n & -n).bit_length() - 1

    def unrank(self, v):
        """The first generation at level v.  Others are (2m+1)*2^v - 1."""
        return (1 << v) - 1

    def siblings(self, v, count=6):
        """Generations sharing level v -- all give the same edge run."""
        return [(2 * m + 1) * (1 << v) - 1 for m in range(count)]

    # ---- the lookup ----------------------------------------------------
    def run(self, g):
        """Edge run at generation g.  Pure lookup, no simulation."""
        v = self.rank(g)
        if v > self.vmax:
            raise KeyError(f"v2(g+1) = {v} exceeds the verified table (<= {self.vmax})")
        return self.table[v]

    def bits(self, g):
        """The edge bits themselves.  Always a solid run of ones, so the
        length from run() is the whole story."""
        return "1" * self.run(g)

    @staticmethod
    def of_value(value):
        """Pull the edge off a value you already have."""
        s = bin(value)[2:]
        return s[len(s.rstrip("1")):]

    # ---- checking it ---------------------------------------------------
    def measure(self, g, W=512):
        """The edge run actually computed, for comparison with run()."""
        return edge_run_at(g, W)[g]

    def frame(self, lo=1, hi=40, rows=True):
        """gen, v2, looked-up run, measured run, the row, and flipstrip.

        `rows=False` drops the row and value columns, which get large fast.
        """
        _need_pd()
        w = ntmc.width_for(hi) if rows else None
        meas = edge_run_at(hi, 512)
        out = []
        for g in range(lo, hi + 1):
            d = {
                "gen": g,
                "v2": self.rank(g),
                "lookup": self.run(g),
                "measured": meas[g],
                "agree": self.run(g) == meas[g],
            }
            if rows:
                r = ntmc.row(g, w)
                v = ntmc.flipstrip(r)
                s = bin(v)[2:]
                k = len(s) - len(s.rstrip("1"))
                d["row"] = r
                d["value"] = v
                d["value_bin"] = s
                d["span"] = 2 * g + 1
                d["edge_bits"] = s[-k:] if k else ""
                d["tail"] = ("..." + s[-(k + 8):-k] + "|" + s[-k:]) if k else s[-12:]
                c = w // 2
                d["centre"] = (r >> c) & 1
                d["pred_centre"] = (None if g == 0
                                    else ntmc.center_bit_sync(r, w))
            out.append(d)
        df = pd.DataFrame(out)
        big = [c for c in ("row", "value") if c in df]
        data = {c: pd.Series(df[c], dtype=object if c in big else None)
                for c in df.columns}
        return pd.DataFrame(data).set_index("gen")

    def verify(self, upto=4095, W=512):
        """Check the lookup against measurement for every generation."""
        meas = edge_run_at(upto, W)
        bad = [g for g in range(1, upto + 1)
               if self.rank(g) <= self.vmax and self.run(g) != meas[g]]
        return {"checked": upto, "mismatches": len(bad), "first_bad": bad[:5]}

    def __repr__(self):
        return (f"<EdgeIndex v2 0..{self.vmax}, "
                f"runs {[self.table[v] for v in range(min(6, self.vmax + 1))]}...>")


# =====================================================================
# forbidden bits: squares and rule 30 side by side
# =====================================================================
def _ks(B, kmax=64):
    """Every k whose set bits lie inside the forbidden set B.

    v ^ k == v + k exactly when k borrows only bits that v never sets,
    so the k-set is the whole sublattice on B: 2^|B| - 1 values.
    """
    return [k for k in range(1, kmax)
            if all(not (k >> b) & 1 or b in B for b in range(k.bit_length()))]


def forbidden(mod=8, gens=2000, nbits=10, kmax=64, w=None):
    """Which low bits a family never sets, and the k that follow.

    Squares appear as a row for comparison: they forbid exactly bit 1,
    which is why k=2 is their only constant.  Rule 30 values forbid more,
    and the count grows with how many 2s divide the generation.
    """
    _need_pd()
    from math import isqrt
    w = w or ntmc.width_for(gens)
    vals = {g: _val(g, w) for g in range(1, gens + 1)}

    rows = []

    sq = [m * m for m in range(1, 4000)]
    B = [b for b in range(nbits) if not any((s >> b) & 1 for s in sq)]
    rows.append({
        "family": "squares",
        "forbidden_bits": B,
        "n_forbidden": len(B),
        "k": _ks(B, kmax),
        "n_k": len(_ks(B, kmax)),
        "members": len(sq),
        "example": 25,
    })

    for r in range(mod):
        gs = [g for g in range(1, gens + 1) if g % mod == r]
        if not gs:
            continue
        B = [b for b in range(nbits)
             if not any((vals[g] >> b) & 1 for g in gs)]
        k = _ks(B, kmax)
        rows.append({
            "family": f"g mod {mod} = {r}",
            "forbidden_bits": B,
            "n_forbidden": len(B),
            "k": k,
            "n_k": len(k),
            "members": len(gs),
            "example": vals[gs[0]],
        })

    df = pd.DataFrame(rows).set_index("family")
    df["lattice_ok"] = df.n_k == (2 ** df.n_forbidden - 1).clip(upper=kmax - 1)
    return df


# =====================================================================
# how near the xor-k square machine does a rule 30 value fall
# =====================================================================
def leak(n):
    """Distance to the nearest perfect square, either side."""
    from math import isqrt
    r = isqrt(n)
    return min(n - r * r, (r + 1) ** 2 - n)


def xor_offsets(m=5, kmax=60000, pmax=256):
    """What `k*k ^ m` really is: `k*k` plus a fixed offset, by `k mod P`.

    The mask only touches low bits, and squares take few values mod 2^b, so
    the xor is not really bitwise at all -- it collapses to a handful of
    fixed +/- offsets on a power-of-two period.  For m=5 that is

        k = 0 (mod 4)   ->   k*k + 5
        k = 1, 3        ->   k*k + 3
        k = 2           ->   k*k - 3

    which is why the xor5 machine sits within 3 of the plain squares in
    `centre_index`, and -- the useful part -- why `isqrt(k*k ^ 5)` comes
    back as k-1 exactly when k = 2 (mod 4).  The offset goes negative there
    and drops the value below its own square.  The 1M centre column has
    big % 4 == 2, which is the whole story of that off-by-one.

    m=2 is the odd one out with period 1: bit 1 is zero in every square, so
    `k*k ^ 2` is `k*k + 2` for every k with no case split at all.  It is the
    only nonzero mask that manages that.

      offset      (k*k ^ m) - k*k for this residue
      keeps_index isqrt(k*k ^ m) is still k, i.e. the offset did not go
                  negative -- holds once k is past the first few

    Verified against every k < kmax.
    """
    _need_pd()
    tab = None
    for P in [1 << j for j in range(pmax.bit_length())]:
        seen, ok = {}, True
        for k in range(kmax):
            d = (k * k ^ m) - k * k
            r = k % P
            if seen.setdefault(r, d) != d:
                ok = False
                break
        if ok:
            tab = seen
            break
    if tab is None:
        raise ValueError(f"no period <= {pmax} for mask {m}")
    big_k = 1 << 20                                    # past the small-k edge
    rows = [{"k_mod_P": r, "offset": tab[r],
             "keeps_index": isqrt((big_k + r) ** 2 ^ m) == big_k + r}
            for r in sorted(tab)]
    df = pd.DataFrame(rows).set_index("k_mod_P")
    df.attrs["mask"] = m
    df.attrs["period"] = P
    df.attrs["checked_to"] = kmax
    return df


def square_family(lo=0, hi=36, m=5, sort=False):
    """The x that are square AND have isqrt(x) ^ m square -- in closed form.

    The brute-force reading is "test every x below the bound twice".  The
    closed form is one line: set x = n*n, then n ^ m must be square, so
    n = i*i ^ m and x = (i*i ^ m)**2.  One expression, no is_square calls,
    exact at any scale -- 1,643,000 candidates become 36 rows of arithmetic.

    Every column is computed, never searched:

      i          the index -- this is the whole address
      i_mod      i % period, which picks the offset (see `xor_offsets`)
      offset     what ^ m really does to i*i:  +5, +3, -3 for m=5
      root       i*i ^ m, i.e. isqrt(x)
      x          root**2 -- the member
      inner      isqrt(x) ^ m, which is i*i again, closing the loop
      i_square   i is itself a square, so x is a double square
      gap        x minus the previous x, in i order
      digits     decimal length of x

    `sort=True` gives ascending x instead of ascending i.  They differ at
    the start -- i = 0, 1, 2 give x = 25, 16, 1 -- because the xor can flip
    a high bit off while i is still small.

    Density is a fourth root: |{x < X}| ~ X**(1/4), which is exactly what
    `numtoolsd.count_below(X, m)` returns without scanning anything.
    """
    _need_pd()
    off = xor_offsets(m)
    P = off.attrs["period"]
    rows, prev = [], None
    for i in range(lo, hi):
        n = i * i
        root = n ^ m
        x = root * root
        rows.append({
            "i": i,
            "i_mod": i % P,
            "offset": root - n,
            "root": root,
            "x": x,
            "inner": isqrt(x) ^ m,
            "i_square": isqrt(i) ** 2 == i,
            "gap": None if prev is None else x - prev,
            "digits": len(str(x)),
        })
        prev = x
    df = pd.DataFrame(rows).set_index("i")
    if sort:
        df = df.sort_values("x")
    df.attrs["m"] = m
    df.attrs["period"] = P
    df.attrs["closed_form"] = "x = (i*i ^ m)**2"
    df.attrs["exact"] = bool((df.inner.values == df.index.values ** 2).all())
    return df


def square_reach(lo=0, hi=36, m=5):
    """`square_family` extended rightward -- each member to the next one.

    The family is sparse and gets sparser: x = (i*i ^ m)**2 is about i**4,
    so consecutive members pull apart like the derivative, 4*i**3.  These
    columns say how far you would have to walk from one member to the next
    if you were still brute-forcing, and what fraction of that walk the
    closed form skips.

      x           the member at i
      next_x      the member at i+1
      gap         next_x - x
      skipped     integers strictly between -- all of them non-members
      four_i3     4*i**3, the asymptotic gap from d(i**4)/di
      gap_ratio   gap / four_i3, which walks to 1 as i grows
      density     1 / gap, members per unit at this scale
      seen        how many members up to and including this one

    The last row has no successor, so its forward columns are None.
    """
    _need_pd()
    rows = []
    for i in range(lo, hi):
        x = (i * i ^ m) ** 2
        nxt = ((i + 1) * (i + 1) ^ m) ** 2 if i + 1 < hi else None
        gap = None if nxt is None else nxt - x
        f = 4 * i ** 3
        rows.append({
            "i": i,
            "root": i * i ^ m,
            "x": x,
            "next_x": nxt,
            "gap": gap,
            "skipped": None if gap is None else abs(gap) - 1,
            "four_i3": f,
            "gap_ratio": None if (gap is None or f == 0) else round(gap / f, 4),
            "density": None if not gap else round(1 / abs(gap), 12),
            "seen": i - lo + 1,
        })
    big = ("x", "next_x", "gap", "skipped", "four_i3")
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("i")
    df.attrs["m"] = m
    df.attrs["total_skipped"] = sum(r["skipped"] or 0 for r in rows)
    df.attrs["span"] = rows[-1]["x"] - rows[0]["x"]
    return df


def square_window(lo=4, hi=16, m=5, reach=3, show="gap", digits=12):
    """The family read sideways: each index with its neighbours either way.

    `square_reach` only looks forward.  This puts the index in the middle
    and opens columns left and right of it, so a row is a window on the
    family rather than a single member.  `reach=32` gives the 65 columns
    you would want to see a whole period of the offset wobble at once.

      show='gap'   x(i+k) - x(i)  -- signed, so left is negative
      show='x'     the member at i+k
      show='root'  the root at i+k, i.e. (i+k)**2 ^ m
      show='off'   the xor offset at i+k, which cycles with `xor_offsets`

    Column 0 is the centre.  Columns are labelled -reach .. +reach, and
    negative indices are dropped rather than reflected.
    """
    _need_pd()
    at = lambda j: (j * j ^ m) ** 2
    rows = []
    for i in range(lo, hi):
        row = {"i": i}
        for k in range(-reach, reach + 1):
            j = i + k
            if j < 0:
                row[k] = None
                continue
            if show == "gap":
                row[k] = at(j) - at(i)
            elif show == "x":
                row[k] = at(j)
            elif show == "root":
                row[k] = j * j ^ m
            elif show == "off":
                row[k] = (j * j ^ m) - j * j
            else:
                raise ValueError("show must be gap, x, root or off")
        rows.append(row)
    df = pd.DataFrame(rows).set_index("i")
    if show in ("gap", "x") and digits:
        df = df.map(lambda v: _brief(v, digits, 512) if isinstance(v, int) else v)
    df.attrs["m"] = m
    df.attrs["reach"] = reach
    df.attrs["show"] = show
    return df


def square_parity(lo=3, hi=24, m=5):
    """Do the spans between members ever balance?  Yes -- every other one.

    Chain the family so each row runs from its own member to the next, and
    count evens against odds inside that span.  No single span balances:
    every one is odd in length, so it always finishes one ahead.  But the
    sign alternates, so **every consecutive pair cancels exactly** and the
    running total returns to zero at every even index, forever.

    The reason is bit 0, the one bit that is exactly balanced across
    squares (see `ntpar.exact`).  Parity of the root is
    `parity(i) XOR parity(m)`, and parity of the member is parity of the
    root, so members alternate even, odd, even, odd whatever the mask is.
    Consecutive members therefore always differ in parity, every span is
    odd, and the +1 and -1 always come in pairs.

      start, end  this member and the next -- rows chain, end == next start
      span        end - start, always odd
      r_start     isqrt(start), i.e. i*i ^ m
      r_end       isqrt(end)
      squares     perfect squares strictly inside the span -- there are
                  r_end - r_start - 1 of them, and their roots are the
                  integers between, which `roots` names
      roots       "13..20", the roots of those squares
      sq_frac     squares / span, how thin the squares are in there
      par         parity of `start`
      evens, odds counted in [start, end) by closed form, not by looping
      diff        evens - odds, always +1 or -1
      cum         running diff
      at_parity   cum == 0 -- true at every even index and no other

    Counting is closed form, so this runs at any scale: the span between
    i = 10**20 and its successor is about 4*10**60 integers and the row
    still costs nothing.
    """
    _need_pd()
    at = lambda j: (j * j ^ m) ** 2
    evens_below = lambda n: (n + 1) // 2          # count of even t in [0, n)
    rows, ce, co = [], 0, 0
    for i in range(lo, hi):
        a, b = at(i), at(i + 1)
        ra, rb = i * i ^ m, (i + 1) * (i + 1) ^ m
        span = b - a
        nsq = rb - ra - 1
        ev = evens_below(b) - evens_below(a)
        od = span - ev
        ce += ev
        co += od
        rows.append({
            "i": i, "start": a, "end": b, "span": span,
            "r_start": ra, "r_end": rb, "squares": nsq,
            "roots": f"{ra + 1}..{rb - 1}" if nsq > 0 else "-",
            "sq_frac": round(nsq / span, 8) if span else None,
            "par": "E" if a % 2 == 0 else "O",
            "evens": ev, "odds": od, "diff": ev - od,
            "cum": ce - co, "at_parity": ce == co,
        })
    big = ("start", "end", "span", "evens", "odds", "r_start", "r_end", "squares")
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("i")
    df.attrs["m"] = m
    df.attrs["total_squares"] = sum(r["squares"] for r in rows)
    df.attrs["parity_at"] = [int(i) for i in df.index[df.at_parity]]
    df.attrs["all_spans_odd"] = bool(all(r["span"] % 2 == 1 for r in rows))
    df.attrs["diff_always_one"] = bool(all(abs(r["diff"]) == 1 for r in rows))
    return df


def walk_to(value, bits=64, steps=None, fast=None):
    """Hand it a big number and say where to stop.  It walks down and lands.

    `square_chain` only descends while everything stays square, so on an
    arbitrary value it stops almost at once.  This is the other walk -- the
    root spine, `v -> isqrt(v)` -- which always descends, shedding a
    remainder at each step.  Bits halve every rung, so the whole trip is
    `log2(start / target)` steps however big the start is: two million bits
    down to sixty-four is fifteen.

        ntfun.walk_to(big*big, bits=64)      # land on a machine word
        ntfun.walk_to(big*big, steps=3)      # or just take three rungs

      bits       bit length at this rung
      value      the rung itself
      rem_bits   what the step shed -- `v - isqrt(v)**2`, the half you
                 would need to climb back
      halved     bits is ceil(previous / 2), which it always is

    Give `bits` to stop at a size, or `steps` to take a fixed number of
    rungs.  `steps` wins if you pass both.  The last row is where you
    landed; `climb` or `deep_pack` are how you would get back.
    """
    _need_pd()
    mp, sq = _fast(fast)
    v = mp(value)
    rows = [{"step": 0, "bits": v.bit_length(), "value": int(v),
             "rem_bits": None, "halved": None}]
    n = 0
    while True:
        if steps is not None and n >= steps:
            break
        if steps is None and (v.bit_length() <= bits or v <= 1):
            break
        r = sq(v)
        rem = v - r * r
        n += 1
        rows.append({"step": n, "bits": r.bit_length(), "value": int(r),
                     "rem_bits": int(rem).bit_length(),
                     "halved": r.bit_length() == -(-rows[-1]["bits"] // 2)})
        v = r
    data = {c: pd.Series([x[c] for x in rows],
                         dtype=object if c == "value" else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("step")
    df.attrs["from_bits"] = rows[0]["bits"]
    df.attrs["landed_bits"] = rows[-1]["bits"]
    df.attrs["landed"] = rows[-1]["value"]
    df.attrs["steps"] = n
    df.attrs["shed_bits"] = sum(r["rem_bits"] or 0 for r in rows)
    return df


def square_chain(x, m=5, cap=64):
    """Keep going down: sqrt, xor m, sqrt, xor m, until it stops being square.

    `441 -> 21 -> 16 -> 4 -> 1 -> 1 -> 4 -> 2 -> 7`, and it only stops
    because 7 is not a square.  One round is "take the root, xor the mask";
    every member of the family survives at least one, and some survive
    several.
    """
    out = [x]
    for _ in range(cap):
        r = isqrt(x)
        if r * r != x:
            break
        q = r ^ m
        out += [r, q]
        x = q
    return out


def square_ladder(lo=0, hi=20, m=5):
    """Every member with its whole descent -- the completish view.

    `square_family` stops at the index.  This keeps taking roots and
    xoring until the ladder runs out, which shows what a single row hides:
    the family is not a flat list, it nests.  Member 22201 descends through
    149 to 144, and 144 is member 3.  So some members stand on top of
    others, and `into` names which.

      x         the member
      rounds    complete sqrt-then-xor steps before a non-square appears
      chain     the whole descent, written out
      lands_on  where it stops -- the first non-square
      into      the index of another family member the chain passes
                through, or None if it descends into fresh ground
    """
    _need_pd()
    at = lambda j: (j * j ^ m) ** 2
    known = {at(j): j for j in range(max(hi * 2, 64))}
    rows = []
    for i in range(lo, hi):
        x = at(i)
        c = square_chain(x, m)
        rounds = (len(c) - 1) // 2
        hit = next((known[v] for v in c[1:] if v in known and known[v] != i), None)
        rows.append({
            "i": i, "x": x, "rounds": rounds,
            "chain": " -> ".join(map(str, c)) if len(c) < 14 else
                     " -> ".join(map(str, c[:6])) + f" ... ({len(c)} long)",
            "lands_on": c[-1], "into": hit,
        })
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in ("x", "lands_on") else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("i")
    df.attrs["m"] = m
    df.attrs["max_rounds"] = int(df.rounds.max())
    df.attrs["nested"] = [int(j) for j in df.index[df["into"].notna()]]
    return df


def square_family_check(hi=1_643_000, m=5):
    """Brute force below a bound vs the closed form -- do they agree?

    The point is that the closed form is *complete*, not merely a subset:
    every x the raw predicate finds is emitted, and nothing else.
    """
    _need_pd()
    brute = [x for x in range(hi)
             if isqrt(x) ** 2 == x
             and isqrt(isqrt(x) ^ m) ** 2 == (isqrt(x) ^ m)]
    i = 0
    while ((i * i ^ m) ** 2) < hi or i < 4:
        i += 1
    closed = sorted({(j * j ^ m) ** 2 for j in range(i + 2)
                     if (j * j ^ m) ** 2 < hi})
    return pd.DataFrame([{
        "bound": hi,
        "m": m,
        "brute_force": len(brute),
        "closed_form": len(closed),
        "identical": brute == closed,
        "fourth_root": round(hi ** 0.25, 1),
    }]).set_index("bound")


def xor_k_near(v, k=5, window=3):
    """Nearest member of {i^2 XOR k} to v.  Returns (index, member, distance).

    ntd.xor_k_index(k).unrank(i) is i*i ^ k, so the index is the root.
    XOR by a small k moves i*i by at most k, so the nearest member sits
    within a couple of roots of isqrt(v).
    """
    from math import isqrt
    r = isqrt(v)
    best = None
    for i in range(max(r - window, 0), r + window + 1):
        m = (i * i) ^ k
        d = abs(v - m)
        if best is None or d < best[2]:
            best = (i, m, d)
    return best


def xor_to_machine(v, k=5, window=64):
    """Smallest XOR that lands v on the xor-k square machine.

    Membership is v ^ K == i*i ^ k, so K = v ^ i*i ^ k.  Minimise K over i.
    A small K means the two agree in their high bits, so the winning i is
    near isqrt(v); the window covers the search.

    Returns (K, index, member).
    """
    from math import isqrt
    r = isqrt(v)
    best = None
    for i in range(max(r - window, 0), r + window + 1):
        m = (i * i) ^ k
        K = v ^ m
        if best is None or K < best[0]:
            best = (K, i, m)
    return best


def xor_index(lo=1, hi=24, k=5, window=64, w=None):
    """What XOR takes each rule 30 value onto the xor-k square machine.

      value      the flipstrip value
      bits       its bit length (always odd, 2g+1)
      K          the smallest XOR that lands it on the machine
      K_bin      that XOR in binary
      K_bits     its bit length
      K_pop      how many bits have to move
      index      the machine index it lands on (the root)
      member     i*i ^ k, the member it becomes
    """
    _need_pd()
    w = w or ntmc.width_for(hi)
    rows = []
    for g in range(lo, hi + 1):
        v = _val(g, w)
        K, i, m = xor_to_machine(v, k, window)
        vc = (v >> (v.bit_length() // 2)) & 1        # exact middle: bits is odd
        mc = (m >> (m.bit_length() // 2)) & 1
        rows.append({
            "gen": g,
            "value": v,
            "bits": v.bit_length(),
            "centre": vc,
            "K": K,
            "K_bin": bin(K)[2:] if K else "0",
            "K_bits": K.bit_length(),
            "K_pop": bin(K).count("1"),
            "index": i,
            "member": m,
            "mem_centre": mc,
            "centre_kept": vc == mc,
        })
    big = ("value", "K", "member")
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    return pd.DataFrame(data).set_index("gen")


def _agree(a, b):
    """Fraction of shared low bits on which two values agree."""
    n = min(a.bit_length(), b.bit_length())
    if n == 0:
        return None
    m = (1 << n) - 1
    return round(1 - bin((a ^ b) & m).count("1") / n, 6)


def centre_index(n=1024, machines=None, square=True, window=3, w=None):
    """Index the centre column in each square machine, hunting a correlation.

    The centre bits are read as one integer, squared, and handed to each
    machine; what comes back is the index, and K is the XOR that would have
    to move to make the value an exact member.  Squaring is the point: the
    machines are all `T(i*i)`, so nothing that is not a square can sit near
    one, and v*v ^ k differs from a square only in the bits k touches.

    That nearness is automatic, though -- a coin flip and a run of ones
    square just as close, K_pop 0/1/2 on m2/xor2/xor5.  So the distance
    carries nothing; only the index can.  Look for `idx_bal` off a half, or
    `agree` off a half, or K_ratio well under it.

      idx_bits    bit length of the index
      idx_bal     ones/bits of the index -- 0.5 is a fair coin
      K_bits      bit length of v XOR member
      K_pop       its popcount
      K_ratio     K_pop / K_bits -- 0.5 is a fair coin
      agree       fraction of low bits where the index matches the value
      exact       K == 0, i.e. v is a member outright

    On m2, xor2, xor5 and m2+2 the index of v*v comes back as v itself, so
    `agree` is 1 by construction and there is nothing to correlate.  nibrev
    is the one machine here whose index is not the identity.
    `square=False` indexes the raw value instead.
    """
    _need_pd()
    bits = centre_bits(n - 1, w) if n else ""
    v = int(bits, 2) if bits else 0
    if square:
        v = v * v
    machines = machines or list(MACHINES)
    res = machine_index(v, machines, window)
    rows = []
    for name in machines:
        K, i, m = res[name]
        rows.append({
            "machine": name,
            "idx_bits": i.bit_length(),
            "idx_bal": round(bin(i).count("1") / max(i.bit_length(), 1), 6),
            "K_bits": K.bit_length(),
            "K_pop": bin(K).count("1"),
            "K_ratio": round(bin(K).count("1") / max(K.bit_length(), 1), 6),
            "agree": _agree(i, v),
            "exact": K == 0,
        })
    df = pd.DataFrame(rows).set_index("machine")
    df.attrs["value_bits"] = v.bit_length()
    df.attrs["centre_bits"] = len(bits)
    return df


def root_tower(n=1024, value=None, w=None):
    """Iterated isqrt -- every index that grows to the centre value.

    Each level is the m2 index of the one above, r[k+1] = isqrt(r[k]), so
    the tower is the whole chain of squares standing under v.  Bits halve
    every step, so it reaches 1 -- the rule 30 seed -- in about
    log2(bit_length) levels: ten of them for a 1024-bit column.

    The descent is lossy on its own.  isqrt discards rem = r - isqrt(r)**2,
    and those remainders add up to the full bit length of v, so nothing is
    saved by going down; keep them and `climb` rebuilds v exactly.  Drop
    them and a root names an interval instead of a number: r covers
    [r**2, (r+1)**2 - 1], which is `span` = 2r+1 values wide.

      root        the index at this level
      root_bits   its bit length
      rem         r - isqrt(r)**2, what the step threw away
      rem_bits    its bit length
      span        2*root + 1, how many values share this root
      exact_sq    whether this level is a perfect square (rem == 0)

    `value` indexes a number you already have instead of generating the
    centre column.
    """
    _need_pd()
    if value is None:
        value = int(centre_bits(n - 1, w), 2)
    rows, r = [], value
    while r > 1:
        nr = isqrt(r)
        rem = r - nr * nr
        rows.append({"level": len(rows), "root": r, "root_bits": r.bit_length(),
                     "rem": rem, "rem_bits": rem.bit_length(),
                     "span": 2 * r + 1, "exact_sq": rem == 0})
        r = nr
    rows.append({"level": len(rows), "root": r, "root_bits": r.bit_length(),
                 "rem": None, "rem_bits": None, "span": 2 * r + 1,
                 "exact_sq": None})
    big = ("root", "rem", "span")
    data = {c: pd.Series([x[c] for x in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("level")
    df.attrs["value"] = value
    df.attrs["rem_bits_total"] = sum(x["rem_bits"] or 0 for x in rows)
    return df


def climb(df):
    """Rebuild the value a `root_tower` came from, bottom to top.

    r[k] = r[k+1]**2 + rem[k], applied from the seed up.  Exact -- the
    remainders are precisely what the descent dropped.
    """
    rows = list(df.itertuples())
    v = rows[-1].root
    for row in reversed(rows[:-1]):
        v = v * v + row.rem
    return v


def save_tower(path, n=1024, machine="xor5", value=None, w=None):
    """Store the centre column as a tower of squares, and return the frame.

    The value is pushed onto `machine` first -- for xor5 that is v*v ^ 5,
    which sits on the machine exactly with index v -- and the tower of that
    is what gets written: the seed, then one remainder per level, bottom
    up, in hex.  `load_tower` climbs it back and takes the root, so the
    centre bits come off disk instead of out of rule 30.

    Hex, not decimal: int-to-string in base 10 is quadratic, and hex packs
    four bits per character instead of the one bit per byte a '0101' file
    spends.  A million centre bits store in about 500 KB.
    """
    if value is None:
        value = int(centre_bits(n - 1, w), 2)
    fwd, _ = MACHINES[machine]
    top = fwd(value)
    df = root_tower(value=top)
    rows = list(df.itertuples())
    with open(path, "w") as fh:
        fh.write("# ntfun tower v1\n")
        fh.write(f"machine {machine}\n")
        fh.write(f"value_bits {value.bit_length()}\n")
        fh.write(f"levels {len(rows) - 1}\n")
        fh.write(f"seed {rows[-1].root:x}\n")
        for row in reversed(rows[:-1]):            # bottom up, so the
            fh.write(f"{row.rem:x}\n")             # climb is a plain loop
    return df


def _read_tower(path):
    """Parse a saved tower into (header dict, remainders bottom-up)."""
    head, rems = {}, []
    keys = ("machine", "value_bits", "levels", "seed")
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            k = line.split(" ", 1)[0]
            if k in keys and k not in head:
                head[k] = line.split(" ", 1)[1]
            else:
                rems.append(int(line, 16))
    return head, rems


def load_tower(path, bits=True):
    """Climb a saved tower and take the root -- the inverse of `save_tower`.

    Returns the centre bits as a string, or the integer with `bits=False`.
    """
    head, rems = _read_tower(path)
    v = int(head["seed"], 16)
    for rem in rems:
        v = v * v + rem
    _, inv = MACHINES[head["machine"]]
    v = isqrt(max(inv(v), 0))
    if not bits:
        return v
    return bin(v)[2:].zfill(int(head["value_bits"]))


def load_frame(path):
    """The `root_tower` frame back from a saved tower -- nothing regenerated.

    Climbing already visits every level, so the roots come for free; this
    just records them on the way up and hands back the same frame
    `save_tower` returned.
    """
    _need_pd()
    head, rems = _read_tower(path)
    v = int(head["seed"], 16)
    roots = [v]
    for rem in rems:
        v = v * v + rem
        roots.append(v)
    roots.reverse()                                # top down, level 0 first
    rems = list(reversed(rems)) + [None]           # the seed has none
    rows = [{"level": k, "root": r, "root_bits": r.bit_length(),
             "rem": m, "rem_bits": None if m is None else m.bit_length(),
             "span": 2 * r + 1, "exact_sq": None if m is None else m == 0}
            for k, (r, m) in enumerate(zip(roots, rems))]
    big = ("root", "rem", "span")
    data = {c: pd.Series([x[c] for x in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("level")
    df.attrs["value"] = roots[0]
    df.attrs["machine"] = head["machine"]
    return df


def xor_tower(seed=2, depth=20, k=5):
    """A tower with no remainders at all -- the xor absorbs every step.

    Climb with v -> v*v ^ k and descend with v -> isqrt(v ^ k).  That pair
    is exact in both directions, because xor is an involution: (v*v ^ k) ^ k
    is v*v on the nose, so isqrt lands on v and never on v-1.  Nothing is
    dropped, so there is nothing to store.

    Contrast `root_tower`, where the descent is plain isqrt and each level
    sheds a remainder of half its bits.  Here the whole number is named by
    three small integers -- seed, depth, k -- however many bits it has.  The
    catch is the obvious one: only numbers built this way have the property.

      value_bits   bit length at this level -- doubles every step
      exact        isqrt(v ^ k) recovered the level below exactly
    """
    _need_pd()
    v = seed
    rows = [{"level": 0, "value": v, "value_bits": v.bit_length(),
             "exact": True}]
    for lvl in range(1, depth + 1):
        nxt = v * v ^ k
        rows.append({"level": lvl, "value": nxt, "value_bits": nxt.bit_length(),
                     "exact": isqrt(nxt ^ k) == v})
        v = nxt
    data = {c: pd.Series([x[c] for x in rows],
                         dtype=object if c == "value" else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("level")
    df.attrs["seed"] = seed
    df.attrs["k"] = k
    df.attrs["depth"] = depth
    return df


def machine_tower(seed=3, depth=12, machine="nibrev"):
    """`xor_tower` for any machine in MACHINES, not just the xor ones.

    Every machine is `i -> T(i*i)` for an invertible T, so all of them
    climb the same way and descend with `isqrt(T(v))`.  For the xor masks
    T is `^ k`; for nibrev it is the nibble reversal.  Both are
    involutions, so both are exact and nothing is dropped.

    The reason to bother is `xor_offsets`: on a square, `^ m` collapses to
    `+ constant`, so m2, m2+2, xor2 and xor5 are one lattice wearing four
    hats and their towers differ only by tiny offsets.  nibrev does not
    collapse -- it permutes bits instead of adding to them, so it is the
    only machine here whose index is not the identity in disguise, and the
    only ladder that looks different on the way up.

      value       the level
      value_bits  bit length -- doubles for the xor machines; nibrev's
                  reversal can shed leading zeros, so its width wanders
      exact       isqrt(T(v)) recovered the level below
      stalled     the level did not grow -- a fixed point of the machine

    Watch the seed with nibrev: nibrev(4) is 2, so seed 2 is a fixed point
    and the tower sits there forever.  Seed 3 climbs.  The xor machines
    have no such trap because adding a constant always moves.
    """
    _need_pd()
    fwd, inv = MACHINES[machine]
    v = seed
    rows = [{"level": 0, "value": v, "value_bits": v.bit_length(),
             "exact": True, "stalled": False}]
    for lvl in range(1, depth + 1):
        nxt = fwd(v)
        rows.append({"level": lvl, "value": nxt,
                     "value_bits": nxt.bit_length(),
                     "exact": isqrt(max(inv(nxt), 0)) == v,
                     "stalled": nxt == v})
        v = nxt
    data = {c: pd.Series([x[c] for x in rows],
                         dtype=object if c == "value" else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("level")
    df.attrs["machine"] = machine
    df.attrs["seed"] = seed
    df.attrs["all_exact"] = bool(all(r["exact"] for r in rows))
    return df


def machine_descend(v, machine="nibrev"):
    """Walk a value back to its seed through `machine`, collecting levels.

    Stops when a step stops being exact, so it doubles as a membership
    test: anything `machine_tower` built walks all the way down.
    """
    fwd, inv = MACHINES[machine]
    out = [v]
    while v > 1:
        nxt = isqrt(max(inv(v), 0))
        if nxt < 1 or nxt >= v or fwd(nxt) != v:
            break                       # no progress, or not a member
        out.append(nxt)
        v = nxt
    return out


def xor_climb(seed=2, depth=20, k=5):
    """The value at the top of `xor_tower`, without building the frame."""
    v = seed
    for _ in range(depth):
        v = v * v ^ k
    return v


def xor_descend(v, k=5):
    """Walk a value back down to its seed, collecting every level.

    Stops when a step stops being exact, so it also works as a membership
    test: a value built by `xor_climb` walks all the way to its seed.
    """
    out = [v]
    while v > 1:
        nxt = isqrt(v ^ k)
        if nxt * nxt ^ k != v:
            break
        out.append(nxt)
        v = nxt
    return out


def _deep_walk(v, floor, path, out, sq=isqrt):
    """Split v into root and remainder and recurse on both.

    Values land in the frame as plain ints even when `sq` is gmpy2's, so
    nothing downstream ever meets an mpz."""
    leaf = v.bit_length() <= floor
    out.append({"path": path, "depth": len(path) - 1, "bits": v.bit_length(),
                "value": int(v), "leaf": leaf})
    if not leaf:
        r = sq(v)
        _deep_walk(r, floor, path + "r", out, sq)
        _deep_walk(v - r * r, floor, path + "m", out, sq)


def deep_tower(n=1024, value=None, floor=32, w=None, fast=None):
    """`root_tower` taken all the way: index the remainders too.

    A flat tower sheds a remainder at every level and has to store it.  This
    recurses on the remainder as well, and on that remainder's remainder,
    until every piece is under `floor` bits -- so nothing is ever stored raw,
    only split.  The result is a binary tree whose every node is
    root**2 + remainder and whose leaves are all small.

    It does not shrink anything, and the measurement is the interesting
    part: isqrt(v) has ceil(n/2) bits and v - isqrt(v)**2 has about n/2, so
    every split conserves the total and the leaves always sum back to the
    original.  1024 bits gives 45 leaves and 1038 leaf bits; 16384 gives 703
    leaves and 16402.  Pushing `floor` lower costs more, never less.

      path    'v' for the whole value, then 'r' down a root, 'm' down a
              remainder -- so 'vrm' is the remainder of the root of v
      depth   how many splits deep
      bits    bit length at this node
      leaf    small enough to stop

    `deep_rebuild` puts it back together.
    """
    _need_pd()
    if value is None:
        value = int(centre_bits(n - 1, w), 2)
    mp, sq = _fast(fast)
    rows = []
    _deep_walk(mp(value), floor, "v", rows, sq)
    data = {c: pd.Series([x[c] for x in rows],
                         dtype=object if c == "value" else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("path")
    leaves = [x for x in rows if x["leaf"]]
    df.attrs["value"] = value
    df.attrs["floor"] = floor
    df.attrs["leaves"] = len(leaves)
    df.attrs["leaf_bits"] = sum(x["bits"] for x in leaves)
    df.attrs["ratio"] = round(df.attrs["leaf_bits"] / value.bit_length(), 4)
    return df


def deep_rebuild(df, fast=None):
    """Rebuild the value a `deep_tower` came from, leaves upward.

    Only the leaves carry numbers; every interior node is recomputed as
    root**2 + remainder.  Exact -- the splits dropped nothing.
    """
    mp, _ = _fast(fast)
    vals = {r.Index: mp(r.value) for r in df.itertuples() if r.leaf}
    inner = [r for r in df.itertuples() if not r.leaf]
    for r in sorted(inner, key=lambda x: -x.depth):
        root, rem = vals[r.Index + "r"], vals[r.Index + "m"]
        vals[r.Index] = root * root + rem
    return int(vals["v"])


_TIGHT_MAGIC = b"NTFUNTG1"
_SNUG_MAGIC = b"NTFUNSN1"


def save_snug(path, packed):
    """Write a `snug_pack` tuple.  32-byte header, then the bit blob."""
    import numpy as np
    raw, depth, value_bits, floor, lo, dw, nbits = packed
    with open(path, "wb") as fh:
        fh.write(_SNUG_MAGIC)
        fh.write(np.array([depth, value_bits, floor, lo, dw, nbits],
                          dtype="<u4").tobytes())
        fh.write(raw)
    return 1 << depth


def load_snug(path, value=False, fast=None):
    """Read a `save_snug` file.  `value=True` rebuilds instead of returning
    the tuple."""
    import numpy as np
    with open(path, "rb") as fh:
        if fh.read(len(_SNUG_MAGIC)) != _SNUG_MAGIC:
            raise ValueError(f"{path} is not an ntfun snug file")
        depth, value_bits, floor, lo, dw, nbits = np.frombuffer(fh.read(24), "<u4")
        raw = fh.read()
    packed = (raw, int(depth), int(value_bits), int(floor),
              int(lo), int(dw), int(nbits))
    return snug_unpack(packed, fast) if value else packed



def tight_pack(value, floor=512, fast=None):
    """`deep_pack` with the scaffolding removed.  Keeps no shape at all.

    Three things `deep_pack` spends bytes on that a complete tree does not
    need.  The shape mask is implied by the depth -- if every path is the
    same length, "internal or leaf" is just "am I at the bottom yet", so one
    integer replaces 8 KB.  A higher floor cuts the node count with it: 512
    gives depth 12 and 4,096 leaves where 64 gives depth 15 and 32,768, so
    the tree builds in half the time.  And the leaves pack into fixed blocks
    of floor/8 bytes instead of uint64 slots.

        deep_pack  floor 64    270,360 bytes   1.0814x the bare integer
        tight_pack floor 512   262,168 bytes   1.0487x
                   floor 4096  262,168 bytes   1.0487x, depth 9

    It buys speed and the mask, not much else -- the leaves are still the
    number, and the number does not compress.  What it cannot do is
    incomplete trees: raise the floor until the tree comes out complete, or
    stay with `deep_pack`, which carries a shape for exactly this reason.

    Returns (leaves, depth, value_bits, floor).  `deep_pack` is still the
    general one; this is the tuned one.
    """
    import numpy as np
    if floor % 8:
        raise ValueError("floor must be a multiple of 8")
    mp, sq = _fast(fast)
    width = floor // 8
    leaves, depths = [], []

    def walk(n, d):
        if n.bit_length() <= floor:
            leaves.append(int(n))
            depths.append(d)
            return
        r = sq(n)
        walk(r, d + 1)
        walk(n - r * r, d + 1)

    walk(mp(value), 0)
    depth = depths[0]
    if len(set(depths)) != 1 or len(leaves) != 1 << depth:
        raise ValueError(
            f"tree is not complete at floor={floor} "
            f"({len(leaves):,} leaves, depths {min(depths)}-{max(depths)}); "
            f"raise the floor or use deep_pack")
    buf = bytearray(len(leaves) * width)
    for i, v in enumerate(leaves):
        buf[i * width:(i + 1) * width] = v.to_bytes(width, "big")
    return (np.frombuffer(bytes(buf), dtype=np.uint8), depth,
            value.bit_length(), floor)


def tight_unpack(packed, fast=None):
    """Rebuild from `tight_pack`.  The shape comes from the depth alone.

    A complete tree needs no stored topology: recurse `depth` times, and
    whatever you land on is a leaf.  Pre-order collection means the leaves
    come off the array in exactly the order the recursion asks for them.
    """
    mp, _ = _fast(fast)
    buf, depth, _value_bits, floor = packed
    width = floor // 8
    raw = buf.tobytes() if hasattr(buf, "tobytes") else bytes(buf)
    pos = [0]

    def build(d):
        if d == 0:
            i = pos[0]
            pos[0] += 1
            return mp(int.from_bytes(raw[i * width:(i + 1) * width], "big"))
        root = build(d - 1)
        rem = build(d - 1)
        return root * root + rem

    return int(build(depth))


def _concat_bits(pairs):
    """Concatenate (value, width) pairs into one integer, divide and conquer.

    Left-folding with `acc = acc << w | v` copies the accumulator every
    step, which is quadratic once the accumulator is a megabit.  Halving
    keeps it n log n.
    """
    if not pairs:
        return 0, 0
    if len(pairs) == 1:
        return pairs[0]
    mid = len(pairs) // 2
    a, aw = _concat_bits(pairs[:mid])
    b, bw = _concat_bits(pairs[mid:])
    return (a << bw) | b, aw + bw


def snug_pack(value, floor=512, fast=None):
    """`tight_pack` with the leaf padding squeezed out too.

    A leaf averages ~488 bits and sits in a 512-bit block, so `tight_pack`
    spends 12,099 bytes on zeros.  This stores every leaf at its own width
    instead, back to back, and keeps the widths as small deltas off the
    smallest one -- they cluster inside a couple of dozen bits, so five or
    six bits each covers it.  The leading bit of a leaf is 1 by definition
    of bit_length, so that is dropped too and re-added on the way back.

        tight_pack floor 512   262,168 bytes   1.0487x
        snug_pack  floor 512   252,7xx bytes   1.011x

    The cost is that leaves are no longer a fixed-width array -- no
    `np.frombuffer`, no `tight_leaves` view without unpacking first.  Use
    `tight_pack` when you want to look at it, `snug_pack` when you want it
    small.  Complete trees only, same as `tight_pack`.

    Returns (blob, depth, value_bits, floor, lo, dw).
    """
    mp, sq = _fast(fast)
    leaves, depths = [], []

    def walk(n, d):
        if n.bit_length() <= floor:
            leaves.append(int(n))
            depths.append(d)
            return
        r = sq(n)
        walk(r, d + 1)
        walk(n - r * r, d + 1)

    walk(mp(value), 0)
    depth = depths[0]
    if len(set(depths)) != 1 or len(leaves) != 1 << depth:
        raise ValueError(
            f"tree is not complete at floor={floor} "
            f"({len(leaves):,} leaves, depths {min(depths)}-{max(depths)}); "
            f"raise the floor or use deep_pack")
    widths = [v.bit_length() for v in leaves]
    lo, hi = min(widths), max(widths)
    dw = max((hi - lo).bit_length(), 1)
    pairs = [(w - lo, dw) for w in widths]
    pairs += [(v & ((1 << (w - 1)) - 1) if w else 0, max(w - 1, 0))
              for v, w in zip(leaves, widths)]
    blob, nbits = _concat_bits(pairs)
    raw = blob.to_bytes((nbits + 7) // 8, "big") if nbits else b""
    return (raw, depth, value.bit_length(), floor, lo, dw, nbits)


def snug_unpack(packed, fast=None):
    """Rebuild from `snug_pack` -- widths first, then the bodies."""
    mp, _ = _fast(fast)
    raw, depth, _value_bits, _floor, lo, dw, nbits = packed
    n = 1 << depth
    blob = int.from_bytes(raw, "big")           # leading pad does not shift it
    m = (1 << dw) - 1
    widths = [((blob >> (nbits - (i + 1) * dw)) & m) + lo for i in range(n)]
    leaves, cur = [], nbits - n * dw            # bodies sit under the widths
    for w in widths:
        b = max(w - 1, 0)
        cur -= b
        tail = (blob >> cur) & ((1 << b) - 1) if b else 0
        leaves.append(((1 << (w - 1)) | tail) if w else 0)
    pos_i = [0]

    def build(d):
        if d == 0:
            v = leaves[pos_i[0]]
            pos_i[0] += 1
            return mp(v)
        root = build(d - 1)
        rem = build(d - 1)
        return root * root + rem

    return int(build(depth))


def tight_leaves(packed, frame=True):
    """The leaves of a `tight_pack` as numbers instead of raw bytes.

    `deep_pack` hands back a uint64 array you can read straight off --
    `2362801676638030348` and so on.  A 512-bit leaf fits no machine type,
    so `tight_pack` stores fixed byte blocks and the array prints as bytes.
    This cuts the blocks back into integers.

    Note the padding: a leaf averages ~488 bits in a 512-bit block, so the
    first two or three bytes of each block are the zeros you see in the
    raw array.

        ntfun.show(ntfun.tight_leaves(t)).head()
        ntfun.tight_leaves(t, frame=False)[:3]     # plain list of ints
    """
    buf, depth, _value_bits, floor = packed
    width = floor // 8
    raw = buf.tobytes() if hasattr(buf, "tobytes") else bytes(buf)
    vals = [int.from_bytes(raw[i * width:(i + 1) * width], "big")
            for i in range(1 << depth)]
    if not frame:
        return vals
    _need_pd()
    rows = [{"i": i, "bits": v.bit_length(), "value": v}
            for i, v in enumerate(vals)]
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c == "value" else None)
            for c in rows[0]}
    df = pd.DataFrame(data).set_index("i")
    df.attrs["floor"] = floor
    df.attrs["padding_bits"] = (1 << depth) * floor - int(df.bits.sum())
    return df


def save_tight(path, packed):
    """Write a `tight_pack` tuple.  Header is 24 bytes and then the leaves."""
    import numpy as np
    buf, depth, value_bits, floor = packed
    with open(path, "wb") as fh:
        fh.write(_TIGHT_MAGIC)
        fh.write(np.array([depth, value_bits, floor, 0], dtype="<u4").tobytes())
        fh.write(buf.tobytes() if hasattr(buf, "tobytes") else bytes(buf))
    return 1 << depth


def load_tight(path, value=False, fast=None):
    """Read a `save_tight` file.  `value=True` rebuilds instead of returning
    the tuple."""
    import numpy as np
    with open(path, "rb") as fh:
        if fh.read(len(_TIGHT_MAGIC)) != _TIGHT_MAGIC:
            raise ValueError(f"{path} is not an ntfun tight file")
        depth, value_bits, floor, _ = np.frombuffer(fh.read(16), "<u4")
        buf = np.frombuffer(fh.read(), dtype=np.uint8)
    packed = (buf, int(depth), int(value_bits), int(floor))
    return tight_unpack(packed, fast) if value else packed


_LEAF_MAGIC = b"NTFUNLV1"


def save_leaves(path, n=1024, value=None, floor=64, w=None, fast=None):
    """Write a `deep_tower` as raw uint64 leaves -- its smallest honest form.

    Every leaf is under `floor` bits, so at floor 64 they are all machine
    words and the whole thing is a flat array plus the tree's shape.  The
    shape is one bit per node in pre-order (1 internal, 0 leaf), which for
    a million-bit column is 8 KB against 256 KB of leaves.

    Layout: magic, four uint32 (nodes, leaves, value_bits, floor), the
    packed shape, then the leaves little-endian in pre-order.  No hex, no
    CSV -- `load_leaves` reads it back with numpy and a loop.
    """
    if value is None:
        value = int(centre_bits(n - 1, w), 2)
    packed = deep_pack(value, floor, fast)
    return save_pack(path, packed, value.bit_length(), floor)


def load_leaves(path, bits=False, fast=None):
    """Rebuild the value from a `save_leaves` file -- squares and adds only."""
    import numpy as np
    with open(path, "rb") as fh:
        if fh.read(len(_LEAF_MAGIC)) != _LEAF_MAGIC:
            raise ValueError(f"{path} is not an ntfun leaf file")
        n_nodes, n_leaves, value_bits, _floor = np.frombuffer(
            fh.read(16), dtype="<u4")
        shape = np.unpackbits(np.frombuffer(
            fh.read((int(n_nodes) + 7) // 8), dtype=np.uint8))[:int(n_nodes)]
        leaves = np.frombuffer(fh.read(int(n_leaves) * 8), dtype="<u8")
    pos = [0, 0]                                   # node cursor, leaf cursor

    def build():
        internal = shape[pos[0]]
        pos[0] += 1
        if not internal:
            v = int(leaves[pos[1]])
            pos[1] += 1
            return v
        root = build()
        rem = build()
        return root * root + rem

    v = build()
    return bin(v)[2:].zfill(int(value_bits)) if bits else v


def deep_pack(value, floor=64, fast=None):
    """A deep tower held in memory with nothing wasted -- the store, not the view.

    `deep_tower` keeps every node as a python int in a frame, which is what
    you want when you are looking at one: 9 MB resident and a 29 MB peak for
    a 2,000,002-bit value.  This walks the same split and keeps only what
    cannot be recomputed -- the shape as packed bits, the leaves as uint64 --
    so the interior nodes live on the recursion stack and die there.  Same
    value costs 0.27 MB resident, 7 MB peak.

    That is 34x, and it grows with the input: the frame is linear in the
    value, so ten million bits is a ~45 MB frame against ~1.3 MB packed.

    Returns (shape, leaves, n_nodes) -- hand the whole tuple to
    `deep_unpack`.  Same layout `save_leaves` writes, without the file.

        p = ntfun.deep_pack(big*big ^ 5)
        ntfun.deep_unpack(p) == big*big ^ 5      # True
    """
    import numpy as np
    if floor > 64:
        raise ValueError("floor must be <= 64 for uint64 leaves")
    mp, sq = _fast(fast)
    shape, leaves = [], []

    def walk(n):
        if n.bit_length() <= floor:
            shape.append(0)
            leaves.append(int(n))
            return
        shape.append(1)
        r = sq(n)
        walk(r)
        walk(n - r * r)

    walk(mp(value))
    return (np.packbits(np.array(shape, np.uint8)),
            np.array(leaves, dtype=np.uint64), len(shape))


def deep_unpack(packed, fast=None):
    """Rebuild the value from `deep_pack` -- squares and adds, leaves upward.

    Peak memory is about 4 MB for a 2,000,002-bit value and that part is
    irreducible: the top nodes have to exist as integers to be added, and
    the root alone is 250 KB.
    """
    import numpy as np
    mp, _ = _fast(fast)
    shape_bits, leaves, n_nodes = packed
    shape = np.unpackbits(shape_bits)[:n_nodes]
    pos = [0, 0]

    def build():
        internal = shape[pos[0]]
        pos[0] += 1
        if not internal:
            v = mp(int(leaves[pos[1]]))
            pos[1] += 1
            return v
        root = build()
        rem = build()
        return root * root + rem

    return int(build())


def _span(shape, i):
    """(nodes, leaves) in the pre-order subtree starting at position i.

    Every node has two children or none, so one counter suffices: an
    internal node consumes its slot and owes two more, a leaf pays one off.
    """
    need, n, lv = 1, 0, 0
    while need:
        internal = shape[i + n]
        n += 1
        if internal:
            need += 1
        else:
            need -= 1
            lv += 1
    return n, lv


def deep_subtree(packed, path):
    """One node's subtree, pulled out of a `deep_pack` without rebuilding it.

    `deep_pack` walks pre-order, so a subtree is a contiguous run in both
    arrays at once -- nodes and leaves slice together, no gather.  The
    result is a packed tuple in its own right, so `deep_unpack` takes it
    and the cost is proportional to the subtree, not the whole tower:
    0.15 s for the root, 0.02 s three levels down, 0 s at the bottom.

        p = ntfun.deep_pack(big*big ^ 5)
        ntfun.deep_unpack(ntfun.deep_subtree(p, 'vrrm'))   # that node's value

    `path` follows the `deep_tower` convention -- 'v', then 'r' into a root
    and 'm' into a remainder.

    Worth being clear about what this is not: a subtree is NOT a region of
    the underlying bits.  `isqrt` is global, so 'vr' is the square root of
    the whole value, not its first half, and flipping one bit of the input
    changes every node.  The tower has no positional locality to index.
    Slice the bit string if that is what you want.
    """
    import numpy as np
    shape_bits, leaves, n_nodes = packed
    shape = np.unpackbits(shape_bits)[:n_nodes]
    i = lf = 0
    for step, c in enumerate(path[1:], 1):
        if not shape[i]:
            raise ValueError(f"{path[:step]!r} is a leaf -- cannot descend")
        i += 1                                   # 'r' child is the next node
        if c == "m":
            n, lv = _span(shape, i)              # skip the whole 'r' subtree
            i += n
            lf += lv
        elif c != "r":
            raise ValueError(f"bad step {c!r} in {path!r} -- use 'r' or 'm'")
    n, lv = _span(shape, i)
    return (np.packbits(shape[i:i + n]), leaves[lf:lf + lv], n)


def save_pack(path, packed, value_bits=0, floor=64):
    """Write a `deep_pack` tuple straight to disk, no rebuilding.

    Same file `save_leaves` produces -- `load_leaves` reads it back as the
    value, `load_pack` reads it back as the tuple.  `value_bits` is only
    used by `load_leaves(bits=True)` to left-pad; leave it at 0 if you do
    not have it and the padding is skipped.

        p = ntfun.deep_pack(big*big ^ 5)
        ntfun.save_pack('centre_1M.leaves', p, value_bits=(big*big ^ 5).bit_length())
    """
    import numpy as np
    shape_bits, leaves, n_nodes = packed
    with open(path, "wb") as fh:
        fh.write(_LEAF_MAGIC)
        fh.write(np.array([n_nodes, len(leaves), value_bits, floor],
                          dtype="<u4").tobytes())
        fh.write(np.asarray(shape_bits, dtype=np.uint8).tobytes())
        fh.write(np.asarray(leaves).astype("<u8").tobytes())
    return len(leaves)


def load_pack(path):
    """Read a leaf file back as a `deep_pack` tuple, without rebuilding.

    Use this when you want the packed form in memory -- 0.27 MB for a
    2,000,002-bit value.  `load_leaves` goes straight to the integer
    instead, which costs the 4 MB of intermediates.
    """
    import numpy as np
    with open(path, "rb") as fh:
        if fh.read(len(_LEAF_MAGIC)) != _LEAF_MAGIC:
            raise ValueError(f"{path} is not an ntfun leaf file")
        n_nodes, n_leaves, value_bits, _floor = np.frombuffer(fh.read(16), "<u4")
        shape_bits = np.frombuffer(fh.read((int(n_nodes) + 7) // 8), np.uint8)
        leaves = np.frombuffer(fh.read(int(n_leaves) * 8), "<u8")
    return (shape_bits, leaves, int(n_nodes))


def _brief(n, digits=30, max_bits=20000):
    """A long integer as head...tail, with its size.  Decimal while that is
    cheap, hex above it -- int to decimal is quadratic and a 350,000-digit
    conversion is not worth waiting for just to look at a frame."""
    if n is None or not isinstance(n, int):
        return n
    if n.bit_length() <= 64:
        return str(n)
    if n.bit_length() <= max_bits:
        s = str(n)
        if len(s) <= 2 * digits + 3:
            return s
        return f"{s[:digits]}...{s[-digits:]} ({len(s):,}d)"
    s = format(n, "x")
    return f"0x{s[:digits]}...{s[-digits:]} ({n.bit_length():,}b)"


def show(df, digits=30, max_bits=20000):
    """Display-safe copy of a tower frame -- big integers abbreviated.

    pandas prints object-dtype ints in full, so a single huge cell pads its
    whole column and the frame turns into megabytes of whitespace.  This
    swaps those cells for head...tail summaries so the shape stays readable.

        ntfun.show(ntfun.xor_tower(2, 21))
        ntfun.show(ntfun.root_tower(1024))
    """
    _need_pd()
    out = df.copy()
    for c in out.columns:
        if out[c].dtype == object:
            out[c] = [_brief(x, digits, max_bits) for x in out[c]]
    return out


def _tree_rows(df, max_depth, digits):
    """Walk a `deep_tower` top-down, collecting drawable rows."""
    rows = {r.Index: r for r in df.itertuples()}
    under = {}
    for r in sorted(rows.values(), key=lambda x: -x.depth):
        under[r.Index] = (1 if r.leaf else
                          under.get(r.Index + "r", 0) + under.get(r.Index + "m", 0))
    out = []

    def walk(path, prefix, last):
        r = rows.get(path)
        if r is None:
            return
        if path == "v":
            head, tag = "", "value"
        else:
            head = "`- " if last else "|- "
            tag = "root" if path[-1] == "r" else "rem"
        cut = r.depth >= max_depth and not r.leaf
        out.append({"path": path, "shape": f"{prefix}{head}{tag}",
                    "bits": r.bits, "under": under[path], "leaf": bool(r.leaf),
                    "cut": cut,
                    "value": _brief(r.value, digits, 512) if r.bits <= 64 else ""})
        if cut or r.leaf:
            return
        kid = prefix + ("" if path == "v" else ("   " if last else "|  "))
        walk(path + "r", kid, False)
        walk(path + "m", kid, True)

    walk("v", "", True)
    return out


def tree(df, max_depth=4, digits=14):
    """A `deep_tower` as an indented tree -- narrow, whatever the size.

    Returns a string, so print it.  `tree_frame` is the same walk as a
    DataFrame if you would rather look at one of those.

        print(ntfun.tree(d))
        print(ntfun.tree(d, max_depth=6))
    """
    out = []
    for r in _tree_rows(df, max_depth, digits):
        note = f"  ... {r['under']:,} leaves" if r["cut"] else ""
        out.append(f"{r['shape']:<24}{r['bits']:>11,}b  {r['value']}{note}")
    return "\n".join(out)


def tree_frame(df, max_depth=4, digits=14):
    """`tree` as a DataFrame -- the drawing kept in a column.

    The shape column is padded to a fixed width so pandas cannot right-align
    the indentation away.  `under` is how many leaves sit beneath each node,
    which is the only thing you lose by cutting the walk off early.

        ntfun.tree_frame(d)
        ntfun.tree_frame(d, max_depth=6)
    """
    _need_pd()
    rows = _tree_rows(df, max_depth, digits)
    w = max(len(r["shape"]) for r in rows)
    for r in rows:
        r["shape"] = r["shape"].ljust(w)
        del r["cut"]
    return pd.DataFrame(rows).set_index("path")


def rungs(df, max_bits=256):
    """The readable bottom of a tower, as the climb written out.

    The top levels are hundreds of thousands of digits and print as a wall
    of padding; this keeps the rows that fit on a screen and shows each
    step as the equation that makes it.
    """
    _need_pd()
    rows = list(df.itertuples())
    out = []
    for k in range(len(rows) - 1, 0, -1):
        up = rows[k - 1]
        if up.root_bits > max_bits:
            break
        out.append({"level": up.Index, "bits": up.root_bits,
                    "step": f"{rows[k].root}^2 + {up.rem} = {up.root}"})
    return pd.DataFrame(out).set_index("level")


def halves(lo=1, hi=24, w=None):
    """Split each value at its centre cell and flipstrip both wings.

    A value spans 2g+1 bits, so the centre is exact and the two wings are
    the same length.  Each wing gets the same treatment as the whole --
    reversed and stripped -- so it is a width-free number in its own right.

      value        the flipstrip value
      bits         2g+1
      hi_bin       the wing above the centre
      centre       the centre cell
      lo_bin       the wing below the centre
      hi_val       flipstrip of the high wing
      lo_val       flipstrip of the low wing
      hi_pop       live cells in the high wing
      lo_pop       live cells in the low wing
      symmetric    do the two wings agree bit for bit
      xor_wings    hi_val ^ lo_val, zero iff they agree after stripping
    """
    _need_pd()
    w = w or ntmc.width_for(hi)
    rows = []
    for g in range(lo, hi + 1):
        v = _val(g, w)
        s = bin(v)[2:]
        n = len(s)
        c = n // 2
        hs, mid, ls = s[:c], s[c], s[c + 1:]
        hv = ntmc.flipstrip(int(hs, 2)) if hs.strip("0") else 0
        lv = ntmc.flipstrip(int(ls, 2)) if ls.strip("0") else 0
        rows.append({
            "gen": g,
            "value": v,
            "bits": n,
            "hi_bin": hs,
            "centre": int(mid),
            "lo_bin": ls,
            "hi_val": hv,
            "lo_val": lv,
            "hi_pop": hs.count("1"),
            "lo_pop": ls.count("1"),
            "symmetric": hs == ls[::-1],
            "xor_wings": hv ^ lv,
        })
    big = ("value", "hi_val", "lo_val", "xor_wings")
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    return pd.DataFrame(data).set_index("gen")


SQUARE_FORMS = ("m2", "m2+2", "m2^2", "m2^5")


def _target(i, form):
    s = i * i
    if form == "m2":
        return s
    if form == "m2+2":
        return s + 2
    if form == "m2^2":
        return s ^ 2
    if form == "m2^5":
        return s ^ 5
    raise ValueError(form)


def near_square(v, forms=SQUARE_FORMS, window=4):
    """Nearest square-family member to v by XOR distance.

    Searches m^2, m^2+2, m^2 XOR 2 and m^2 XOR 5 over indices around
    isqrt(v).  Note m^2+2 and m^2 XOR 2 coincide for every square -- bit 1
    is never set in a square -- so they agree wherever the target is
    exactly a square plus two.

    Returns (K, index, form, target, arithmetic_distance).
    """
    from math import isqrt
    r = isqrt(v)
    best = None
    for i in range(max(r - window, 0), r + window + 1):
        for f in forms:
            t = _target(i, f)
            K = v ^ t
            if best is None or K < best[0]:
                best = (K, i, f, t, abs(v - t))
    return best


def regimes(lo=1, hi=24, w=None):
    """Score each generation's three regions, and hunt squares in the wings.

    Under the flipstrip convention a value has three distinct regimes:

      hi wing   periodic -- the tip repeats with period 4, verified 20
                cells deep; cells are functions of g mod 1, 2 or 4
      lo wing   the light-cone edge, carrying the ruler law -- structured
                but not periodic; needs period 64 by depth nine
      centre    the chaotic middle, where the centre column lives

    For each wing it reports the nearest member of the square family by
    XOR distance, so `K == 0` marks an exact hit.

      hi_val / lo_val     flipstrip of each wing
      hi_K / lo_K         XOR distance to the nearest square-family member
      hi_form / lo_form   which form matched
      hi_idx / lo_idx     the index m
      hi_exact / lo_exact K == 0
    """
    _need_pd()
    w = w or ntmc.width_for(hi)
    rows = []
    for g in range(lo, hi + 1):
        v = _val(g, w)
        s = bin(v)[2:]
        n = len(s)
        c = n // 2
        hs, mid, ls = s[:c], s[c], s[c + 1:]
        hv = ntmc.flipstrip(int(hs, 2)) if hs.strip("0") else 0
        lv = ntmc.flipstrip(int(ls, 2)) if ls.strip("0") else 0
        d = {"gen": g, "bits": n, "centre": int(mid),
             "hi_val": hv, "lo_val": lv}
        for tag, x in (("hi", hv), ("lo", lv)):
            if x < 1:
                d.update({f"{tag}_K": None, f"{tag}_pop": None,
                          f"{tag}_idx": None, f"{tag}_form": None,
                          f"{tag}_exact": None})
                continue
            K, i, f, t, _ = near_square(x)
            d.update({f"{tag}_K": K, f"{tag}_pop": bin(K).count("1"),
                      f"{tag}_idx": i, f"{tag}_form": f,
                      f"{tag}_exact": K == 0})
        rows.append(d)
    big = ("hi_val", "lo_val", "hi_K", "lo_K")
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    return pd.DataFrame(data).set_index("gen")


def _nibrev(n):
    """Reverse each 4-bit block in place -- the nibrev machine's transform."""
    out = sh = 0
    while n:
        b = n & 15
        out |= (((b & 1) << 3) | ((b & 2) << 1) |
                ((b & 4) >> 1) | ((b & 8) >> 3)) << sh
        n >>= 4
        sh += 4
    return out


MACHINES = {
    "m2":     (lambda i: i * i,          lambda n: n),
    "m2+2":   (lambda i: i * i + 2,      lambda n: n - 2 if n >= 2 else 0),
    "xor2":   (lambda i: (i * i) ^ 2,    lambda n: n ^ 2),
    "xor5":   (lambda i: (i * i) ^ 5,    lambda n: n ^ 5),
    "nibrev": (lambda i: _nibrev(i * i), _nibrev),
}


def machine_index(v, machines=None, window=3):
    """Nearest machine member to v, by XOR distance.  Cheap for huge v.

    Every machine here is `i -> T(i*i)` for an invertible T, so the index
    is recovered by isqrt(T(v)) -- one big-integer square root.  On a
    200,001-bit value that costs ~22 ms for the plain machines and ~125 ms
    for nibrev.  The index lands within a step or two of the truth, so a
    small window around it finds the nearest member exactly.

    Returns dict of {name: (K, index, member)}.
    """
    from math import isqrt
    machines = machines or list(MACHINES)
    out = {}
    for name in machines:
        fwd, inv = MACHINES[name]
        r = isqrt(max(inv(v), 0))
        best = None
        for i in range(max(r - window, 0), r + window + 1):
            m = fwd(i)
            K = v ^ m
            if best is None or K < best[0]:
                best = (K, i, m)
        out[name] = best
    return out


def index_frame(lo=1, hi=20, wing=None, machines=None, w=None):
    """Nearest machine member per generation, with the index and XOR cost.

    `wing` selects what to index: None for the whole value, 'hi' or 'lo'
    for one side of the centre split.

      index      the machine index, essentially isqrt of the untransformed
                 value -- computable directly, no search
      K          XOR distance to that member
      K_pop      bits that must move
      best       which machine sits nearest
    """
    _need_pd()
    w = w or ntmc.width_for(hi)
    machines = machines or list(MACHINES)
    rows = []
    for g in range(lo, hi + 1):
        v = _val(g, w)
        if wing:
            s = bin(v)[2:]
            c = len(s) // 2
            part = s[:c] if wing == "hi" else s[c + 1:]
            v = ntmc.flipstrip(int(part, 2)) if part.strip("0") else 0
        if v < 1:
            continue
        res = machine_index(v, machines)
        d = {"gen": g, "value": v, "bits": v.bit_length()}
        for name in machines:
            K, i, m = res[name]
            d[f"{name}_K"] = K
            d[f"{name}_pop"] = bin(K).count("1")
            d[f"{name}_idx"] = i
        d["best"] = min(machines, key=lambda n: res[n][0])
        d["best_K"] = res[d["best"]][0]
        d["exact"] = d["best_K"] == 0
        rows.append(d)
    big = tuple(["value", "best_K"] + [f"{n}_K" for n in machines] +
                [f"{n}_idx" for n in machines])
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    return pd.DataFrame(data).set_index("gen")


def xor5(lo=1, hi=40, k=5, w=None):
    """Rule 30 values against the xor-k square machine.

      value      the flipstrip value
      bits       its bit length (always odd: 2g+1)
      sq_leak    distance to the nearest perfect square
      xk_index   which index of {i^2 XOR k} it lands nearest
      xk_value   that member
      xk_leak    distance to it
      closer     which machine it sits nearer to
    """
    _need_pd()
    from math import isqrt
    w = w or ntmc.width_for(hi)
    rows = []
    for g in range(lo, hi + 1):
        v = _val(g, w)
        i, m, d = xor_k_near(v, k)
        sl = leak(v)
        rows.append({
            "gen": g,
            "value": v,
            "bits": v.bit_length(),
            "sq_root": isqrt(v),
            "sq_leak": sl,
            "xk_index": i,
            "xk_value": m,
            "xk_leak": d,
            "closer": "xor%d" % k if d < sl else ("square" if sl < d else "tie"),
        })
    big = ("value", "xk_value")
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in big else None)
            for c in rows[0]}
    return pd.DataFrame(data).set_index("gen")


# =====================================================================
# centre column, laid out for reading and pasting
# =====================================================================
def centres(n=1024, per_line=128, label=True, w=None):
    """The rule 30 centre column as a block of text.

    Generation 0 (the seed) is the leading 1.  `per_line` chooses the
    wrap; 128 gives 8 lines for 1024 bits, 64 gives 16.
    """
    bits = centre_bits(n, w)
    out = []
    for i in range(0, len(bits), per_line):
        chunk = bits[i:i + per_line]
        out.append(f"{i:>7}  {chunk}" if label else chunk)
    return "\n".join(out)


def centre_bits(n=1024, w=None):
    """The centre column as a plain string, no formatting.

    Streamed, so n in the hundreds of thousands costs time and not
    memory.  Note the length is n+1: generation 0 is the leading 1.
    """
    return ntmc.center_string(n, w)


def centre_stats(n=4096, w=None):
    """Balance and run statistics for the centre column.

    The prize questions are whether it is periodic and whether it is
    balanced; these are the descriptive numbers, not an answer to either.
    """
    _need_pd()
    s = centre_bits(n, w)
    runs, cur = [], 1
    for a, b in zip(s, s[1:]):
        if a == b:
            cur += 1
        else:
            runs.append(cur)
            cur = 1
    runs.append(cur)
    ones = s.count("1")
    return pd.DataFrame([{
        "bits": len(s),
        "ones": ones,
        "zeros": len(s) - ones,
        "balance": round(ones / len(s), 6),
        "longest_run": max(runs),
        "mean_run": round(sum(runs) / len(runs), 4),
        "distinct_16mers": len({s[i:i + 16] for i in range(len(s) - 15)}),
        "of_possible": 1 << 16,
    }]).set_index("bits")


# =====================================================================
# everything at once
# =====================================================================
def findings(gens=2000):
    """Every result from the exploration, as a dict of frames.

        f = ntfun.findings()
        f['edge_law']      the v2(g+1) -> edge-run table
        f['parity']        mod 4, the k with v^k == v+k, square candidates
        f['periods']       low bits as functions of generation mod 2^j
        f['gaps']          gap distribution vs geometric p=1/2
        f['spikes']        zero runs, flagged at powers of two
        f['squares']       values that are perfect squares
    """
    _need_pd()
    return {
        "edge_law": edge_law(1, min(gens, 4095)),
        "parity": parity(1, min(gens, 40)),
        "periods": periods(16, gens),
        "gaps": gaps(1, min(gens, 500)),
        "spikes": spikes(max(1, gens - 24), gens),
        "squares": squares(1, gens),
    }


def nest_tower(seed=7, m=5, height=6):
    """Build a NESTED family tower: feed each member back in as the next index.

    One member is `(k**2 ^ m)**2`.  Do it again to the whole member, and again --
    every level is a genuine member of the level below, so the tower is real all
    the way up.  Bits roughly quadruple per level, so a small seed reaches
    thousands of bits fast: seed 7, height 6 -> an 11,675-bit top.

    Walk it back down with `nest_walk`; the descent lands exactly on `seed`.

      level    0 is the seed, `height` is the top
      value    the member at this level (Python int, exact)
      bits      its bit length
      digits    its decimal length
    """
    _need_pd()
    fam = lambda k: (k * k ^ m) ** 2
    v, rows = int(seed), []
    for lvl in range(height + 1):
        rows.append({"level": lvl, "value": v,
                     "bits": v.bit_length(), "digits": len(str(v))})
        v = fam(v)
    df = pd.DataFrame(
        {c: pd.Series([r[c] for r in rows],
                      dtype=object if c == "value" else None) for c in rows[0]}
    ).set_index("level")
    df.attrs["m"] = m
    df.attrs["seed"] = int(seed)
    df.attrs["top"] = rows[-1]["value"]
    df.attrs["top_bits"] = rows[-1]["bits"]
    return df


def nest_walk(top, m=5, fast=None):
    """Square-walk a nested member all the way down to its seed index.

    Each level of the tower is `(k**2 ^ m)**2`, so ONE step down is two integer
    square roots with the mask xored between them -- undo the square, undo the
    mask, undo the square.  (That is the difference from `square_chain`, which
    xors after every root and so only clears one flat level.)  It descends until
    a level is no longer a nested member, and reports where it landed.

        ntfun.nest_walk(ntfun.nest_tower(7).attrs["top"])   # -> lands on 7

      level        the level it just dropped ONTO (top-1 down to 0)
      root_bits    bit length of the outer root at this level
      value        the index recovered at this level
      exact        both roots were perfect squares (a real nested member)
    """
    _need_pd()
    mp, sq = _fast(fast)
    v = mp(top)
    rows, lvl = [], 0
    while True:
        r = sq(v)
        if r * r != v:                     # outer not a square -> not a member
            break
        q = r ^ m
        k = sq(q)
        if k * k != q:                     # inner (after xor) not a square
            break
        lvl += 1
        rows.append({"level": lvl, "root_bits": int(r).bit_length(),
                     "value": int(k), "exact": True})
        if k <= 1:
            break
        v = k
    # number the rows so the last one is level 0 (the seed)
    depth = len(rows)
    for i, row in enumerate(rows):
        row["level"] = depth - 1 - i
    df = pd.DataFrame(
        {c: pd.Series([r[c] for r in rows],
                      dtype=object if c == "value" else None) for c in rows[0]}
    ).set_index("level") if rows else pd.DataFrame()
    df.attrs["m"] = m
    df.attrs["height"] = depth
    df.attrs["landed"] = int(rows[-1]["value"]) if rows else int(top)
    return df


def lift_root(s, bits, tag=None):
    """Square roots of `s` modulo `2**bits`, built from the LOW bit up.

    The mirror of `isqrt`: instead of one top-down root, this lifts the whole
    solution set one bit at a time (`mod 2**j` -> `mod 2**(j+1)`).  A square has
    SEVERAL roots mod `2**bits`, so the rebuild BRANCHES -- that branch tree is
    the result.  `tag=None` returns every root (sorted); an int `tag` returns
    just that branch (mod the count), or None if `s` is not a square there.
    Verify any root by squaring it back.
    """
    s &= (1 << bits) - 1
    roots = {0}
    for j in range(bits):
        mod = (1 << (j + 1)) - 1
        roots = {c for r in roots for c in (r, r | (1 << j))
                 if (c * c) & mod == s & mod}
        if not roots:
            return None if tag is not None else []
    ordered = sorted(roots)
    return ordered if tag is None else ordered[tag % len(ordered)]


def lift_index(x, m=5, bits=None):
    """Rebuild the index of a `(k**2 ^ m)**2` member from its LOW `bits` bits.

    The mirror of the top-down inverse `isqrt(isqrt(x) ^ m)`: lift the outer
    root with `lift_root`, undo the mask, lift the inner root, and return EVERY
    index residue `k` (mod `2**bits`) whose member agrees with `x` in those low
    bits.  The true index is always in the set; squaring each candidate forward
    (`(k*k ^ m)**2 == x`) picks it out.

    Reaches where the top-down inverse cannot: a truncated / streamed `x` (only
    its low bits known) still pins the index's low bits to this small branch
    set.  `bits` defaults to the full width, where one branch is the exact index.
    """
    x = int(x)
    if bits is None:
        bits = x.bit_length() + 2
    mask = (1 << bits) - 1
    out = set()
    for r in lift_root(x, bits):
        out.update(lift_root((r ^ int(m)) & mask, bits))
    return sorted(out)


def lift_verify(x, m=5):
    """Low-bits-up recovery, resolved to one answer: rebuild the index branches
    with `lift_index`, square each forward, and return the single `k` with
    `(k*k ^ m)**2 == x` -- or None if `x` is not a `(k**2 ^ m)**2` member."""
    x = int(x)
    for k in lift_index(x, m):
        if (k * k ^ int(m)) ** 2 == x:
            return k
    return None


def qnr2_filter(hi, tight=True, rows=False):
    """Lars Rocha's Mersenne-offset prime filter -- fully modular, ~85x faster.

    Scans odd x in [3, hi) using  offset = (2**x - 1) - 2**(x//2)  and keeps x with
    `offset % x == 2`.  By Fermat's little theorem (2**x == 2 mod x) plus Euler's
    criterion (2**((x-1)/2) == the Legendre symbol (2|x) == +-1), for prime x this
    holds exactly when 2 is a quadratic NON-residue mod x -- i.e. x == 3 or 5 (mod 8).

    So condition 1 alone is a base-2 Euler-Jacobi primality test; it lets a few Euler
    pseudoprimes through (325, 3277, 3531, 3725, 6063, 29341, 31525, 49141, ...).
    `tight=True` (default) ALSO requires  x % (offset % (x//2)) == 2, which removes
    every one of those pseudoprimes -- empirically CLEAN (no composite passes) to at
    least 3e7.  The tight filter keeps a slowly-decaying slice of the non-residue
    primes (~4.6% near 1e4 -> ~2.1% near 1e7) and tilts ~2.9:1 toward x == 3 (mod 8),
    a parity-of-(x-1)/2 effect.

    Never materialises 2**x (uses pow(2, e, m)).  Returns a list of matching x, or,
    with rows=True, the original 7-column diagnostic tuples
    (x, (off%x)+1, x%(off%r), x%((off*off-2)%r), ((s%off)-1)%x, isprime(3), isprime(x)).
    """
    out = []
    for x in range(3, int(hi), 2):
        r = x // 2
        off_x = (pow(2, x, x) - 1 - pow(2, r, x)) % x
        if off_x != 2:
            continue
        off_r = (pow(2, x, r) - 1 - pow(2, r, r)) % r
        if not tight:
            if rows:
                out.append((x, off_x + 1, None, None, (pow(2, r, x) - 1) % x, True,
                            bool(_gmpy2.is_prime(x)) if _gmpy2 else None))
            else:
                out.append(x)
            continue
        if off_r and x % off_r == 2:
            if rows:
                m = (off_r * off_r - 2) % r
                col4 = x % m if m else None
                out.append((x, off_x + 1, x % off_r, col4, (pow(2, r, x) - 1) % x, True,
                            bool(_gmpy2.is_prime(x)) if _gmpy2 else None))
            else:
                out.append(x)
    return out


_QNR2_PRIMORIAL = {}
def _qnr2_primorial(bound):
    """Product of the odd primes up to `bound` (cached) -- for the gcd pre-screen."""
    P = _QNR2_PRIMORIAL.get(bound)
    if P is None:
        sieve = bytearray([1]) * (bound + 1)
        for i in range(2, int(bound ** 0.5) + 1):
            if sieve[i]:
                sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
        P = 1
        for p in range(3, bound + 1):
            if sieve[p]:
                P *= p
        _QNR2_PRIMORIAL[bound] = P
    return P


def qnr2_test(x, bound=2000):
    """Fast base-2 Euler-Jacobi probable-prime test -- the 'sports-car' of qnr2_filter.

    Same criterion as qnr2_filter's condition 1 (`offset % x == 2`, i.e. 2 is a
    quadratic NON-residue mod x, so x is a probable prime == 3 or 5 mod 8), but ~100x
    faster PER CANDIDATE via two tricks:
      * a small-prime gcd pre-screen (reject anything sharing a factor <= bound -- kills
        easy composites, and most Euler pseudoprimes, for pennies), and
      * `2**x == 2*(2**(x//2))**2`, so it needs ONE modexp instead of two.
    Benchmarks at ~gmpy2.is_prime speed (0.9-1.0x) from 100 to 3000 digits.

    HEURISTIC, not a proof.  A single-base test has infinitely many pseudoprimes;
    'smooth' composites whose prime factors ALL exceed `bound` (e.g. 3277 = 29*113 when
    bound<29) can still slip through.  Raise `bound` to catch more, or confirm survivors
    with gmpy2.is_prime / ECPP for a guarantee.  Returns True/False.
    """
    x = int(x)
    if x < 3 or x % 2 == 0:
        return False
    from math import gcd as _gcd
    g = _gcd(x, _qnr2_primorial(bound))
    if g != 1 and g != x:                       # a proper small factor -> composite
        return False
    r = x >> 1
    a = int(_gmpy2.powmod(2, r, x)) if _gmpy2 else pow(2, r, x)   # one modexp
    return (2 * a * a - 1 - a) % x == 2          # offset % x == 2, via 2^x = 2*(2^r)^2


# --- the walkable, indexable Lucas-Lehmer square-ladder ---------------------
# Not a shortcut -- an admiration piece.  The LL test is a PERFECT bijection mod a
# prime: every rung s -> s*s-2 (mod M) un-walks via sqrt_mod(s+2), because s+2 is
# always a perfect square (s_n + 2 = s_{n-1}^2).  So the whole chain is transparent
# and fully indexable both directions; its hardness is *distance*, not *disguise*,
# and the step where it lands on 0 (rung p-2) IS the primality certificate.

def ll_step(s, M):
    """One rung forward: the Lucas-Lehmer square-and-mod."""
    return (s * s - 2) % M

def ll_unstep(s, M, prev_hint=None):
    """One rung BACKWARD (only exact mod a prime): un-square s+2, which is always
    a perfect square. Returns the predecessor; with prev_hint picks the matching root."""
    import sympy as _sp
    roots = _sp.sqrt_mod((s + 2) % M, M, all_roots=True)
    if not roots:
        return None
    roots = [int(r) for r in roots]
    if prev_hint is not None and prev_hint in roots:
        return prev_hint
    return min(roots)

def ll_ladder(p, start=4, mode="normal", steps=None):
    """A DataFrame of the whole walkable ladder for M_p = 2**p - 1.

    mode -- which map to walk (see LL_MODES): 'normal' is real Lucas-Lehmer;
            the others are XOR/scale twists for indexing experiments.
    steps -- how many rungs (default p-2, the Lucas-Lehmer length).

    Columns:
      n         rung index (the step number == how far along the chain)
      s         the value at this rung
      bits      s in binary, p wide
      fwd       s -> next rung under this mode
      back      sqrt_mod(s+2) -- the previous rung recovered (ONLY meaningful for
                'normal'; NaN for other modes, whose inverse isn't a modular sqrt)
      back_ok   did back land on the true predecessor?  (normal mode only)
    For 'normal', the last rung is 0 exactly when M_p is prime, at n = p-2.
    """
    _need_pd()
    M = (1 << p) - 1
    f = LL_MODES[mode]
    nsteps = (p - 2) if steps is None else steps
    s = start
    seq = [s]
    for _ in range(nsteps):
        s = f(s, M)
        seq.append(s)
    rows = []
    for n, v in enumerate(seq):
        nxt = f(v, M)
        prev = seq[n - 1] if n > 0 else None
        back = ll_unstep(v, M, prev_hint=prev) if mode == "normal" else None
        rows.append({
            "n": n,
            "s": v,
            "bits": format(v, f"0{p}b"),
            "fwd": nxt,
            "back": back if back is not None else pd.NA,
            "back_ok": (n == 0) or (back == prev),
        })
    df = pd.DataFrame(rows).set_index("n")
    # keep big-int columns as clean Python ints (object dtype) -- no float/sci-notation,
    # and no precision loss past 2**53 for large p.
    for col in ("s", "fwd", "back"):
        df[col] = df[col].astype(object)
    df.attrs["M"] = M
    df.attrs["p"] = p
    df.attrs["prime"] = (seq[-1] == 0)
    df.attrs["count_to_zero"] = (p - 2) if seq[-1] == 0 else None
    df.attrs["fully_indexable"] = bool(df["back_ok"].all())
    if seq[-1] != 0:
        df.attrs["note"] = (f"p={p} is composite, so M_p is a composite Mersenne -- "
                            f"this ladder never reaches 0 (use ll_orbit to see the cycle).")
    return df


def ll_zero_row(p):
    """Summon the gem: for prime M_p, where the ladder touches 0 and what the last
    rung before it is.  The predecessor is always +-2**((p+1)/2) -- a pure power of 2
    (a shifted single bit) -- because it is a square root of 2 mod M_p."""
    M = (1 << p) - 1
    s = 4
    prev = None
    for n in range(1, p - 1):
        prev = s
        s = (s * s - 2) % M
        if s == 0:
            pw = pow(2, (p + 1) // 2, M)
            return {
                "M": M, "prime": True,
                "zero_at": n, "is_p_minus_2": (n == p - 2),
                "before_zero": prev,
                "power_of_two": prev in (pw, M - pw),
                "fwd_from_zero": (0 * 0 - 2) % M,   # == M-2
            }
    return {"M": M, "prime": False, "zero_at": None,
            "note": "never reaches 0 -- composite Mersenne"}


def ll_sqrt2(p):
    """The power-of-2 square root of 2 mod M_p:  2**((p+1)/2).  Exists for EVERY
    Mersenne (prime or composite) since 2**p == 1 mod (2**p - 1), so its square is
    2**(p+1) == 2.  Returns (pw, M - pw) -- the pure-power-of-2 root and its mate."""
    M = (1 << p) - 1
    pw = pow(2, (p + 1) // 2, M)
    return pw, M - pw


def ll_pow2(ps=(5, 7, 11, 13, 17, 19, 23, 29, 31), cap=5_000_000):
    """Do composite ladders share the primes' power-of-2 structure?

    For each p: the power-of-2 sqrt-of-2 (which EXISTS for every Mersenne), and
    whether the LL orbit (start 4) ever VISITS it.  Key fact: touching +-2**((p+1)/2)
    forces the very next rung to 0 (since (sqrt2)**2 - 2 == 0).  So an orbit can visit
    it ONLY if it zeroes out -- i.e. only primes.  Composites HAVE the gem but their
    orbit is forbidden from ever touching it."""
    _need_pd()
    import sympy as _sp
    rows = []
    for p in ps:
        M = (1 << p) - 1
        pw, mate = ll_sqrt2(p)
        # walk the orbit; note if we ever hit pw or its mate
        s, seen, visit = 4, {4}, None
        for n in range(1, cap):
            s = (s * s - 2) % M
            if s in (pw, mate):
                visit = n
            if s == 0 or s in seen:
                break
            seen.add(s)
        rows.append({
            "p": p, "prime": bool(_sp.isprime(M)),
            "sqrt2_pow2": pw, "is_pure_pow2": (pw & (pw - 1)) == 0,
            "orbit_visits_it": visit is not None,
            "visit_step": visit,
        })
    return pd.DataFrame(rows).set_index("p")


# --- LL map variants (for indexing experiments, not primality) --------------
# Each takes (s, M) -> next s.  'normal' is the real Lucas-Lehmer; the rest are
# Lars's XOR/scale twists, kept here so we can compare their ORBIT geometry
# (do they cycle cleanly = invertible/indexable, or merge = not).
LL_MODES = {
    "normal":       lambda s, M: (s * s - 2) % M,           # real Lucas-Lehmer
    "xor_insquare": lambda s, M: ((s * s ^ 5) - 2) % M,     # (s^2 XOR 5) - 2
    "xor_premod":   lambda s, M: ((s * s - 2) ^ 5) % M,     # (s^2 - 2) XOR 5
    "xor_times2":   lambda s, M: ((s * s ^ 5) * 2) % M,     # (s^2 XOR 5) * 2
    "sq_xor":       lambda s, M: (((s * s) ** 2) ^ 5) % M,  # (s^4) XOR 5
    "pow5":         lambda s, M: pow(s * s - 2, 5, M),      # (s^2 - 2) ^ 5  (5th power)
}

def ll_step_mode(s, M, mode="normal"):
    """One rung under any named LL variant (see LL_MODES)."""
    return LL_MODES[mode](s, M)

def ll_orbit_mode(M, mode="normal", start=4, cap=2_000_000):
    """Walk any LL variant on modulus M until it hits 0 or repeats.
    Returns kind / steps / tail / period / min_popcount / indexable.
      tail == 0  -> the orbit is a clean cycle from the start = INVERTIBLE (indexable)
      tail  > 0  -> a 'rho' with a merge point = some value has 2 preimages = NOT indexable
    """
    f = LL_MODES[mode]
    s = start
    seen = {s: 0}
    minpc = bin(s).count("1")
    for n in range(1, cap):
        s = f(s, M)
        pc = bin(s).count("1")
        if pc < minpc:
            minpc = pc
        if s == 0:
            return {"kind": "zero", "steps": n, "tail": None, "period": None,
                    "min_popcount": minpc, "indexable": None}
        if s in seen:
            tail = seen[s]
            return {"kind": "cycle", "steps": None, "tail": tail, "period": n - tail,
                    "min_popcount": minpc, "indexable": (tail == 0)}
        seen[s] = n
    return {"kind": "capped", "steps": None, "tail": None, "period": None,
            "min_popcount": minpc, "indexable": None}

def ll_modes(p, start=4, modes=None, cap=2_000_000):
    """Run ALL the LL variants on M_p side by side -- the indexing comparison.
    Shows, per mode: does it zero / cycle, the tail (0 = invertible), period,
    min popcount reached, and whether the orbit is indexable (clean cycle)."""
    _need_pd()
    M = (1 << p) - 1
    modes = modes or list(LL_MODES)
    rows = []
    for mode in modes:
        o = ll_orbit_mode(M, mode, start=start, cap=cap)
        rows.append({
            "mode": mode,
            "kind": o["kind"],
            "steps_to_0": o["steps"],
            "tail": o["tail"],
            "period": o["period"],
            "min_popcount": o["min_popcount"],
            "indexable": o["indexable"],
        })
    df = pd.DataFrame(rows).set_index("mode")
    df.attrs["M"] = M
    df.attrs["p"] = p
    return df


def ll_preimages(s, M):
    """The TWO predecessors of s under the LL map (the two square roots of s+2 mod M).
    Forward LL is 2-to-1, so backward it branches -- this is the indexing tree's fork.
    Returns a sorted list (1 or 2 roots; mod a prime a nonzero QR has exactly 2)."""
    import sympy as _sp
    return sorted(int(r) for r in _sp.sqrt_mod((s + 2) % M, M, all_roots=True))

def ll_back_tree(s, M, depth=4):
    """Walk the LL indexing tree BACKWARD from s as a binary tree of square-root
    choices.  Returns a DataFrame: each row is a value with the L/R path (address)
    that reaches it from s.  path '' is s itself; '0'/'1' pick the low/high root.
    Every value's forward LL walk climbs this address back up to s."""
    _need_pd()
    rows = [{"path": "", "value": int(s) % M, "popcount": bin(int(s) % M).count("1")}]
    frontier = [("", int(s) % M)]
    for _ in range(depth):
        nxt = []
        for path, v in frontier:
            pre = ll_preimages(v, M)
            for bit, r in zip("01", pre + pre[:1] * (2 - len(pre))):
                rows.append({"path": path + bit, "value": r,
                             "popcount": bin(r).count("1")})
                nxt.append((path + bit, r))
        frontier = nxt
    return pd.DataFrame(rows)

def ll_address(v, M, cap=100000):
    """The binary tree-address of value v under LL mod M: the L/R square-root path
    from the root 0 down to v.  Forward LL = drop the last bit (shift right); walking
    backward = append a bit.  Returns a bit-string, or None if v never reaches 0."""
    v = int(v) % M
    bits = []
    u = v
    seen = set()
    while u != 0:
        w = (u * u - 2) % M
        pre = ll_preimages(w, M)
        bits.append(str(pre.index(u)))
        u = w
        if u in seen or len(bits) > cap:
            return None
        seen.add(u)
    return "".join(reversed(bits))

def ll_from_address(addr, M):
    """Recover the value at a binary tree-address: walk backward from 0, taking the
    low root for '0' and the high root for '1' at each step.  Inverse of ll_address."""
    v = 0
    for b in addr:
        v = ll_preimages(v, M)[int(b)]
    return v

def ll_clean_landing(p):
    """For a prime M_p: does the rung-before-0 land on the CLEAN popcount-1 power of
    two (+2**((p+1)/2)) or on its high-popcount mate (M - 2**((p+1)/2))?  Both are
    square roots of 2; only the '+' one is a single bit.  Which side the ladder picks
    is an open little puzzle -- go hunt the rule.  Returns a dict (None if composite)."""
    M = (1 << p) - 1
    s = 4
    prev = None
    for n in range(1, p - 1):
        prev = s
        s = (s * s - 2) % M
        if s == 0:
            pw = pow(2, (p + 1) // 2, M)
            return {
                "p": p, "before_zero": prev,
                "popcount": bin(prev).count("1"),
                "clean_pow2": (prev == pw),         # True = the single-bit beam
                "landed_on": "+sqrt2 (clean pow2)" if prev == pw else "-sqrt2 (shadow)",
            }
    return None


def ll_clean_scan(exps=(3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279)):
    """Split the Mersenne PRIME exponents into 'clean' (land on the popcount-1 power
    of two before 0) vs 'shadow' (land on the high-popcount mate).  Feed it Mersenne
    prime exponents and hunt the pattern for which side wins."""
    clean, shadow = [], []
    for p in exps:
        r = ll_clean_landing(p)
        if r is None:
            continue
        (clean if r["clean_pow2"] else shadow).append(p)
    return {"clean": clean, "shadow": shadow}


def ll_orbit(M, start=4, cap=5_000_000):
    """Walk the LL map s -> s*s-2 (mod M) on ANY modulus until it hits 0 or cycles.

    Returns a dict:
      kind        'zero' (landed on 0) or 'cycle' (fell into a loop) or 'capped'
      steps       steps taken to reach 0 (if kind=='zero')
      tail        steps before the cycle starts (if kind=='cycle')
      period      cycle length (if kind=='cycle')
      hit_zero    convenience bool
    Prime Mersennes zero out at step p-2; composites loop forever.
    """
    s = start
    seen = {s: 0}
    for n in range(1, cap):
        s = (s * s - 2) % M
        if s == 0:
            return {"kind": "zero", "steps": n, "hit_zero": True,
                    "tail": None, "period": None}
        if s in seen:
            return {"kind": "cycle", "steps": None, "hit_zero": False,
                    "tail": seen[s], "period": n - seen[s]}
        seen[s] = n
    return {"kind": "capped", "steps": None, "hit_zero": False,
            "tail": None, "period": None}


def ll_compare(ps=(5, 7, 11, 13, 17, 19, 23, 29, 31)):
    """Prime vs composite Mersenne ladders side by side -- the fun comparison.
    Shows who zeroes (prime, at p-2) and who just loops forever (composite)."""
    _need_pd()
    import sympy as _sp
    rows = []
    for p in ps:
        M = (1 << p) - 1
        o = ll_orbit(M)
        rows.append({
            "p": p, "M": M, "prime": bool(_sp.isprime(M)),
            "kind": o["kind"],
            "zeroed_at": o["steps"], "expect_p-2": p - 2,
            "cycle_tail": o["tail"], "cycle_period": o["period"],
        })
    return pd.DataFrame(rows).set_index("p")


if __name__ == "__main__":
    if pd is not None:
        pd.set_option("display.width", 170)
        pd.set_option("display.max_columns", 20)
    print(decay(2048))
    print()
    print(where(2048))
