"""edgeframe -- the edge law as frames you can read down.

The period spectrum behind `edgelaw`, laid out row by row.  Every bit of
rule 30's light-cone edge is exactly periodic in the generation, with a
power-of-two period; the edge run at level v is the bits whose period
divides 2^v.  See EDGE_LAW.md.

    from sitoolz import edgeframe

    edgeframe.bits()          # one row per bit: its period, its plateau
    edgeframe.levels()        # one row per v2 rung: run, plateau, first gen
    edgeframe.ladder(1, 40)   # one row per generation: v2 and the edge run
    edgeframe.plateaus()      # just the flat stretches, longest first
    edgeframe.show()          # all of it, printed, without pandas

Frames are built from one evolution and cached; the first call pays for
it and the rest are free.  `depth()` rebuilds deeper.
"""
from collections import Counter

from sitoolz import edgelaw

try:
    import pandas as pd
except ImportError:                                   # pragma: no cover
    pd = None


def _need_pd():
    if pd is None:
        raise ImportError("edgeframe frames need pandas; use show() without it")


# ------------------------------------------------------------------ machine
_CACHE = {}


def _build(G, nbits, kmax, W=128):
    """Evolve once, extract every column, resolve every period.

    Bit j depends only on bits j-1 and j-2, so information flows strictly
    low-to-high and the W-bit window is exact for every bit it keeps.
    """
    m = (1 << W) - 1
    E = 1
    rows = [E]
    for _ in range(G - 1):
        E = (E ^ ((E << 1) | (E << 2))) & m
        rows.append(E)

    cols = [bytes((r >> b) & 1 for r in rows) for b in range(nbits)]

    P = []
    for C in cols:
        found = None
        for k in range(kmax + 1):
            p = 1 << k
            if C[p:] == C[:-p]:
                found = k
                break
        P.append(found)

    return {"G": G, "nbits": nbits, "kmax": kmax, "cols": cols, "P": P}


def depth(G=1 << 16, nbits=48, kmax=12):
    """Rebuild the cache at a given depth.  Larger G resolves more bits.

    G = 2^16, kmax=12   ~1 min, resolves ~34 bits      (default)
    G = 2^19, kmax=15   ~3 min, resolves ~39 bits
    """
    _CACHE["d"] = _build(G, nbits, kmax)
    return _CACHE["d"]


def _data():
    if "d" not in _CACHE:
        depth()
    return _CACHE["d"]


# ------------------------------------------------------------------- frames
def bits(nbits=None):
    """One row per bit of the edge.

      p          log2 of the bit's period in the generation
      period     the period itself, 2^p
      plateau    how many consecutive bits share this period
      first_on   least generation at which the bit is 1, = 2^p - 1
      in_phase   bit is 1 at every g == -1 mod 2^p        (claim C3)
      driver_odd its driver has odd weight over one period (claim C2)
      doubles    period is twice the bit below it -- C2 predicts this
    """
    _need_pd()
    d = _data()
    P, cols, G = d["P"], d["cols"], d["G"]
    res = [(b, p) for b, p in enumerate(P) if p is not None]
    if nbits is not None:
        res = res[:nbits]

    plateau_len = Counter(p for _, p in res)
    rows = []
    for b, p in res:
        step = 1 << p
        in_phase = all(cols[b][g] == 1 for g in range(step - 1, G, step))
        if b >= 2 and P[b - 1] is not None and P[b - 2] is not None:
            pd_exp = max(P[b - 1], P[b - 2])
            Pd = 1 << pd_exp
            driver = [x | y for x, y in zip(cols[b - 1][:Pd], cols[b - 2][:Pd])]
            odd = sum(driver) % 2 == 1
            doubles = p == pd_exp + 1
        else:
            odd = doubles = None
        rows.append({"bit": b, "p": p, "period": 1 << p,
                     "plateau": plateau_len[p], "first_on": (1 << p) - 1,
                     "in_phase": in_phase, "driver_odd": odd,
                     "doubles": doubles})
    df = pd.DataFrame(rows).set_index("bit")
    df.attrs["generations"] = G
    df.attrs["C2_agree"] = bool((df.driver_odd == df.doubles).where(
        df.driver_odd.notna(), True).all())
    df.attrs["C3_all_in_phase"] = bool(df.in_phase.all())
    return df


def levels(vmax=None):
    """One row per v2 rung -- the edge law itself, with its spectrum.

      run        edgelaw.TABLE[v], length of the edge run at this level
      plateau    run - previous run, = how many bits joined at this rung
      bit_lo     first bit that came in phase here
      bit_hi     last one
      first_gen  least generation on this rung, = 2^v - 1
      spectrum   #{b : period(b) divides 2^v} -- should equal run  (C5)
      digits     decimal digits of the edge read as 2^run - 1
    """
    _need_pd()
    d = _data()
    P = d["P"]
    T = edgelaw.TABLE
    vs = sorted(T)
    if vmax is not None:
        vs = [v for v in vs if v <= vmax]

    rows, prev = [], 0
    for v in vs:
        run = T[v]
        resolved = all(p is not None for p in P[:run])
        spec = (sum(1 for p in P if p is not None and p <= v)
                if resolved and v < d["kmax"] else None)
        rows.append({"v": v, "run": run, "plateau": run - prev,
                     "bit_lo": prev, "bit_hi": run - 1,
                     "first_gen": (1 << v) - 1, "spectrum": spec,
                     "agree": None if spec is None else spec == run,
                     "digits": len(str((1 << run) - 1))})
        prev = run
    df = pd.DataFrame(rows).set_index("v")
    df["spectrum"] = df.spectrum.astype("Int64")     # keep it an integer column
    checked = df.agree.dropna()
    df.attrs["C5_levels_checked"] = int(len(checked))
    df.attrs["C5_all_agree"] = bool(checked.all())
    return df


def ladder(g0=1, g1=40):
    """One row per generation -- what the law says, and what it measures.

      n          g + 1, the number whose valuation matters
      v2         v2(n), the rung
      run        the edge run, by lookup
      measured   the edge run, by evolving the edge
      same_rung  first other generation sharing this rung
    """
    _need_pd()
    rows = []
    for g in range(g0, g1 + 1):
        v = edgelaw.level(g)
        run = edgelaw.edge(g) if edgelaw.covers(g) else None
        rows.append({"g": g, "n": g + 1, "v2": v, "run": run,
                     "measured": edgelaw.measure(g),
                     "same_rung": edgelaw.same_as(g, 2)[1]})
    df = pd.DataFrame(rows).set_index("g")
    df.attrs["all_agree"] = bool((df.run == df.measured).all())
    return df


def plateaus(minlen=1):
    """The flat stretches of the spectrum, longest first.

    A plateau of length L at v means L consecutive bits share period 2^v
    and come in phase together, so the edge law jumps by L at that rung.
    """
    _need_pd()
    df = levels()
    out = (df[df.plateau >= minlen][["plateau", "bit_lo", "bit_hi", "run"]]
           .sort_values(["plateau", "bit_lo"], ascending=[False, True]))
    out.attrs["longest"] = int(df.plateau.max())
    out.attrs["distribution"] = dict(sorted(Counter(df.plateau).items()))
    return out


# -------------------------------------------------------------------- plain
def show(vmax=30):
    """The whole picture, printed.  No pandas needed."""
    T = edgelaw.TABLE
    print("  v   run   plateau   bits        first gen")
    prev = 0
    for v in sorted(T):
        if v > vmax:
            break
        run = T[v]
        print(f" {v:2}   {run:3}   {run - prev:7}   "
              f"{prev}..{run - 1:<6} {(1 << v) - 1:>12,}   "
              f"{'#' * (run - prev)}")
        prev = run
    dist = Counter(T[v] - (T[v - 1] if v else 0) for v in sorted(T))
    print()
    print("plateau lengths:", dict(sorted(dist.items())))
    print("longest        :", max(dist))


if __name__ == "__main__":
    show()
