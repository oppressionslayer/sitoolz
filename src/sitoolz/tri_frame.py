"""tri_frame -- dataframes over tri_xor's sixteen drivers.

`tri_xor.py` is deliberately standalone: "one file, stdlib only, no imports
from anywhere else."  That promise is worth keeping, so the frames live here
instead and import it rather than modifying it.

This is also where the 2026-08-23 shoulder analysis lives, because the two
questions turned out to be the same question:

    a SOLID RUN of r ones, read as an integer, is  2^r - 1   -- Mersenne-form
    a TWO-BIT state,       read as an integer, is 2^a + 2^b  -- Fermat-form

and the periodicity theorem says these events cannot be sporadic.  Every
coordinate is exactly periodic with a power-of-two period, so if a shoulder
state ever appears it must recur on a schedule the theory already constrains.
It does, and the schedule is the interesting part.

    from sitoolz import tri_frame as tf

    tf.survey_frame()          # all 16 drivers: periods, delta density, opacity
    tf.spectrum_frame(tx.OR)   # one driver, coordinate by coordinate
    tf.shoulder_frame()        # solid runs and two-bit states, all 16
    tf.event_frame(tx.A)       # WHEN the shoulder events happen -- the schedule
    tf.anf_frame()             # algebraic normal form vs shoulder capability
    tf.edge_frame()            # the v2 law, tri_xor's own table
    tf.report()                # the lot

Findings this module records, all reproducible from it:

  * driver `a` puts solid runs at n = 2^k - 1 and two-bit states at n = 2^k --
    the two shoulders of a power of two -- and the pattern is self-similar,
    repeating shifted by the window.  Gaps read 1,2,4,8,16 -- the ruler.
  * `OR`, rule 30's edge, produces SIX solid runs in 8192 generations and
    ZERO two-bit states.  Every other non-degenerate driver produces hundreds
    or thousands of both.  It is a measured outlier.
  * no predicate on the algebraic normal form separates the seven drivers that
    make sparse states from the nine that do not -- searched to three literals
    and all linear combinations.  `1+a` makes 128 of them and `1+b` makes none,
    so the property is DEPTH-SENSITIVE: it lives in how the driver couples to
    a = x_(j-1) and b = x_(j-2), not in its algebra.

Uses pandas when available, dfkit otherwise.
"""
from sitoolz import tri_xor as tx


def _frame(rows, index=None):
    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        return df.set_index(index) if index else df
    except ImportError:
        from sitoolz.dfkit import frame
        return frame(rows, index=index)


# ------------------------------------------------------------------- states
def states(driver=tx.OR, nbits=32, steps=1 << 13, seed=None):
    """The whole run as integers, one per generation.  Bit j is coordinate j."""
    cols = tx.evolve(driver, nbits, steps, seed)
    out = []
    for t in range(steps):
        v = 0
        for j in range(nbits):
            if cols[j][t]:
                v |= 1 << j
        out.append(v)
    return out


def is_solid(v):
    """v is 2^r - 1 for some r >= 1 -- a solid run of ones."""
    return bool(v) and (v & (v + 1)) == 0


def is_two_bit(v):
    """v is 2^a + 2^b -- the sparse shoulder."""
    return bin(v).count("1") == 2


# ------------------------------------------------------------------- frames
def survey_frame(nbits=40, steps=1 << 12, kmax=11):
    """tri_xor.survey as a frame: what resolves, max period, doubling density."""
    return _frame(tx.survey(nbits, steps, kmax), index="code")


def spectrum_frame(driver=tx.OR, nbits=40, steps=1 << 13, kmax=12):
    """One driver's spectrum, coordinate by coordinate.

      j          the coordinate
      p          period exponent, or None if unresolved in the window
      period     2^p
      doubled    True where the period doubled from j-1
    """
    P = tx.spectrum(driver, nbits, steps, kmax)
    rows = []
    for j, p in enumerate(P):
        d = None
        if j and p is not None and P[j - 1] is not None:
            d = p == P[j - 1] + 1
        rows.append({"j": j, "p": p, "period": None if p is None else 1 << p,
                     "doubled": d})
    return _frame(rows, index="j")


def shoulder_frame(nbits=32, steps=1 << 13):
    """The shoulder census: how often each driver makes each form.

      solid      generations whose state is 2^r - 1   (Mersenne-form)
      two_bit    generations whose state is 2^a + 2^b (Fermat-form)
      pairs      a solid run IMMEDIATELY followed by a two-bit state
      first_pair the first such transition, as "t: solid -> sparse"
    """
    rows = []
    for code in range(16):
        S = states(code, nbits, steps)
        solid = [(t, v) for t, v in enumerate(S) if is_solid(v)]
        two = [(t, v) for t, v in enumerate(S) if is_two_bit(v)]
        pairs = [(t, S[t], S[t + 1]) for t, _ in solid
                 if t + 1 < len(S) and is_two_bit(S[t + 1])]
        rows.append({"code": code, "driver": tx.NAMES[code],
                     "solid": len(solid), "two_bit": len(two),
                     "pairs": len(pairs),
                     "makes_sparse": len(two) > 0,
                     "first_pair": (f"t={pairs[0][0]}: {pairs[0][1]} -> {pairs[0][2]}"
                                    if pairs else None)})
    return _frame(rows, index="code")


def event_frame(driver=tx.A, nbits=32, steps=1 << 13, limit=24):
    """WHEN the shoulder events happen -- the schedule the theorem forces.

    For driver `a` the solid runs land on n = 2^k - 1 and the two-bit states on
    n = 2^k, and the whole pattern repeats shifted by the window.  That is the
    two shoulders of a power of two, on a self-similar schedule.
    """
    S = states(driver, nbits, steps)
    solid = [t for t, v in enumerate(S) if is_solid(v)][:limit]
    two = [t for t, v in enumerate(S) if is_two_bit(v)][:limit]
    n = max(len(solid), len(two))
    rows = []
    for i in range(n):
        st = solid[i] if i < len(solid) else None
        tw = two[i] if i < len(two) else None
        rows.append({"i": i, "solid_at": st,
                     "solid_gap": None if i == 0 or st is None else st - solid[i - 1],
                     "two_bit_at": tw,
                     "two_gap": None if i == 0 or tw is None else tw - two[i - 1]})
    return _frame(rows, index="i")


def anf(code):
    """(c0, c1, c2, c3) of G = c0 + c1*a + c2*b + c3*ab over GF(2)."""
    g = lambda a, b: (code >> ((a << 1) | b)) & 1
    c0 = g(0, 0)
    return (c0, g(1, 0) ^ c0, g(0, 1) ^ c0,
            g(0, 0) ^ g(0, 1) ^ g(1, 0) ^ g(1, 1))


def anf_text(code):
    c0, c1, c2, c3 = anf(code)
    t = [s for s, c in (("1", c0), ("a", c1), ("b", c2), ("ab", c3)) if c]
    return " + ".join(t) or "0"


def anf_frame(nbits=32, steps=1 << 13):
    """Algebraic normal form against shoulder capability.

    No predicate on (c0,c1,c2,c3) separates the seven drivers that make sparse
    states from the nine that do not -- searched to three literals and all
    linear combinations, 2026-08-23.  The reason is visible in the frame:
    `1 + a` makes them and `1 + b` does not, so the property depends on the
    SHIFT DEPTH the variable sits at, not on the Boolean function.
    """
    counts, rows = {}, []
    for code in range(16):
        S = states(code, nbits, steps)
        counts[code] = (sum(1 for v in S if is_solid(v)),
                        sum(1 for v in S if is_two_bit(v)))
    for code in range(16):
        c0, c1, c2, c3 = anf(code)
        solid, two = counts[code]
        rows.append({"code": code, "driver": tx.NAMES[code], "anf": anf_text(code),
                     "c0": c0, "c1": c1, "c2": c2, "c3": c3,
                     "degree": 2 if c3 else (1 if (c1 or c2) else 0),
                     "solid": solid, "two_bit": two,
                     "makes_sparse": two > 0})
    return _frame(rows, index="code")


def edge_frame(vmax=12, **kw):
    """tri_xor's own v2 law: run(g) = table[v2(g+1)], and the first generation."""
    tbl = tx.edge_table(vmax, **kw)
    rows, prev = [], 0
    for v, r in enumerate(tbl):
        rows.append({"v": v, "edge_run": r, "plateau": r - prev,
                     "first_gen": (1 << v) - 1,
                     "measured": tx.edge_measure((1 << v) - 1)})
        prev = r
    return _frame(rows, index="v")


def report():
    print("all sixteen drivers -- shoulder census")
    print(shoulder_frame()); print()
    print("driver `a` -- when the shoulders happen")
    print(event_frame(tx.A, limit=14)); print()
    print("algebraic normal form vs sparse capability")
    print(anf_frame())


if __name__ == "__main__":
    report()
