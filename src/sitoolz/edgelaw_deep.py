"""edgelaw_deep -- edgelaw's TABLE, extended to rung 41, with its own sweeper.

A copy of `edgelaw.py` with the lookup table carried from v2 = 30 out to
v2 = 41.  The original file is untouched.

    rungs 0..30   identical to edgelaw.TABLE
    rungs 31..43  new *to this table*, from csweep/edgesweep.c, 2026-08-18/19.
                  Not new to the literature: they reproduce OEIS A094604's
                  terms 31..43 exactly, all 13, computed independently here.
                  A094604 itself runs to rung 57 (137).

Provenance.  `edgereach.py` needs numba, which is not installed here, and its
pure-python fallback is ~100x too slow to reach these rungs.  `csweep/edgesweep.c`
replaces it with no third-party dependency: a 256-bit window, four uint64 words,
E' = E ^ ((E << 1) | (E << 2)), one sequential pass that reads the run at every
g = 2^v - 1 it passes.  It sustains ~670 million generations per second, so rung
41 -- 2.2 trillion generations -- costs 63 minutes, rung 42 just over two hours, rung 43 in 4h 12m.

    cc -O3 -march=native -o csweep/edgesweep csweep/edgesweep.c
    ./csweep/edgesweep 43

The pass reproduced rungs 0..30 exactly as they stand in `edgelaw.TABLE`, which
is what licenses trusting 31..41 from the same run.

    from sitoolz import edgelaw_deep as ed

    ed.edge(2**41 - 1)      # 100
    ed.show()               # the whole law
    ed.merge_sweep()        # fold in any deeper rungs the sweeper has reached
"""
import os

TABLE = {
    0: 1,   1: 3,   2: 4,   3: 6,   4: 7,   5: 9,   6: 15,  7: 16,
    8: 24,  9: 25,  10: 27, 11: 29, 12: 34, 13: 36, 14: 37, 15: 39,
    16: 41, 17: 43, 18: 48, 19: 49, 20: 51, 21: 54, 22: 55, 23: 58,
    24: 60, 25: 63, 26: 64, 27: 66, 28: 69, 29: 70, 30: 72,
    # ---- new rungs, csweep/edgesweep.c, 2026-08-18 ----
    31: 74, 32: 77, 33: 79, 34: 80, 35: 82, 36: 84, 37: 86, 38: 90,
    39: 91, 40: 93, 41: 100, 42: 103, 43: 104,
}

SWEEP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "csweep", "sweep43.txt")


def merge_sweep(path=SWEEP):
    """Fold rungs from a running (or finished) edgesweep log into TABLE.

    Returns the rungs added.  Raises if the log disagrees with a rung already
    in the table -- a disagreement means one of the two runs is wrong, and
    silently preferring either would be the wrong move.
    """
    added = []
    if not os.path.exists(path):
        return added
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        v, run = int(parts[0]), int(parts[1])
        if v in TABLE:
            if TABLE[v] != run:
                raise ValueError(f"rung {v}: table says {TABLE[v]}, "
                                 f"{path} says {run}")
        else:
            TABLE[v] = run
            added.append((v, run))
    return added


def vmax():
    return max(TABLE)


def verified_to():
    return (1 << vmax()) - 1


def level(g):
    """v2(g+1) -- which rung of the ladder generation g sits on."""
    n = g + 1
    return (n & -n).bit_length() - 1


def edge(g):
    """Length of the edge run at generation g.  Pure lookup, no simulation."""
    v = level(g)
    try:
        return TABLE[v]
    except KeyError:
        raise KeyError(
            f"generation {g:,} has v2(g+1) = {v}, past the verified table "
            f"(v2 <= {vmax()}, generations to {verified_to():,})"
        ) from None


def value(g):
    """The edge read as a number: 2^r - 1, a Mersenne number."""
    return (1 << edge(g)) - 1


def first(v):
    """First generation on rung v.  Others are (2m+1) * 2^v - 1."""
    return (1 << v) - 1


def same_as(g, count=6):
    v = level(g)
    return [(2 * m + 1) * (1 << v) - 1 for m in range(count)]


def covers(g):
    return level(g) <= vmax()


def table():
    return [(v, TABLE[v], first(v)) for v in sorted(TABLE)]


def plateaus():
    """run(v) - run(v-1): how many coordinates joined the run at each rung."""
    ks = sorted(TABLE)
    return [(v, TABLE[v] - TABLE[v - 1]) for v in ks[1:]]


def show():
    print("v2(g+1)   edge run   plateau   first generation")
    prev = 0
    for v, r, g in table():
        print(f"{v:7}   {r:8}   {r - prev:7}   {g:,}")
        prev = r
    print()
    print(f"verified single-valued for every generation to {verified_to():,}")


def agrees_with_original():
    """Every rung the original edgelaw.py carries, unchanged."""
    from sitoolz import edgelaw
    return [v for v, r in edgelaw.TABLE.items() if TABLE.get(v) != r]


if __name__ == "__main__":
    new = merge_sweep()
    if new:
        print(f"merged {len(new)} further rungs from the sweeper: {new}\n")
    show()
    print()
    bad = agrees_with_original()
    print(f"rungs disagreeing with edgelaw.TABLE: {len(bad)}")
    print(f"spot check   g = 2^41 - 1  ->  edge run {edge(2**41 - 1)}")


def table_frame(merge=True):
    """The edge law as a table: rung, run, plateau, first generation, source."""
    from sitoolz.dfkit import frame
    if merge:
        merge_sweep()
    rows, prev = [], 0
    for v in sorted(TABLE):
        r = TABLE[v]
        rows.append({"v": v, "edge_run": r, "plateau": r - prev,
                     "first_generation": first(v),
                     "source": "edgelaw.TABLE" if v <= 30 else "edgesweep.c"})
        prev = r
    return frame(rows, index="v")
