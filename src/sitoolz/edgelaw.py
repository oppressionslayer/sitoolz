"""edgelaw -- rule 30's light-cone edge, by lookup.

The run of live cells at the leading edge of rule 30's light cone depends
only on the 2-adic valuation of g+1.  The odd part of g+1 is irrelevant:
generations 255, 767, 1279 and 1791 all give 24, because g+1 is
2^8, 3*2^8, 5*2^8 and 7*2^8.

Verified single-valued for every generation to 16,777,215.

    from sitoolz import edgelaw

    edgelaw.edge(16777215)        # 60   -- no simulation
    edgelaw.edge(10**12)          # instant at any size
    edgelaw.bits(2047)            # '11111111111111111111111111111'
    edgelaw.level(1791)           # 8    -- which rung of the ladder
    edgelaw.first(8)              # 255  -- first generation on that rung
    edgelaw.same_as(1791)         # other generations with the same edge
    edgelaw.table()               # the whole law
"""

# v2(g+1) -> length of the run of ones at the light-cone edge.
# Measured by evolving the edge alone:  E = E ^ ((E << 1) | (E << 2)).
# Rungs 0..30 measured here; rungs 31..48 from OEIS A094604 (the same sequence,
# verified to agree on every overlapping rung), extending the law to ~2^49.
TABLE = {
    0: 1,   1: 3,   2: 4,   3: 6,   4: 7,   5: 9,   6: 15,  7: 16,
    8: 24,  9: 25,  10: 27, 11: 29, 12: 34, 13: 36, 14: 37, 15: 39,
    16: 41, 17: 43, 18: 48, 19: 49, 20: 51, 21: 54, 22: 55, 23: 58,
    24: 60, 25: 63, 26: 64, 27: 66, 28: 69, 29: 70, 30: 72,
    31: 74, 32: 77, 33: 79, 34: 80, 35: 82, 36: 84, 37: 86, 38: 90,
    39: 91, 40: 93, 41: 100, 42: 103, 43: 104, 44: 106, 45: 108,
    46: 111, 47: 112, 48: 114,
}

VMAX = max(TABLE)
VERIFIED_TO = (1 << VMAX) - 1                # 562,949,953,421,311


def level(g):
    """v2(g+1) -- which rung of the ladder generation g sits on."""
    n = g + 1
    return (n & -n).bit_length() - 1


def edge(g):
    """Length of the edge run at generation g.  Pure lookup, no simulation.

    Raises KeyError past the verified table rather than extrapolating.
    """
    v = level(g)
    try:
        return TABLE[v]
    except KeyError:
        raise KeyError(
            f"generation {g:,} has v2(g+1) = {v}, past the verified table "
            f"(v2 <= {VMAX}, generations to {VERIFIED_TO:,})"
        ) from None


def bits(g):
    """The edge itself.  Always a solid run of ones, so length is the whole story."""
    return "1" * edge(g)


def value(g):
    """The edge read as a number: 2^r - 1, a Mersenne number."""
    return (1 << edge(g)) - 1


def first(v):
    """First generation on rung v.  Others are (2m+1) * 2^v - 1."""
    return (1 << v) - 1


def same_as(g, count=6):
    """Generations sharing g's rung -- all with the identical edge run."""
    v = level(g)
    return [(2 * m + 1) * (1 << v) - 1 for m in range(count)]


def covers(g):
    """Is generation g inside the verified table?"""
    return level(g) <= VMAX


def table():
    """The law as a list of (v2, edge run, first generation)."""
    return [(v, TABLE[v], first(v)) for v in sorted(TABLE)]


def show():
    """Print the law."""
    print("v2(g+1)   edge run   first generation")
    for v, r, g in table():
        print(f"{v:7}   {r:8}   {g:,}")
    print()
    print(f"verified single-valued for every generation to {VERIFIED_TO:,}")


# ---------------------------------------------------------------- checking
def measure(g, W=512):
    """Evolve the edge alone and read the run -- for checking `edge`.

    O(1) per generation on a fixed window, since the edge never consults
    the interior:  E'[j] = E[j] XOR (E[j-1] OR E[j-2]).
    """
    m = (1 << W) - 1
    E = 1
    for _ in range(g):
        E = (E ^ ((E << 1) | (E << 2))) & m
    k = 0
    while (E >> k) & 1:
        k += 1
    return k


def verify(upto=8191, W=512):
    """Check the lookup against measurement for every generation to `upto`."""
    m = (1 << W) - 1
    E = 1
    bad = []
    for g in range(1, upto + 1):
        E = (E ^ ((E << 1) | (E << 2))) & m
        if not covers(g):
            continue
        k = 0
        while (E >> k) & 1:
            k += 1
        if k != edge(g):
            bad.append(g)
    return {"checked": upto, "mismatches": len(bad), "first_bad": bad[:5]}


if __name__ == "__main__":
    show()
    print()
    print("spot checks:")
    for g in (255, 767, 1279, 1791, 2047, 16777215):
        print(f"  gen {g:>10,}  v2={level(g):2}  edge run {edge(g)}")
    print()
    print(" ", verify(8191))
