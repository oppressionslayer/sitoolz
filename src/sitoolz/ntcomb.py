"""ntcomb -- combinations addressed by bits.

The same claim as the square towers, on a different set: an exact
bidirectional index beats a search.  C(m, j) was always just "how many
m-bit numbers have j ones", so a popcount census rebuilds Pascal's row.
Run it the other way and C(n, j) is odd exactly when j is a submask of n,
which is what draws Sierpinski -- no subdivision rule anywhere.  Gosper's
hack walks the j-subsets in order in four operations each, and the
combinatorial number system skips the walk entirely and addresses the i-th
subset straight from binomials.

    from sitoolz import ntcomb

    ntcomb.census(6)          # Pascal row 6, rebuilt by counting popcounts
    ntcomb.sierpinski(64)     # C(n,j) mod 2 -- the triangle, from submasks
    ntcomb.walk(6, 3)         # Gosper's hack, the masks in order
    ntcomb.unrank(7, 3)       # the 7th 3-subset, no enumeration
    ntcomb.verify(6, 3)       # walk == jump, and rank(unrank(i)) == i
"""
from collections import Counter
from math import comb

try:
    import pandas as pd
except ImportError:
    pd = None


def _need_pd():
    if pd is None:
        raise ImportError("pandas required")


def popcount(n):
    """Set bits in n, by Kernighan's method -- one iteration per set bit."""
    c = 0
    while n:
        n &= n - 1
        c += 1
    return c


# =====================================================================
# machine 03 -- Pascal, addressed by bits
# =====================================================================
def census(m=6):
    """Pascal's row m, rebuilt by counting popcounts of 0 .. 2^m - 1.

    C(m, j) is not a formula that happens to count subsets -- it *is* the
    number of m-bit patterns with j ones.  Counting them recovers the row
    exactly, which is the whole content of the binomial coefficient.

      j        ones in the pattern
      counted  how many of 0 .. 2^m-1 have that popcount
      binomial C(m, j)
      agree    the two are equal
    """
    _need_pd()
    h = Counter(popcount(n) for n in range(1 << m))
    rows = [{"j": j, "counted": h.get(j, 0), "binomial": comb(m, j),
             "agree": h.get(j, 0) == comb(m, j)} for j in range(m + 1)]
    df = pd.DataFrame(rows).set_index("j")
    df.attrs["m"] = m
    df.attrs["total"] = 1 << m
    return df


def submasks(n):
    """Every j with j & n == j, descending.  The standard subset trick."""
    out, j = [], n
    while True:
        out.append(j)
        if j == 0:
            return out
        j = (j - 1) & n


def parity_row(n):
    """C(n, j) mod 2 for j = 0 .. n, and whether j is a submask of n.

    Kummer's theorem in its cheapest form: C(n, j) is odd exactly when j's
    bits are a subset of n's.  No factorials, no division -- one AND.
    """
    _need_pd()
    rows = [{"j": j, "binomial_mod2": comb(n, j) % 2,
             "j_submask_of_n": (j & n) == j} for j in range(n + 1)]
    df = pd.DataFrame(rows).set_index("j")
    df.attrs["n"] = n
    df.attrs["n_bin"] = bin(n)[2:]
    df.attrs["all_agree"] = bool((df.binomial_mod2 == df.j_submask_of_n).all())
    return df


def sierpinski(depth=32, ascii_art=False):
    """C(n, j) mod 2 as a triangle -- Sierpinski, drawn by an AND.

    Nothing here is recursive and nothing subdivides.  Each cell is a
    single `(j & n) == j` test, and the gasket is what that predicate
    looks like when you plot it.

    Returns a 0/1 frame, or the picture as a string with `ascii_art`.
    """
    _need_pd()
    grid = [[1 if (j & n) == j else 0 for j in range(depth)]
            for n in range(depth)]
    if ascii_art:
        return "\n".join("".join("#" if c else "." for c in row[:n + 1])
                         for n, row in enumerate(grid))
    df = pd.DataFrame(grid)
    df.index.name = "n"
    df.columns.name = "j"
    df.attrs["ones"] = sum(map(sum, grid))
    return df


# =====================================================================
# machine 04 -- walk it, or jump straight to it
# =====================================================================
def next_mask(x):
    """Gosper's hack: the next integer with the same popcount, in four ops.

    Add the lowest set bit to carry the block up, then put the bits that
    fell off back at the bottom.  No searching, no filtering -- the masks
    come out in ascending order because the arithmetic keeps them there.
    """
    if x == 0:
        return 0
    c = x & -x
    r = x + c
    return (((r ^ x) >> 2) // c) | r


def walk(m=6, j=3, limit=None):
    """Every m-bit mask with j ones, in ascending order, via `next_mask`.

      mask    the integer
      bits    its binary, m wide
      subset  the set bit positions
      i       position in the walk -- and, as `verify` shows, the colex rank
    """
    _need_pd()
    if j == 0:
        vals = [0]
    else:
        x, vals = (1 << j) - 1, []
        while x < (1 << m):
            vals.append(x)
            if limit and len(vals) >= limit:
                break
            x = next_mask(x)
    rows = [{"i": i, "mask": x, "bits": format(x, f"0{m}b"),
             "subset": [t for t in range(m) if (x >> t) & 1]}
            for i, x in enumerate(vals)]
    df = pd.DataFrame(rows).set_index("i")
    df.attrs["m"], df.attrs["j"] = m, j
    df.attrs["expected"] = comb(m, j)
    df.attrs["complete"] = len(vals) == comb(m, j)
    return df


def unrank(i, j):
    """The i-th j-subset in colex order, straight from binomials.

    Greedy from the top: the largest c with C(c, j) <= i is the biggest
    element, subtract and repeat.  j binomial lookups, no enumeration --
    C(52, 5) is 2.6 million poker hands and this reaches any one of them
    in five steps.
    """
    out, rem = [], i
    for k in range(j, 0, -1):
        c = k - 1
        while comb(c + 1, k) <= rem:
            c += 1
        out.append(c)
        rem -= comb(c, k)
    return sorted(out)


def rank(subset):
    """Colex rank of a subset -- the exact inverse of `unrank`."""
    return sum(comb(c, k) for k, c in enumerate(sorted(subset), 1))


def jump(lo=0, hi=12, j=3):
    """`unrank` over a range, with the round trip checked on every row.

      subset  the i-th j-subset
      mask    the same thing as an integer
      rank    rank(subset) -- must come back as i
      ok      it did
    """
    _need_pd()
    rows = []
    for i in range(lo, hi):
        s = unrank(i, j)
        rows.append({"i": i, "subset": s,
                     "mask": sum(1 << c for c in s), "rank": rank(s)})
    df = pd.DataFrame(rows).set_index("i")
    df["ok"] = df.index == df["rank"]
    df.attrs["j"] = j
    df.attrs["all_ok"] = bool(df.ok.all())
    return df


def verify(m=6, j=3):
    """Walk and jump agree entry for entry, so each one checks the other.

    Ascending mask order and colex subset order are the same order -- that
    is the fact that makes the two halves of this module a pair rather
    than two separate tricks.

      walk_len   masks Gosper produced
      binomial   C(m, j)
      same_order unrank(i) == the walk's i-th subset, for every i
      round_trip rank(unrank(i)) == i, for every i
    """
    _need_pd()
    w = walk(m, j)
    subs = list(w.subset)
    same = all(unrank(i, j) == s for i, s in enumerate(subs))
    trip = all(rank(unrank(i, j)) == i for i in range(len(subs)))
    return pd.DataFrame([{
        "m": m, "j": j,
        "walk_len": len(subs),
        "binomial": comb(m, j),
        "complete": len(subs) == comb(m, j),
        "same_order": same,
        "round_trip": trip,
    }]).set_index(["m", "j"])
