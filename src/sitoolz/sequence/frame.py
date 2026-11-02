"""The sequence as a DataFrame, plus the XOR-square selector used as an index.

Two pieces:

  sequence_frame()    the bounded-gap sequence from better.py, tabulated
  xor_square_index()  your loop: x a perfect square AND isqrt(x) ^ 5 a
                      perfect square.  ^ is XOR, so the real condition is
                      "n XOR 5 is square", reported at x = n^2.

The second is an index of a sort -- the surviving n are sparse and they land
on rows of the first, so you can use them to select.
"""

import math

import gmpy2
import pandas as pd
from sympy import isprime, primepi

from .better import futurex, odd_part, pf, step

XOR_K = 5


# ------------------------------------------------------------------ the terms
def sequence_frame(lo=2, hi=60):
    rows, gap = [], 0
    for x in range(lo, hi):
        W, nly = step(x, gap)
        gap -= nly
        odd = odd_part(W)
        rows.append({
            "x": x,
            "futurex": futurex(x),
            "pf": pf(x),
            "W": W,
            "nly": nly,
            "odd": odd,
            "v2": (W & -W).bit_length() - 1,
            "gap": gap,
            "pidx": int(primepi(odd)),          # A'-style: index of the odd part
            "xor5": x ^ XOR_K,
            "xor5_sq": bool(gmpy2.is_square(x ^ XOR_K)),
        })
    return pd.DataFrame(rows).set_index("x")


# ------------------------------------------------------------------ the index
def xor_square_index(limit=16430, k=XOR_K):
    """Your loop, tabulated. x must be square and isqrt(x) ^ k must be square."""
    rows = []
    for x in range(limit):
        if gmpy2.is_square(x) and gmpy2.is_square(math.isqrt(x) ^ k):
            n = math.isqrt(x)
            rows.append({
                "i": len(rows),
                "n": n,                          # the index itself
                "n_xor_k": n ^ k,
                "root": math.isqrt(n ^ k),       # its exact square root
                "x": x,                          # = n^2
                "xm2_xor_k": (x - 2) ^ k,
            })
    return pd.DataFrame(rows).set_index("i")


def xor_square_closed(m, k=XOR_K):
    """The same index without searching. Every m yields exactly one n.

    n = m^2 ^ k, and for k = 5 (binary 101) that XOR only touches bits 0 and 2,
    so it collapses to three cases -- m^2 mod 8 is only ever 0, 1 or 4:

        m odd        ->  n = m^2 + 3
        m = 0 mod 4  ->  n = m^2 + 5
        m = 2 mod 4  ->  n = m^2 - 3

    The loop over 16430 numbers finds nothing a formula doesn't give directly.
    """
    if k != XOR_K:
        return m * m ^ k
    return m * m + (3 if m % 2 else (5 if m % 4 == 0 else -3))


def xor_square_frame(limit=16430, k=XOR_K):
    """The index built from the closed form rather than by searching."""
    rows, m = [], 0
    while xor_square_closed(m, k) ** 2 < limit:
        n = xor_square_closed(m, k)
        rows.append({"m": m, "n": n, "n_xor_k": n ^ k, "root": m, "x": n * n,
                     "offset": n - m * m})
        m += 1
    return pd.DataFrame(rows).set_index("m")


if __name__ == "__main__":
    pd.set_option("display.width", 120)

    seq = sequence_frame()
    idx = xor_square_index()

    print("=== the sequence ===")
    print(seq.head(12).to_string())

    print(f"\n=== the XOR-{XOR_K} square index (n such that n ^ {XOR_K} is square) ===")
    print(idx.to_string())

    hits = idx[idx["n"].isin(seq.index)]["n"].tolist()
    print(f"\n=== index applied: seq.loc{hits} ===")
    print(seq.loc[hits].to_string())

    closed = xor_square_frame()
    print(f"\n=== same index, closed form (no search) ===")
    print(closed.to_string())
    print(f"\n  matches the brute-force loop: "
          f"{sorted(closed['n']) == sorted(idx['n'])}")

    print(f"\n  {len(idx)} survivors below {16430};  {len(hits)} land inside the "
          f"sequence's range x={seq.index.min()}..{seq.index.max()}")
    print(f"  all odd parts prime: {seq['odd'].map(isprime).all()}   "
          f"gap range {seq['gap'].min()} .. {seq['gap'].max()}")
