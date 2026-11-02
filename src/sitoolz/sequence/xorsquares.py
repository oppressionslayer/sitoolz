"""Your XOR-square scan, cleaned up so it pastes without breaking.

The original failed for two reasons, neither of them logic:
  1. a blank line inside the `for` body -- IPython ends the block there
  2. the long print() wrapped mid-token, splitting `gmpy2` from `.is_square`

Nothing about the math changed. Note `^` is XOR, not exponentiation.
"""

import math

import gmpy2

K = 5
LIMIT = 16430


def scan(limit=LIMIT, k=K, show=True):
    """x must be a perfect square AND isqrt(x) ^ k must be a perfect square."""
    hits = []
    for x in range(limit):
        n = math.isqrt(x)
        if not (gmpy2.is_square(n ^ k) and gmpy2.is_square(x)):
            continue
        hits.append(x)
        if show:
            print(f"{len(hits)-1:3}  x={x:6}  n={n:4}  "
                  f"n^{k}={n ^ k:5}  sqrt={math.isqrt(n ^ k):4}  "
                  f"(x-2)^{k}={(x - 2) ^ k:6}  sqrt(x)={math.isqrt(x):4}")
    return hits


def scan_closed(limit=LIMIT, k=K):
    """Same answer with no search. Every m yields exactly one n = m^2 ^ k.

    For k = 5 the XOR only touches bits 0 and 2, and m^2 mod 8 is only ever
    0, 1 or 4 -- so it collapses to three cases:

        m odd        ->  n = m^2 + 3
        m % 4 == 0   ->  n = m^2 + 5
        m % 4 == 2   ->  n = m^2 - 3
    """
    out, m = [], 0
    while (m * m ^ k) ** 2 < limit:
        out.append((m * m ^ k) ** 2)
        m += 1
    return sorted(out)


if __name__ == "__main__":
    hits = scan()
    print(f"\n{len(hits)} hits below {LIMIT}")
    print(f"closed form agrees: {sorted(hits) == scan_closed()}")
