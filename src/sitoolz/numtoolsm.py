"""numtoolsm - perfect machines. Import as ntm.

    n(i) = T(i^2)   for an involution T.

A machine never fails. T(T(n)) = n unconditionally, so:

    guarantee   T(n) is a perfect square, for EVERY i
    distinct    involutions are bijections, so no collisions, free
    invertible  rank(n) = isqrt(T(n)), you get it rather than design it

The only design rule found so far:

    LOCAL involutions  (act on fixed-size bit blocks, never look at the top
                        bit) have NO HORIZON -- valid at 10^50 as at 3.
    GLOBAL involutions (need to know the width) always have a horizon, and
                        they fail SILENTLY at it.

The `leak` column measures how crackable a machine is: the distance from the
output to the nearest perfect square. XOR-with-5 leaks <= 5 on every term,
which is why that one is findable. Nibble-reversal leaks enormously.

DESIGN RULE, found 2026-08-03 -- ROTATIONS LEAK, REVERSALS DO NOT.

    A rotation is a SHIFT whenever nothing wraps around the block edge.
    A shift is multiplication by a power of 2.
    Even powers of 2 are perfect squares.
    So a rotation maps some squares straight to other squares, and on those
    rows the disguise vanishes completely (leak = 0).

    Reversal has no k for which "reverse the bits" equals "multiply by k",
    so it never coincides with arithmetic and never leaks this way.

Measured over i = 1..3000, fraction of rows with leak = 0:

    swap_pairs  0.00%     xor255      0.00%     xor5        0.07%
    nibrev      0.20%     byterev     0.33%     bitrev*     0.33%
    bestlocal   0.53%     nibswap     2.07%  <- rotate bytes by 4
                          rot2        2.50%  <- rotate nibbles by 2

The two worst are the two rotations. rot2 IS multiplication by 4 whenever
every nibble is 0-3 (and division by 4 when every nibble is in {0,4,8,C}),
which is why its failures come in swapping pairs: 1<->4, 16<->64, 49<->196.

THE SCREEN. Before trusting a new involution, check whether T(x)/x is ever a
perfect-square rational:

    from fractions import Fraction
    bad = [i for i in range(1, 3000)
           if (lambda f: isqrt(f.numerator)**2 == f.numerator
                     and isqrt(f.denominator)**2 == f.denominator)
              (Fraction(T(i*i), i*i))]

Note this is an identity, not a prediction -- T(i^2)/i^2 = (p/q)^2 means
T(i^2) = (p*i/q)^2, a square, by definition. Useful as a fast check, not as
insight. The insight is which FAMILIES of involution have the property, and
the answer is: rotations do, reversals do not.

(A wrong guess worth recording: counting DISTINCT values of T(x)/x does not
predict leaking. swap_pairs has the fewest distinct ratios of any local
machine, 1479 of 2999, and leaks zero percent.)
"""

from math import isqrt

import pandas as pd


# ------------------------------------------------------------- involutions
def xor_k(k):
    f = lambda n: n ^ k
    f.__name__ = f"xor_{k}"
    return f, True, None                      # local, no horizon


def nibrev(n):
    """Reverse each 4-bit block in place. Local -> horizon-free."""
    out = sh = 0
    while n:
        b = n & 15
        out |= (((b & 1) << 3) | ((b & 2) << 1) | ((b & 4) >> 1) | ((b & 8) >> 3)) << sh
        n >>= 4
        sh += 4
    return out


def swap_pairs(n):
    """Swap bits (0,1), (2,3), ... Local -> horizon-free."""
    out = sh = 0
    while n:
        b = n & 3
        out |= (((b & 1) << 1) | ((b & 2) >> 1)) << sh
        n >>= 2
        sh += 2
    return out


def byterev(n):
    """Reverse each 8-bit block in place. Local -> horizon-free."""
    out = sh = 0
    while n:
        b, r = n & 255, 0
        for _ in range(8):
            r = (r << 1) | (b & 1)
            b >>= 1
        out |= r << sh
        n >>= 8
        sh += 8
    return out


def bitrev(width):
    """Reverse the low `width` bits. GLOBAL -> horizon at i = 2^(width/2)."""
    def f(n):
        r = 0
        m = n
        for _ in range(width):
            r = (r << 1) | (m & 1)
            m >>= 1
        return r
    f.__name__ = f"bitrev_{width}"
    return f, False, isqrt(1 << width)



def rot2(n):
    """Rotate each 4-bit block by 2. Twice = rotate by 4 = identity. Local."""
    out = sh = 0
    while n:
        b = n & 15
        out |= (((b << 2) | (b >> 2)) & 15) << sh
        n >>= 4
        sh += 4
    return out


def nibswap(n):
    """Swap the two nibbles of each byte. Local."""
    out = sh = 0
    while n:
        b = n & 255
        out |= (((b & 15) << 4) | (b >> 4)) << sh
        n >>= 8
        sh += 8
    return out



def bestlocal(n):
    """Best local involution found by exhaustive search over all 764
    involutive bit-permutations on 8-bit blocks: perm=(6,7,5,4,3,2,0,1).
    Median leak 1112 -- 13x better than byterev, 200x better than xor5,
    and still horizon-free.
    """
    P = (6, 7, 5, 4, 3, 2, 0, 1)
    out = sh = 0
    while n:
        b, w = n & 255, 0
        for i in range(8):
            if (b >> i) & 1:
                w |= 1 << P[i]
        out |= w << sh
        n >>= 8
        sh += 8
    return out


REGISTRY = {
    "xor5":       (lambda n: n ^ 5,   True,  None),
    "xor255":     (lambda n: n ^ 255, True,  None),
    "swap_pairs": (swap_pairs,        True,  None),
    "nibrev":     (nibrev,            True,  None),
    "rot2":       (rot2,              True,  None),
    "nibswap":    (nibswap,           True,  None),
    "byterev":    (byterev,           True,  None),
    "bestlocal":  (bestlocal,         True,  None),
    "bitrev24":   bitrev(24),
    "bitrev32":   bitrev(32),
}


# ----------------------------------------------------------------- machine
class Machine:
    def __init__(self, name, T=None, gen=None, ginv=None):
        """name from REGISTRY, or supply your own involution T.

        gen  : the indexed generator. defaults to i -> i*i.
        ginv : its inverse, needed for rank() if gen is not squaring.
               e.g. gen=lambda i: i**3, ginv=lambda v: round(v ** (1/3))
        """
        if T is None:
            T, self.local, self.horizon = REGISTRY[name]
        else:
            self.local, self.horizon = True, None
        self.name, self.T = name, T
        self.gen = gen or (lambda i: i * i)
        self.ginv = ginv or (isqrt if gen is None else None)

    def unrank(self, i):
        return self.T(self.gen(i))

    def rank(self, n):
        if self.ginv is None:
            raise ValueError(f"{self.name}: pass ginv= to invert a custom gen")
        return self.ginv(self.T(n))

    def holds(self, i):
        back = self.T(self.unrank(i))
        return back == self.gen(i) and isqrt(back) ** 2 == back

    def take(self, n, start=1):
        """The first n outputs, as a plain list."""
        return [self.unrank(i) for i in range(start, start + n)]

    def frame_at(self, indices):
        """Same columns as frame(), but for any list of indices you like --
        so you can look at i = 10**12 without generating a trillion rows."""
        rows, prev = [], None
        for i in indices:
            g = self.gen(i)
            n = self.T(g)
            back = self.T(n)
            r = isqrt(n)
            rows.append({
                "i": i, "gen": g, "n": n,
                "ok": back == g and isqrt(back) ** 2 == back,
                "rank": isqrt(back), "rank_ok": isqrt(back) == i,
                "digits": len(str(n)), "bits": n.bit_length(),
                "popcount": bin(n).count("1"),
                "leak": min(n - r * r, (r + 1) ** 2 - n),
                "rising": None if prev is None else n > prev,
            })
            prev = n
        return pd.DataFrame(rows).set_index("i")

    def frame(self, lo=1, hi=40):
        """The machine as a DataFrame. Slice it however you like."""
        rows, prev = [], None
        for i in range(lo, hi):
            g = self.gen(i)
            n = self.T(g)
            back = self.T(n)
            r = isqrt(n)
            leak = min(n - r * r, (r + 1) ** 2 - n)      # distance to nearest square
            rows.append({
                "i": i, "gen": g, "n": n, "back": back,
                "ok": back == g and isqrt(back) ** 2 == back,
                "rank": isqrt(back), "rank_ok": isqrt(back) == i,
                "bits": n.bit_length(), "popcount": bin(n).count("1"),
                "leak": leak, "rising": None if prev is None else n > prev,
            })
            prev = n
        return pd.DataFrame(rows).set_index("i")

    def __repr__(self):
        h = "no horizon" if self.horizon is None else f"horizon i={self.horizon}"
        return (f"<Machine {self.name!r} [{'local' if self.local else 'GLOBAL'}, "
                f"{h}] {[self.unrank(i) for i in range(1, 6)]}>")


def compare(names=None, lo=1, hi=200):
    """Side-by-side: which machines hold, and how well they hide."""
    names = names or list(REGISTRY)
    rows = []
    for nm in names:
        m = Machine(nm)
        f = m.frame(lo, hi)
        rows.append({
            "machine": nm,
            "local": m.local,
            "horizon": m.horizon,
            "all_ok": bool(f.ok.all()),
            "rank_ok": bool(f.rank_ok.all()),
            "distinct": f.n.nunique() == len(f),
            "monotonic": bool(f.rising.dropna().all()),
            "leak_max": int(f.leak.max()),
            "leak_median": float(f.leak.median()),
        })
    return pd.DataFrame(rows).set_index("machine")


if __name__ == "__main__":
    pd.set_option("display.width", 150)
    print("=== side by side, i = 1..200 ===")
    print(compare().to_string(), "\n")
    print("=== nibrev, first 12 ===")
    print(Machine("nibrev").frame(1, 13).to_string(), "\n")
    print("=== xor5 leaks tiny (that is why it is crackable) ===")
    print(Machine("xor5").frame(1, 8)[["n", "leak", "popcount"]].to_string())
