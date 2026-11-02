"""numtoolsd - XOR/Diophantine tools. Import as ntd.

    from sitoolz import numtools  as nt  # the main library
    from sitoolz import numtoolz  as ntz  # the disjoint pylib half
    from sitoolz import numtoolsd as ntd  # this one

Built from an evening of poking at `n XOR k is a square`. The organizing
idea that came out of it:

    ONE unknown   -> a construction wearing a search's clothes.
                     n ^ k = m^2  solves to  n = m^2 ^ k.  Bijection.
                     Infinitely many hits, closed form, quasi-polynomial.

    TWO unknowns  -> a real Diophantine condition.
                     m^2 ^ m = square  has m on both sides. Nothing to
                     solve for. Usually almost no solutions.

Everything here is either exploring one side of that line or measuring
where the line is.

Key facts the code encodes:

  * XOR is nim-addition / addition in GF(2)^n / carry-less addition.
  * XOR equals ordinary addition exactly when no carries occur, i.e. when
    a & b == 0.  By KUMMER'S THEOREM that happens iff C(a+b, b) is odd;
    equivalently by LUCAS, iff b's bits are a subset of (a+b)'s bits.
  * XOR with a FIXED k touches only bounded low bits, and m^2 mod 2^j is
    periodic in m -- so m^2 ^ k is always m^2 plus a periodic correction.
    A quasi-polynomial. That is why such sequences are not interesting
    enough for OEIS: nothing is left to discover.
  * XOR with a GROWING partner reaches unbounded bits and the periodicity
    dies -- but if the partner is small relative to m^2, the reachable
    window is too tight to contain any square but the obvious ones.
"""

from math import comb, isqrt

import pandas as pd

__all__ = [
    "is_square", "square_index",
    "xor_k_closed", "xor_k_squares", "xor_k_period", "xor_k_survey",
    "carryfree", "kummer_carries", "binom_odd", "kummer_check",
    "xor_window", "xor_self_squares", "near_square_xor",
    "consecutive_xor_squares",
]


# ----------------------------------------------------------------- basics
def is_square(v):
    return v >= 0 and isqrt(v) ** 2 == v


def square_index(n):
    """Where n sits relative to the squares.

    Bands between consecutive squares are always 2k+1 wide -- the gnomon,
    since n^2 = 1 + 3 + 5 + ... + (2n-1).
    """
    k = isqrt(n)
    lo, hi = k * k, (k + 1) ** 2
    return {"n": n, "k": k, "lo": lo, "hi": hi, "offset": n - lo,
            "width": hi - lo, "gnomon": 2 * k + 1, "is_square": lo == n,
            "nearest": lo if n - lo <= hi - n else hi}


# ------------------------------------------- one unknown: the tame family
def xor_k_closed(m, k=5):
    """n = m^2 ^ k. Every m gives exactly one n -- no search, no gaps."""
    return m * m ^ k


COLUMNS = ["m", "n", "n_xor_k", "root", "x", "offset", "carryfree"]


def xor_k_squares(k=5, limit=16430, squares_only=True, count=None):
    """n such that n ^ k is a perfect square. Returns a DataFrame.

    count : produce exactly this many rows, m = 0 .. count-1. Ignores limit.
            This is usually what you want -- `limit` bounds x, and since
            x = n^2 you would need limit=10**8 just to get 100 rows.
    limit : if count is None, include every term whose x is below this.
    squares_only : also require n itself to be a square (the original
            loop's condition). When False, x is n rather than n^2.

    Careful: m -> n is NOT monotonic for small m (n = 5, 4, 1 at m = 0, 1, 2),
    so a limit scan must not stop at the first term over the bound. It walks
    a computed ceiling and filters instead.
    """
    if count is not None:
        ms = range(max(count, 0))
    else:
        # n >= m^2 - k, and we need n <= isqrt(limit) (or limit), so
        # m^2 <= that bound + k. A small margin covers the XOR wobble.
        bound = isqrt(limit) if squares_only else limit
        ms = range(isqrt(max(bound + k, 0)) + 2)

    rows = []
    for m in ms:
        n = xor_k_closed(m, k)
        x = n * n if squares_only else n
        if count is None and x >= limit:
            continue
        rows.append({"m": m, "n": n, "n_xor_k": n ^ k, "root": isqrt(n ^ k),
                     "x": x, "offset": n - m * m,
                     "carryfree": carryfree(m * m, k)})

    if not rows:
        return pd.DataFrame(columns=COLUMNS).set_index("m")
    return pd.DataFrame(rows).set_index("m").sort_values("n")


def xor_k_period(k=5, probe=400):
    """The quasi-polynomial period of m^2 ^ k - m^2, and its offsets.

    Returns (period, offsets) or (None, ...) if none found -- a None here
    means the partner is effectively growing and the family is not tame.
    """
    off = [(m * m ^ k) - m * m for m in range(probe * 3)]
    for p in range(1, probe):
        if all(off[i] == off[i + p] for i in range(probe * 2 - p)):
            return p, off[:p]
    return None, off[:16]


def xor_k_survey(ks=range(1, 33)):
    """Which constants collapse cleanly? Small period == clean collapse."""
    rows = []
    for k in ks:
        p, off = xor_k_period(k)
        rows.append({"k": k, "bits": bin(k)[2:], "popcount": bin(k).count("1"),
                     "period": p, "offsets": off if p else None})
    return pd.DataFrame(rows).set_index("k")


# ------------------------------------------------------ carries / Kummer
def carryfree(a, b):
    """True iff adding a and b in base 2 produces no carries."""
    return (a & b) == 0


def kummer_carries(a, b):
    """Number of base-2 carries when adding a and b.

    Kummer: this equals the exponent of 2 dividing C(a+b, b).
    """
    carries, c = 0, 0
    while a or b or c:
        s = (a & 1) + (b & 1) + c
        c = s >> 1
        carries += c
        a >>= 1
        b >>= 1
    return carries


def binom_odd(a, b):
    """C(a+b, b) odd? Equivalent to carryfree(a, b) by Kummer/Lucas."""
    return comb(a + b, b) % 2 == 1


def kummer_check(a, b):
    """All four views of the same fact, for one pair."""
    return {"a": a, "b": b, "a_xor_b": a ^ b, "a_plus_b": a + b,
            "xor_eq_sum": (a ^ b) == a + b,
            "carryfree": carryfree(a, b),
            "carries": kummer_carries(a, b),
            "binom_odd": binom_odd(a, b)}


# --------------------------------------- two unknowns: the Diophantine side
def xor_window(m, partner=None):
    """How far can XOR move m^2, and which squares are reachable?

    XOR with p only flips bits inside p's width, so it moves you by less
    than 2^bit_length(p).  When the partner is small next to m^2 the window
    is too tight to hold any square but m^2 itself and (m-1)^2 on the edge
    -- which is exactly why `m^2 ^ m is a square` has only trivial answers.
    """
    p = m if partner is None else partner
    reach = 1 << p.bit_length()
    lo, hi = m * m - reach, m * m + reach
    inside = [j * j for j in range(max(isqrt(max(lo, 0)) - 1, 0), isqrt(hi) + 2)
              if lo < j * j < hi]
    return {"m": m, "m2": m * m, "partner": p, "reach": reach,
            "window": (lo, hi), "squares_in_window": inside,
            "actual": m * m ^ p, "hit": is_square(m * m ^ p)}


def xor_self_squares(limit=200000, partner=lambda m: m):
    """m such that m^2 ^ partner(m) is a perfect square.

    partner defaults to m itself -- the wild version, which has m on both
    sides and therefore almost no solutions. Pass any callable to explore.
    """
    return [m for m in range(limit) if is_square(m * m ^ partner(m))]


def near_square_xor(limit=20000, partner=lambda m: m, tol=8):
    """Widen the window: m^2 ^ partner(m) within tol of a square.

    The strict version is nearly empty, so this is where the life is.
    """
    rows = []
    for m in range(limit):
        v = m * m ^ partner(m)
        r = isqrt(v)
        d = min(v - r * r, (r + 1) ** 2 - v, key=abs)
        near = r * r if v - r * r <= (r + 1) ** 2 - v else (r + 1) ** 2
        if abs(d) <= tol:
            rows.append({"m": m, "m2": m * m, "value": v, "nearest_sq": near,
                         "root": isqrt(near), "delta": v - near})
    return pd.DataFrame(rows).set_index("m")


def consecutive_xor_squares(limit=300000, step=1):
    """k such that k^2 ^ (k+step)^2 is a square.

    The A221643 family (step=1). Both operands are LARGE here, so the
    window argument does not apply, carries reach high bits, and you get
    genuinely ragged output -- a real Diophantine condition.
    """
    return [k for k in range(1, limit) if is_square(k * k ^ (k + step) ** 2)]


# ---------------------------------------------------- ranking / unranking
#
# rank(v)   : value -> its index      ("which one is this?")
# unrank(i) : index -> its value      ("what is number i?")
#
# A sequence with a closed form does both in O(1). A sequence that is a real
# Diophantine condition can do NEITHER without walking -- that is the same
# tame/wild line as everywhere else in this module, seen from the other side.


class Indexer:
    """A sequence you can query in both directions.

    closed=True  -> unrank/rank are O(1) formulas
    closed=False -> both walk a membership test, with a cache
    """

    def __init__(self, name, unrank=None, rank=None, member=None, start=0):
        self.name, self._u, self._r, self._m = name, unrank, rank, member
        self.start, self.closed = start, unrank is not None
        self._cache = []

    def unrank(self, i):
        """The i-th term."""
        if self._u is not None:
            return self._u(i)
        self._grow(i + 1)
        return self._cache[i] if i < len(self._cache) else None

    def rank(self, v):
        """The index of v, or None if v is not a member."""
        if self._r is not None:
            return self._r(v)
        self._grow_to(v)
        return self._cache.index(v) if v in self._cache else None

    def take(self, n):
        return [self.unrank(i) for i in range(n)]

    def _grow(self, n):
        m = self._cache[-1] + 1 if self._cache else self.start
        while len(self._cache) < n:
            if self._m(m):
                self._cache.append(m)
            m += 1

    def _grow_to(self, v):
        m = self._cache[-1] + 1 if self._cache else self.start
        while m <= v:
            if self._m(m):
                self._cache.append(m)
            m += 1

    def __getitem__(self, i):
        return self.unrank(i)

    def __repr__(self):
        kind = "closed form, O(1)" if self.closed else "search, cached"
        return f"<Indexer {self.name!r} [{kind}] {self.take(6)}...>"


#: the squares themselves -- the simplest possible indexer
SQUARES = Indexer("squares", unrank=lambda i: i * i,
                  rank=lambda v: isqrt(v) if is_square(v) else None)


def xor_k_index(k=5):
    """Indexer for `n such that n ^ k is a square`. O(1) both directions.

    unrank(i) = i^2 ^ k        rank(n) = isqrt(n ^ k)

    Note this is GENERATION order (by i), not sorted order -- for k=5 the
    first three terms are 5, 4, 1, going down before they go up, because
    while i is small the XOR can flip a high bit off.
    """
    return Indexer(
        f"n ^ {k} is a square",
        unrank=lambda i: i * i ^ k,
        rank=lambda n: isqrt(n ^ k) if is_square(n ^ k) else None,
    )


def xor_k_x_index(k=5):
    """Indexer for the x-form: x square AND isqrt(x) ^ k square.

    x = (i^2 ^ k)^2, so |S below X| ~ X^(1/4) -- see count_below().
    """
    return Indexer(
        f"x square and isqrt(x) ^ {k} square",
        unrank=lambda i: (i * i ^ k) ** 2,
        rank=lambda x: (isqrt(isqrt(x) ^ k)
                        if is_square(x) and is_square(isqrt(x) ^ k) else None),
    )


def count_below(X, k=5):
    """Exact size of {x < X : x square, isqrt(x) ^ k square}.

    Closed form -- does NOT scan X candidates. Takes the fourth root and
    corrects by one or two for the XOR wobble near the boundary, so it costs
    the same whether X is 10^4 or 10^1000.

        count_below(10**100)  ->  10**25
    """
    u = lambda i: (i * i ^ k) ** 2
    m = isqrt(isqrt(X))
    while u(m) < X:
        m += 1
    while m > 0 and u(m - 1) >= X:
        m -= 1
    return m


def search_index(member, name="search", start=0):
    """Indexer for a sequence with no closed form. Walks, but caches.

    Use for the wild families -- e.g. m^2 ^ (m*isqrt(m)) being a square.
    """
    return Indexer(name, member=member, start=start)


if __name__ == "__main__":
    pd.set_option("display.width", 120)

    print("=== tame: n ^ 5 is a square (one unknown -> bijection) ===")
    print(xor_k_squares(5).to_string(), "\n")

    print("=== which k collapse cleanly? (small period = clean) ===")
    s = xor_k_survey(range(1, 17))
    print(s[["bits", "popcount", "period"]].to_string(), "\n")

    print("=== Kummer, one pair ===")
    for m in (4, 5):
        print("  ", kummer_check(m * m, 5))
    print()

    print("=== wild: m^2 ^ m is a square (two unknowns) ===")
    print("  solutions below 200000:", xor_self_squares())
    print("  why -- the window at m=1000:")
    w = xor_window(1000)
    print(f"    reach {w['reach']}, window {w['window']}, "
          f"squares inside: {w['squares_in_window']}\n")

    print("=== widened window, |delta| <= 8 ===")
    print(near_square_xor(limit=5000).head(12).to_string(), "\n")

    print("=== A221643 family: k^2 ^ (k+1)^2 square (both operands large) ===")
    print("  ", consecutive_xor_squares(300000)[:14])


# ------------------------------------------------ the square rail
def master_square(lo=0, hi=130, k=5, restart=True):
    """master_rail, but walking the XOR-SQUARE sequence instead of the primes.

    Members are n with n ^ k a perfect square -- i.e. n = m^2 ^ k. The rails
    lag exactly as in nt.master_rail: they flip only when a member goes by.

        member   is x itself in the sequence?
        count    (members strictly below x) - 1, clamped at 0   [lag 0]
        counn    count - 1                                      [lag 1]
        railn    count - 2                                      [lag 2]
        rodd     2*railn + 1  -- the odd-only rail 1,3,5,7,...
        run      plain running count of members seen so far (no lag)
        to_next  how many integers to step before the NEXT member lands
        since    how many since the last one
        m        the index, isqrt(x ^ k), when x is a member
    """
    members = set()
    m = 0
    while True:
        v = m * m ^ k
        if v > hi + 4096:
            break
        members.add(v)
        m += 1
    ms = sorted(members)

    rows, count, run = [], 0, 0
    start = int(lo) if restart else 0
    last = None
    for x in range(start, int(hi)):
        is_mem = x in members
        railn = max(0, count - 2)
        nxt = next((v for v in ms if v >= x), None)
        if x >= lo:
            rows.append(dict(
                count=count, counn=max(0, count - 1), railn=railn,
                rodd=2 * railn + 1, run=run, x=x, member=is_mem,
                to_next=(None if nxt is None else nxt - x),
                since=(None if last is None else x - last),
                m=(isqrt(x ^ k) if is_mem else None),
                x_xor_k=x ^ k, sq=is_square(x ^ k),
            ))
        if is_mem:
            count = run            # lags like nt.master_rail
            run += 1
            last = x
    return pd.DataFrame(rows)


def square_walk(hi=16430, k=5):
    """Your original loop as a DataFrame -- INCLUDING the rows it skipped,
    and keeping every column your print() had, even where the test fails.

        x         j^2, the candidate
        xm2_xor   (x-2) ^ k
        j_xor     isqrt(x) ^ k   <- the thing that must be a square
        sqrt_j    math.sqrt(j_xor)  -- a FLOAT, so a miss shows as 2.236...
        sq_j      is j_xor a square?   <- the only test that ever fails
        sq_x      is x a square?       <- always True here, by construction
        sqrt_x    math.sqrt(x)
        count     running hit index; blank between hits
        to_next   squares to step before the next hit
    """
    from math import sqrt
    rows, count = [], 0
    cand, j = [], 0
    while j * j < hi:
        cand.append(j)
        j += 1
    hits = {t for t in cand if is_square(t ^ k)}
    coun, nbelow, last_m = 0, 0, None
    for t in cand:
        x, v = t * t, t ^ k
        ok = is_square(v)
        railn = max(0, coun - 2)
        nxt = next((u for u in cand if u >= t and u in hits), None)
        rows.append({"j": t, "coun": coun, "counn": max(0, coun - 1),
                     "railn": railn, "rodd": 2 * railn + 1,
                     "x": x, "xm2_xor": (x - 2) ^ k, "j_xor": v,
                     "sqrt_j": round(sqrt(v), 3), "sq_j": ok,
                     "sq_x": is_square(x), "sqrt_x": float(t),
                     "delta": v - t,
                     "miss": v - isqrt(v) ** 2 if not ok else 0,
                     "since": None if last_m is None else t - last_m,
                     "last_mem": last_m, "next_mem": nxt,
                     "nth": count if ok else None,
                     "to_next": (None if nxt is None else
                                 sum(1 for u in cand if t <= u < nxt))})
        if ok:
            count += 1
            coun = nbelow          # lags exactly like nt.master_rail
            nbelow += 1
            last_m = t
    return pd.DataFrame(rows).set_index("j")


def master_rail_sq(lo=0, hi=128, k=5, restart=True):
    """nt.master_rail with the SQUARE sequence joined on -- count used as an
    index into our squares, and the six number columns pinned to the left.

    Each row's `count` looks up the count-th member of the XOR-square list:

        s_nth   the index used (== count)
        s_x     that member's x        1, 16, 25, 144, 441, ...
        s_xm2   (x-2) ^ k             -6, 11, 18, 139, 434, ...
        s_jxor  isqrt(x) ^ k           4,  1,  0,   9,  16, ...
        s_sqrtj sqrt(s_jxor)          2.0, 1.0, 0.0, 3.0, 4.0
        s_sqrtx sqrt(x)               1.0, 4.0, 5.0, 12.0, 21.0

    Values repeat while count sits still -- that is the point, it shows how
    the prime rail and the square rail line up.

    restart=True  (default) numbers the rails LOCALLY -- count starts at 0 at
                  lo, so a window like (96, 132) restarts the square lookup
                  from the first member.
    restart=False numbers them GLOBALLY from 0, so the window shows where you
                  actually are in both sequences.
    """
    from math import sqrt
    from sitoolz import numtools as _nt

    df = _nt.master_rail(lo, hi, restart=restart).copy()

    # the square sequence, in hit order
    hits, j = [], 0
    while len(hits) < int(df["count"].max()) + 2:
        if is_square(j ^ k):
            x = j * j
            hits.append((x, (x - 2) ^ k, j ^ k, round(sqrt(j ^ k), 3), float(j)))
        j += 1

    def pick(c, i):
        return hits[c][i] if 0 <= c < len(hits) else None

    cols = ["s_x", "s_xm2", "s_jxor", "s_sqrtj", "s_sqrtx"]
    df.insert(0, "s_nth", df["count"])
    for n, i in zip(cols, range(5)):
        df.insert(1 + i, n, [pick(c, i) for c in df["count"]])
    return df


# ================================================================ fast square wheel (opt-in)
# A prime-stacked residue pre-filter for perfect squares.  A square can only land on certain
# leftovers (quadratic residues) mod M, so if n's residue is not one of them, n is definitely
# NOT a square -- reject it for free, before the costly isqrt.  Each small prime in M roughly
# halves the non-squares that slip through.  These are OPT-IN (the `_fast` names) so you always
# know the wheel is on; the plain is_square() above is unchanged.
import functools as _functools
from math import isqrt as _isqrt

_WHEEL_M = 64 * 9 * 5 * 7 * 11          # = 221760 ; rejects ~98% of non-squares with one % M


@_functools.lru_cache(maxsize=8)
def square_wheel(M=_WHEEL_M):
    """Build (and cache) the perfect-square residue filter: (M, table), table[r]==1 iff r is
    a quadratic residue mod M.  Pass a bigger M (more small-prime factors) for a stronger one."""
    table = bytearray(M)
    for j in range(M):
        table[(j * j) % M] = 1
    return M, table


def is_square_fast(n, wheel=None):
    """Perfect-square test with the prime-wheel pre-filter (opt-in).  Same answer as
    is_square(), ~3x faster when scanning and much more on large n (the % M stays cheap
    while isqrt does not).  Pass a custom `wheel` from square_wheel() to tune the filter."""
    n = int(n)
    if n < 0:
        return False
    M, table = wheel or square_wheel()
    if not table[n % M]:                 # not a quadratic residue -> not a square, no isqrt
        return False
    r = _isqrt(n)
    return r * r == n


def is_member_fast(x, k=2, wheel=None):
    """Membership in the (j^2 XOR k)^2 family with the wheel (opt-in): reject both square
    checks by residue before any isqrt.  Returns the index j if x is a member, else None."""
    x = int(x)
    if x < 0:
        return None
    M, table = wheel or square_wheel()
    if not table[x % M]:                 # x can't be a perfect square
        return None
    r = _isqrt(x)
    if r * r != x:
        return None
    w = r ^ int(k)
    if w < 0 or not table[w % M]:         # r ^ k can't be a perfect square
        return None
    j = _isqrt(w)
    return j if j * j == w else None
