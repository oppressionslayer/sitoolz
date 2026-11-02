"""How many primes -- exact counts, the octave law, and the Li estimate.

Three tiers of "how many primes in a range":

* :func:`pi` / :func:`count_range` -- **exact**. Instant from the base up to the
  largest stored prime; a Lucy_Hedgehog count (``O(x^0.75)``, no full prime list)
  beyond that.
* :func:`octaves` -- exact primes per power-of-two octave ``[2^k, 2^(k+1))``, which
  grow ~``2^k/(k ln2)`` and very nearly double each step.
* :func:`li` -- the logarithmic integral, an estimate so tight the error is well
  under 0.01% by k~18.
"""
from __future__ import annotations

import math

from ._base import REACH, count, count_upto, largest
from .query import primes

_GAMMA = 0.5772156649015328606
_FLOAT_Z = 700.0   # above this, li(x) ~ x/ln(x) itself overflows a float (e**710 > 1.8e308)


def _li_big(x: int, z: float) -> int:
    """``li(x)`` for ``x`` past float range, as an int, via the asymptotic series.

    ``li(x) ~ (x/z) * sum_{k>=0} k!/z**k`` with ``z = ln x``. That series is divergent
    -- its terms only start growing again at ``k ~ z`` -- so it is truncated at the
    smallest term, which leaves a relative error around ``e**-z``: utterly negligible
    once ``z`` is in the hundreds. For a 2441-digit ``x``, ``z`` is 5620 and the sum
    settles after about six terms.

    Only the small factor ``C/z`` is computed in float; multiplying it into ``x``
    happens in integer arithmetic, so the answer keeps ``x``'s full magnitude and
    carries ~16 significant digits.
    """
    C, term = 1.0, 1.0
    for k in range(1, int(z) + 1):
        term *= k / z                       # k!/z**k
        if term < 1e-18:                    # smallest term reached: stop
            break
        C += term
    return (x * int(C / z * (1 << 64))) >> 64


def li(x):
    """Logarithmic integral ``Li(x) = Ei(ln x)`` -- the classic prime-count estimate.

    Returns a ``float`` while one can hold the answer, and an ``int`` beyond that:
    ``pi`` of a 2000-digit number is a 2000-digit number, and no float will take it.

    The convergent (Ramanujan) series used below sums ``z**k/(k*k!)``, whose terms
    *peak* at ``k ~ z = ln x`` and only then decay -- so the loop must run well past
    ``z`` or it truncates mid-climb and silently returns a low answer. The bound
    scales with ``z`` for exactly that reason.
    """
    if isinstance(x, float):
        if x <= 1.0:
            return 0.0
        z = math.log(x)
    else:
        x = int(x)
        if x <= 1:
            return 0.0
        z = math.log(x)                     # math.log takes big ints natively

    if z >= _FLOAT_Z:
        return _li_big(int(x), z)

    total, term = 0.0, 1.0
    for k in range(1, max(400, int(4 * z) + 100)):
        term *= z / k                 # z**k / k!
        total += term / k             # z**k / (k * k!)
        if k > z and abs(term / k) < 1e-16 * (abs(total) + 1):
            break                     # only after the peak -- see the docstring
    return _GAMMA + math.log(z) + total


def prime_count(x: int) -> int:
    """Exact ``pi(x)`` by the Lucy_Hedgehog method -- ``O(x^0.75)``, no prime list.

    Counts primes with no table at all, so it works past the stored base (feasible
    to ~1e12 in Python). :func:`pi` calls this only when the base can't answer.
    """
    x = int(x)
    if x < 2:
        return 0
    r = math.isqrt(x)
    vals = [x // i for i in range(1, r + 1)]
    vals += list(range(vals[-1] - 1, 0, -1))
    S = {v: v - 1 for v in vals}
    for p in range(2, r + 1):
        if S[p] > S[p - 1]:            # p is prime
            sp, p2 = S[p - 1], p * p
            for v in vals:
                if v < p2:
                    break
                S[v] -= S[v // p] - sp
    return S[x]


def pi(x: int) -> int:
    """Exact number of primes ``<= x``.

    Counted off the base's half-gaps while ``x`` is within the stored primes (see
    :func:`primebase._base.count_upto` -- it walks, so cost tracks ``x``, and ``x``
    at or past the last stored prime is instant from the metadata); an exact
    Lucy_Hedgehog count beyond. Use :func:`li` for a fast estimate on very large ``x``.
    """
    x = int(x)
    if x < 2:
        return 0
    if x == largest():
        return count()                  # the whole base, straight from the metadata
    if x < largest():
        return count_upto(x)
    return prime_count(x)


def count_range(lo: int, hi: int) -> int:
    """Exact number of primes in the half-open range ``[lo, hi)``."""
    lo, hi = max(2, int(lo)), int(hi)
    if lo >= hi:
        return 0
    if hi <= largest() + 1:                          # both inside the stored base
        return count_upto(hi - 1) - count_upto(lo - 1)
    if hi <= REACH and (hi - lo) <= 5_000_000:       # narrow window -> segmented sieve
        return int(primes(lo, hi).size)
    return pi(hi - 1) - pi(lo - 1)                    # wide range -> Lucy


def octaves(k_lo: int, k_hi: int) -> list:
    """Primes per octave ``[2^k, 2^(k+1))`` for ``k`` in ``[k_lo, k_hi]``.

    Each record: ``{k, exact, li_est, law, ratio}`` where ``law = 2^k/(k ln2)`` and
    ``ratio`` is exact/previous (trends to 2). ``exact`` is ``None`` if the octave
    runs past the exact-count range and only the estimate is available.
    """
    out, prev = [], None
    for k in range(int(k_lo), int(k_hi) + 1):
        lo, hi = 1 << k, 1 << (k + 1)
        try:
            exact = count_range(lo, hi)
        except Exception:
            exact = None
        out.append({
            "k": k,
            "exact": exact,
            "li_est": li(hi) - li(lo),
            "law": lo / (k * math.log(2)),
            "ratio": (exact / prev) if (exact and prev) else None,
        })
        prev = exact
    return out
