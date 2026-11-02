"""Prove any number up to the reach prime, or hand back a factor.

This is a *proof*, not a probabilistic test: a composite ``n <= REACH`` must have
a prime factor ``<= sqrt(n) <= sqrt(REACH)``, and every such prime is in the base.
So if nothing in the base divides ``n``, ``n`` is prime -- guaranteed.

For ``n > REACH`` use a Miller-Rabin / BPSW tester (``sympy.isprime``,
``gmpy2.is_prime``); those need no table and run in microseconds. This module's
edge is that it also returns the *factor* of a composite, which those don't.
"""
from __future__ import annotations

import math

import numpy as np

from ._base import REACH, mod, stream

_TOO_BIG = (
    "n={n:,} exceeds the base's provable reach {reach:,}. "
    "Grow the base (primebase.build) or use a BPSW tester (sympy/gmpy2) for a "
    "yes/no answer with no table."
)


def factor(n: int):
    """Smallest prime factor of ``n``, or ``None`` if ``n`` is prime.

    ``n`` must be ``2 <= n <= REACH``. Vectorized trial division against only the
    base primes up to ``sqrt(n)``, streamed a chunk at a time and abandoned the
    moment a factor turns up -- so RAM is flat even when ``sqrt(n)`` reaches the far
    end of a 12 GB base, and composites usually exit in the first chunk.
    """
    n = int(n)
    if n < 2:
        raise ValueError("n must be >= 2")
    if n > REACH:
        raise ValueError(_TOO_BIG.format(n=n, reach=REACH))
    for c in stream(math.isqrt(n)):
        hits = c[mod(n, c) == 0]
        if hits.size:
            return int(hits[0])
    return None


def is_prime(n: int) -> bool:
    """``True`` iff ``n`` is prime, proven, for ``2 <= n <= REACH``."""
    n = int(n)
    if n < 2:
        return False
    if n > REACH:
        raise ValueError(_TOO_BIG.format(n=n, reach=REACH))
    return factor(n) is None


def factorize(n: int) -> dict:
    """Full prime factorization of ``n`` as ``{prime: exponent}``, for ``n <= REACH``.

    Peels every prime factor ``<= sqrt(n)`` from the base; the leftover cofactor is
    then necessarily a single prime (two factors above ``sqrt(n)`` would exceed
    ``n``). So the result is complete and exact.
    """
    n = int(n)
    if n < 1:
        raise ValueError("n must be >= 1")
    if n > REACH:
        raise ValueError(_TOO_BIG.format(n=n, reach=REACH))
    f: dict = {}
    for c in stream(math.isqrt(n)):
        if n == 1 or int(c[0]) ** 2 > n:   # nothing left, or leftover has no factor
            break                          # <= its own sqrt -> it is prime: stop early
        for p in c[mod(n, c) == 0]:
            p = int(p)
            while n % p == 0:
                f[p] = f.get(p, 0) + 1
                n //= p
    if n > 1:
        f[n] = f.get(n, 0) + 1        # remaining cofactor is prime
    return f


def largest_factor(n: int) -> int:
    """The largest prime factor of ``n`` (``1`` for ``n < 2``), for ``n <= REACH``."""
    n = int(n)
    if n < 2:
        return 1
    if n > REACH:
        raise ValueError(_TOO_BIG.format(n=n, reach=REACH))
    big = 1
    for c in stream(math.isqrt(n)):
        if n == 1 or int(c[0]) ** 2 > n:   # leftover is prime (or gone): stop early
            break
        for p in c[mod(n, c) == 0]:
            p = int(p)
            while n % p == 0:
                n //= p
            big = p                    # chunks ascend: last small divisor is the biggest
    return max(big, n) if n > 1 else big
