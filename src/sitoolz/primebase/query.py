"""Ask the base for primes: ranges, the first N, the k-th, and twin pairs.

``primes(lo, hi)`` is the powerful one: it doesn't just read stored primes, it
*generates* every prime in ``[lo, hi)`` with a segmented sieve seeded by the base.
Because the base holds every prime up to ``sqrt(REACH)``, that window can sit
anywhere up to the reach -- e.g. a million primes starting at a trillion.
"""
from __future__ import annotations

import math

import numpy as np

from ._base import _CHUNK, REACH, count, stream


def first(n: int) -> np.ndarray:
    """The first ``n`` primes: ``first(31) -> [2, 3, 5, ..., 127]``.

    Decodes only the ``n``-prime prefix, so this stays cheap on a base of any size.
    """
    n = int(n)
    if n < 0:
        raise ValueError("n must be >= 0")
    if n > count():
        raise ValueError(
            f"base has {count():,} primes; for more, build a bigger base "
            f"(see primebase.build.build_base / pack)"
        )
    if n == 0:
        return np.empty(0, dtype=np.int64)
    out, total = [], 0
    for c in stream(chunk=max(1, min(n, _CHUNK))):
        out.append(c)
        total += c.size
        if total >= n:
            break
    return np.concatenate(out)[:n]


def nth(k: int) -> int:
    """The ``k``-th prime, 1-indexed: ``nth(1) == 2``, ``nth(1_000_000)``.

    Walks to the ``k``-th prime and stops; nothing past it is decoded.
    """
    k = int(k)
    if not 1 <= k <= count():
        raise ValueError(f"k must be in 1..{count():,}")
    seen = 0
    for c in stream(chunk=max(1, min(k, _CHUNK))):
        if seen + c.size >= k:
            return int(c[k - seen - 1])
        seen += c.size
    raise AssertionError(f"base metadata claims {count():,} primes but ran out at {seen:,}")


def primes(lo: int, hi: int) -> np.ndarray:
    """Every prime in the half-open range ``[lo, hi)``.

    Uses a segmented sieve seeded by the base, so ``hi`` may be anywhere up to
    ``REACH`` (72 quadrillion for the shipped base). The seeds are streamed, so RAM
    is set by the window ``hi - lo``, not by the base::

        primes(10**12, 10**12 + 1000)     # primes just past a trillion

    Time scales with ``pi(sqrt(hi))`` -- one sieve pass per seed prime -- so a window
    near a big base's reach costs billions of passes and is best kept modest.
    """
    lo, hi = max(2, int(lo)), int(hi)
    if hi > REACH:
        raise ValueError(f"hi={hi:,} exceeds the base's provable reach {REACH:,}")
    if lo >= hi:
        return np.empty(0, dtype=np.int64)

    root = math.isqrt(hi - 1)
    sieve = np.ones(hi - lo, dtype=bool)
    for c in stream(root):
        for p in c:
            p = int(p)
            if p * p >= hi:
                break
            start = max(p * p, ((lo + p - 1) // p) * p)   # first multiple of p in [lo,hi), >= p^2
            if start < hi:
                sieve[start - lo::p] = False
    return np.flatnonzero(sieve).astype(np.int64) + lo


def twins(lo: int, hi: int) -> np.ndarray:
    """Twin-prime pairs ``(p, p+2)`` with ``lo <= p`` and ``p + 2 < hi``.

    Returns an ``(m, 2)`` array; ``twins(0, 100)`` -> ``[[3,5],[5,7],[11,13],...]``.
    """
    ps = primes(lo, min(hi + 2, REACH + 1))
    gap2 = np.flatnonzero(np.diff(ps) == 2)
    pairs = np.stack([ps[gap2], ps[gap2] + 2], axis=1)
    return pairs[pairs[:, 1] < hi]
