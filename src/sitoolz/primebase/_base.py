"""Loading the packed base primes.

The base ships as ``data/base_gaps.npy`` -- a uint8 array of *half-gaps* between
the odd primes, 8x smaller than storing the primes themselves. Three access paths,
in order of increasing safety:

* :func:`base`        -- decode the whole array once (cached in RAM).
* :func:`base_upto`   -- decode only the prefix up to a limit.
* :func:`stream`      -- yield the prefix in chunks and never materialize it.

Only :func:`stream` is safe at every base size. A decoded prime costs 8 bytes, so
a 12 GB base of 11.8 billion half-gaps is **94 GB** decoded -- :func:`base` and
:func:`base_upto` therefore refuse to hand back more than ``PRIMEBASE_MAX_DECODE``
primes (default 2e8, ~1.6 GB) rather than let the OOM killer take the session.
Library internals all use :func:`stream`, so they stay flat regardless.

Nothing is loaded at import time; the file is read on first use.
"""
from __future__ import annotations

import functools
import json
import os
from pathlib import Path

import numpy as np

_DATA = Path(__file__).resolve().parent / "data"

# Swap in a bigger base you built with primebase.build.pack by pointing these at
# your <stem>.npy / <stem>.meta.json files.
_GAPS = Path(os.environ.get("PRIMEBASE_GAPS", _DATA / "base_gaps.npy"))
_META = Path(os.environ.get("PRIMEBASE_META", _DATA / "meta.json"))
_CHUNK = 1 << 20  # max half-gaps decoded per step (8 MB of int64 -- keeps RAM tiny)
_RAMP = 1 << 9    # first block; doubles up to _CHUNK so shallow calls stay cheap

# Most primes base()/base_upto() will decode before refusing. 2e8 primes = 1.6 GB.
_MAX_DECODE = int(os.environ.get("PRIMEBASE_MAX_DECODE", 200_000_000))

_TOO_MANY = (
    "decoding {n:,}+ primes would need ~{gb:.1f} GB of RAM (8 bytes each) and the "
    "cap is {cap:,} (raise it with PRIMEBASE_MAX_DECODE if you really have the RAM). "
    "Iterate primebase.stream({lim}) instead -- same primes, in chunks, flat RAM."
)


@functools.lru_cache(maxsize=1)
def meta() -> dict:
    """The packed base's metadata: largest prime, count, provable reach."""
    return json.loads(_META.read_text())


@functools.lru_cache(maxsize=1)
def _halfgaps() -> np.ndarray:
    """The uint8 half-gap array, memory-mapped (never read into RAM).

    A streamed build (:mod:`primebase.stream_build`) over-allocates the file, so
    trust the count in the metadata rather than the file length.
    """
    arr = np.load(_GAPS, mmap_mode="r")
    m = meta()
    n = m.get("halfgaps", m["count"] - 2)     # shipped files hold exactly count-2
    return arr[:n]


def stream(limit: int = None, chunk: int = _CHUNK):
    """Yield base primes ``<= limit`` as int64 chunks -- flat RAM at any base size.

    This is the only decoder with no size ceiling: it holds one chunk (8 MB) at a
    time regardless of whether the base is 15 MB or 12 GB, and stops the moment
    ``limit`` is passed. Break out early and the rest is never touched::

        for c in stream(isqrt(n)):          # trial division, flat RAM
            hits = c[mod(n, c) == 0]
            if hits.size:
                return int(hits[0])

    ``limit=None`` walks the whole base. ``chunk`` is a *ceiling*: the first block is
    small and doubles up to it, so a shallow call (``is_prime`` on a small n needs
    only a few dozen primes) decodes a few hundred half-gaps rather than a million,
    while a full walk still reaches full-size blocks after a dozen doublings.
    """
    limit = None if limit is None else int(limit)
    if limit is not None and limit < 2:
        return

    head = np.array([2, 3], dtype=np.int64)
    if limit is not None:
        head = head[head <= limit]
    if head.size:
        yield head
    if limit is not None and limit < 5:
        return

    hg = _halfgaps()
    running, pos = 3, 0
    size = min(chunk, _RAMP)
    while pos < hg.size:
        c = np.asarray(hg[pos:pos + size], dtype=np.int64)
        decoded = running + 2 * np.cumsum(c)           # odd primes in this chunk
        running = int(decoded[-1])
        pos += size
        size = min(size * 2, chunk)                    # ramp up to the ceiling
        if limit is not None and decoded[-1] > limit:  # last chunk: trim and stop
            decoded = decoded[decoded <= limit]
            if decoded.size:
                yield decoded
            return
        yield decoded


def mod(n: int, primes: np.ndarray) -> np.ndarray:
    """``n % primes``, vectorized, for ``n`` of **any** size.

    Plain ``n % arr`` raises ``OverflowError`` once ``n`` passes int64 (9.2e18, 19
    digits) -- which the big base's 23-digit reach does. Above that this switches to
    Horner's rule over 20-bit limbs of ``n``: each step is ``r*2**20 + limb (mod p)``,
    and since base primes are ``< 2**40`` the widest intermediate is ``2**60``, so it
    stays inside int64 while ``n`` itself never has to.
    """
    n = int(n)
    if n < (1 << 62):
        return n % primes                     # fast path: numpy handles it natively

    limbs, t = [], n
    while t:
        limbs.append(t & 0xFFFFF)             # 20 bits at a time, little end first
        t >>= 20
    r = np.zeros_like(primes)
    for limb in reversed(limbs):              # Horner, most-significant limb first
        r = (r * (1 << 20) + limb) % primes
    return r


@functools.lru_cache(maxsize=1)
def base() -> np.ndarray:
    """Every base prime as one int64 array (decoded once, then cached).

    ~117 MB in RAM for the shipped base -- but 94 GB for a 3e11 base, so this
    raises past ``PRIMEBASE_MAX_DECODE`` rather than OOM the process. Prefer
    :func:`stream` (no ceiling) or :func:`base_upto` when you can.
    """
    n = count()
    if n > _MAX_DECODE:
        raise MemoryError(_TOO_MANY.format(n=n, gb=n * 8 / 1e9, cap=_MAX_DECODE, lim=""))
    return np.concatenate(list(stream()))


def base_upto(limit: int) -> np.ndarray:
    """Base primes ``<= limit``, decoding only the prefix needed (low RAM).

    For ``limit`` well below the base's largest prime this touches only a few MB,
    which is what keeps ``primes()`` and ``is_prime()`` cheap for small inputs. Bails
    with ``MemoryError`` past ``PRIMEBASE_MAX_DECODE`` primes -- use :func:`stream`
    for prefixes that big.
    """
    limit = int(limit)
    if limit < 2:
        return np.empty(0, dtype=np.int64)

    out, total = [], 0
    for c in stream(limit):
        total += c.size
        if total > _MAX_DECODE:
            raise MemoryError(
                _TOO_MANY.format(n=total, gb=total * 8 / 1e9, cap=_MAX_DECODE, lim=limit)
            )
        out.append(c)
    return np.concatenate(out) if out else np.empty(0, dtype=np.int64)


def count_upto(x: int) -> int:
    """How many base primes are ``<= x`` -- counted straight off the half-gaps.

    Never decodes a prime *value*: a chunk whose half-gaps sum to a running total
    still under ``x`` is entirely below ``x``, so it contributes its full length and
    only its ``sum`` is needed. Just one chunk -- the one straddling ``x`` -- gets a
    ``cumsum`` and a binary search. Flat RAM, and no int64 array is ever built.
    """
    x = int(x)
    if x < 2:
        return 0
    if x < 3:
        return 1
    if x < 5:
        return 2                            # 2 and 3, before any half-gap applies

    hg = _halfgaps()
    target = (x - 3) // 2                   # primes are 3 + 2*cumsum: cumsum <= target
    total, n, pos = 0, 2, 0                 # n counts 2 and 3 up front
    while pos < hg.size:
        blk = hg[pos:pos + _CHUNK]
        s = int(blk.sum(dtype=np.int64))    # uint8 summed straight from the memmap
        if total + s <= target:             # whole chunk sits at or below x
            total += s
            n += blk.size
            pos += _CHUNK
            continue
        cs = np.cumsum(np.asarray(blk, dtype=np.int64)) + total    # the straddling chunk
        return n + int(np.searchsorted(cs, target, side="right"))
    return n


def largest() -> int:
    """The largest prime actually stored (= sqrt of the provable reach)."""
    return int(meta()["largest"])


def count() -> int:
    """How many primes are in the base."""
    return int(meta()["count"])


REACH = int(meta()["reach"])  # largest number this base can PROVE prime (= largest**2)
