"""Make your own base -- bigger reach for more RAM, or a raw array for mmap.

The shipped base reaches 72 quadrillion in 15 MB. If you want more, sieve a
bigger limit and pack it. Reach grows as ``limit**2`` while storage grows only as
``limit``, so each doubling of the limit quadruples the reach.

    from primebase.build import build_base, pack
    p = build_base(1_800_000_000)          # ~2 GB RAM, ~30 s
    pack(p, "big_base")                    # writes big_base.npy + big_base.meta.json

Use it by pointing primebase at the two files before importing::

    export PRIMEBASE_GAPS=/path/big_base.npy
    export PRIMEBASE_META=/path/big_base.meta.json

or load it directly with :func:`unpack`.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def build_base(limit: int) -> np.ndarray:
    """Primes ``<= limit`` by a vectorized Sieve of Eratosthenes.

    RAM is ~``limit`` bytes for the boolean sieve (e.g. 1e10 -> ~10 GB). For
    limits that won't fit in memory, sieve in segments instead.
    """
    limit = int(limit)
    if limit < 2:
        return np.empty(0, dtype=np.int64)
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = False
    return np.flatnonzero(sieve).astype(np.int64)


def pack(primes: np.ndarray, stem: str) -> tuple[str, str]:
    """Write ``primes`` to ``<stem>.npy`` (uint8 half-gaps) + ``<stem>.meta.json``.

    Returns the two paths. Half-gaps stay in a byte while gaps are < 512, which
    holds well past 10^12; a larger gap raises rather than silently truncating.
    """
    primes = np.asarray(primes, dtype=np.int64)
    if primes[0] != 2 or primes[1] != 3:
        raise ValueError("expected a prime list starting 2, 3, 5, ...")
    halfgaps = np.diff(primes[1:]) // 2
    if halfgaps.max() >= 256:
        raise ValueError("gap too large for uint8 half-gaps; use a larger base only up to ~1e12")
    npy = f"{stem}.npy"
    meta_path = f"{stem}.meta.json"
    np.save(npy, halfgaps.astype(np.uint8))
    json.dump(
        {"largest": int(primes[-1]), "count": int(primes.size),
         "reach": int(primes[-1]) ** 2, "format": "uint8-halfgaps-from-3",
         "note": "primes[0]=2, primes[1:]=3+2*cumsum([0]+halfgaps)"},
        open(meta_path, "w"), indent=2,
    )
    return npy, meta_path


def unpack(stem: str) -> np.ndarray:
    """Decode a packed ``<stem>.npy`` back to the full int64 prime array."""
    hg = np.load(f"{stem}.npy").astype(np.int64)
    odds = 3 + 2 * np.cumsum(hg)
    return np.concatenate(([2, 3], odds))


def export_raw(primes: np.ndarray, path: str) -> str:
    """Save the plain prime array (int64) for ``np.load(path, mmap_mode='r')``.

    8x bigger than the packed form, but memory-mappable: sessions touch only the
    pages they actually read.
    """
    np.save(path, np.asarray(primes, dtype=np.int64))
    return path
