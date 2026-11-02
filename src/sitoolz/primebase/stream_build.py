"""Segmented, streaming base builder -- flat RAM no matter how big the base.

The monolithic :func:`primebase.build.build_base` needs ``limit`` bytes for the
sieve *plus* an int64 array of every prime (``8 * pi(limit)``) -- that's what put a
2e10 base at ~17 GB. This builder never holds either:

* only the primes up to ``sqrt(limit)`` stay resident (41k of them even at 3e11),
* one segment buffer is sieved at a time,
* half-gaps are written straight into a memory-mapped ``.npy`` as they're produced.

So RAM is flat (~segment size) and the ceilings become **disk** (~1 byte/prime) and
the uint8 encoding (gaps must stay < 512, which holds until ``limit ~ 3.03e11``).

    from primebase.stream_build import build_streamed
    build_streamed(3 * 10**11, "big_base")     # -> big_base.npy + big_base.meta.json
"""
from __future__ import annotations

import json
import math
import os
from math import isqrt

import numpy as np

from .build import build_base
from .counting import li


def build_streamed(limit: int, stem: str, seg: int = 1 << 30, verbose: bool = True) -> tuple:
    """Sieve to ``limit`` in segments, streaming uint8 half-gaps to ``<stem>.npy``.

    ``seg`` is the segment width in integers (default ~1e9 -> a 1 GB bool buffer).
    Returns ``(npy_path, meta_path)``. Raises if a gap reaches 512 (uint8 overflow),
    which starts happening past ``limit ~ 3.03e11``.
    """
    limit = int(limit)
    if limit < 3:
        raise ValueError("limit must be >= 3")
    root = isqrt(limit) + 1
    small = build_base(root)                      # primes <= sqrt(limit): stays tiny
    if verbose:
        print(f"base primes up to sqrt({limit:.3e}) = {root:,}: {small.size:,} primes", flush=True)

    # upper bound on pi(limit): x/ln(x) UNDERSHOOTS pi(x) (by ~8% at 1e6), so use the
    # logarithmic integral, which overshoots pi(x) everywhere in reach, plus a margin.
    upper = int(li(limit) * 1.01) + 1000
    npy, meta_path = f"{stem}.npy", f"{stem}.meta.json"
    out = np.lib.format.open_memmap(npy, mode="w+", dtype=np.uint8, shape=(upper,))

    pos = count = 0
    last = None                                   # last odd prime emitted
    largest = 2
    for lo in range(0, limit + 1, seg):
        hi = min(lo + seg, limit + 1)
        buf = np.ones(hi - lo, dtype=bool)
        if lo == 0:
            buf[0] = buf[1] = False
        for p in small:
            p = int(p)
            if p * p >= hi:
                break
            start = max(p * p, ((lo + p - 1) // p) * p)
            if start < hi:
                buf[start - lo::p] = False
        ps = np.flatnonzero(buf).astype(np.int64) + lo
        if ps.size == 0:
            continue
        if ps[0] == 2:                            # the only even prime: counted, no gap
            count += 1
            ps = ps[1:]
        if ps.size == 0:
            continue
        if last is None:                          # seed the odd chain at 3
            last = int(ps[0])
            count += 1
            ps = ps[1:]
            if ps.size == 0:
                continue
        hg = np.diff(np.concatenate(([last], ps))) // 2
        if hg.max() >= 256:
            raise ValueError(f"prime gap >= 512 near {int(ps[hg.argmax()]):,}: uint8 half-gaps overflow")
        out[pos:pos + hg.size] = hg.astype(np.uint8)
        pos += hg.size
        count += ps.size
        last = largest = int(ps[-1])
        if verbose:
            print(f"  {hi/limit*100:5.1f}%  primes {count:,}  largest {largest:,}", flush=True)

    out.flush()
    del out
    meta = {"largest": largest, "count": count, "reach": largest ** 2,
            "format": "uint8-halfgaps-from-3", "halfgaps": pos,
            "note": "primes[0]=2, primes[1:]=3+2*cumsum([0]+halfgaps[:count-2])"}
    json.dump(meta, open(meta_path, "w"), indent=2)
    if verbose:
        print(f"done: {count:,} primes, largest {largest:,}, reach {largest**2:.3e}", flush=True)
        print(f"  {npy} ({os.path.getsize(npy)/1e9:.2f} GB, {pos:,} half-gaps used of {upper:,})", flush=True)
    return npy, meta_path
