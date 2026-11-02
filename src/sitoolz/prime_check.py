#!/usr/bin/env python3
"""prime_check.py -- prove a number prime (or find a factor) using a saved base.

The base is a set of primes built ONCE by elimination (no oracle) and saved to
disk. With primes up to P you can PROVE any number up to P**2, because any
composite <= P**2 must have a factor <= P sitting in the base. The default base
reaches ~2.8 trillion; it loads from disk if present, else builds and saves it.

Usage:
    python prime_check.py                 # tests 6777213
    python prime_check.py 2800000000099   # tests any number up to the reach
    python prime_check.py 6777213 --grow  # rebuild a bigger base first
"""
import os
import sys
import time

import numpy as np

BASE_FILE = os.path.expanduser("~/sitoolz/base_primes.npy")
BASE_LIMIT = 1_700_000          # primes up to here -> reach ~2.89 trillion (~1 MB base)


def build_base(limit):
    """Sieve of Eratosthenes: primes up to `limit`, built by elimination (no test)."""
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if s[i]:
            s[i * i :: i] = False
    return np.flatnonzero(s)


def load_or_build_base(limit=BASE_LIMIT, path=BASE_FILE, force=False):
    """Load the saved base if it exists and is big enough, else build + save it."""
    if not force and os.path.exists(path):
        base = np.load(path)
        if int(base[-1]) >= limit - 1:
            return base
    base = build_base(limit)
    np.save(path, base)
    return base


def check(n, base):
    """Return (verdict, witness):
       verdict = True (prime), False (composite), or None (base too small to prove).
       witness = the base prime that divides n, if composite."""
    if n < 2:
        return False, None
    for p in base:
        if p * p > n:          # got past sqrt(n): nothing left can divide it
            return True, None
        if n % p == 0:         # a base prime divides n
            return False, int(p)
    return None, None          # ran out of base before sqrt(n) -> can't prove


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    n = int(args[0]) if args else 6_777_213
    force = "--grow" in sys.argv

    t = time.perf_counter()
    base = load_or_build_base(force=force)
    reach = int(base[-1]) ** 2
    print(f"base: {len(base):,} primes up to {int(base[-1]):,}   reach = {reach:,}   "
          f"({(time.perf_counter()-t)*1000:.0f} ms to load)")

    if n > reach:
        need = int(n ** 0.5) + 1
        print(f"  {n:,} is past the reach; rerun with a bigger BASE_LIMIT (>= {need:,}) and --grow")
        return

    t = time.perf_counter()
    verdict, witness = check(n, base)
    dt = (time.perf_counter() - t) * 1000
    if verdict is True:
        print(f"  {n:,} is PRIME   (no base prime <= {int(n**0.5):,} divides it)   [{dt:.2f} ms]")
    elif verdict is False:
        print(f"  {n:,} is COMPOSITE   {n:,} = {witness} x {n // witness:,}   [{dt:.2f} ms]")
    else:
        print(f"  couldn't prove {n:,} -- base too small")


if __name__ == "__main__":
    main()
