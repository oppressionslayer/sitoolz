"""primebase -- a whole prime range from a tiny net.

You don't store the primes; you store the *factors*. The shipped 15 MB base holds
every prime up to 268,499,993 (= sqrt of the reach), and from that it will:

* **prove** any number up to REACH = 72,092,246,241,000,049 prime -- or factor it,
* **generate** every prime in any window up to REACH (a segmented sieve seeded by
  the base -- e.g. the primes just past a trillion),
* hand you the first N primes, the k-th prime, and twin-prime pairs.

Quickstart::

    import primebase as pb

    pb.first(31)                       # array([2, 3, 5, ..., 127])
    pb.primes(10**12, 10**12 + 100)    # primes just past a trillion
    pb.twins(0, 100)                   # twin pairs [[3,5],[5,7],[11,13],...]
    pb.is_prime(9_999_999_999_999_917) # True (proven)
    pb.factor(72_092_246_240_999_997)  # 3
    pb.REACH                           # 72092246241000049

Nothing loads at import; the base is read on first call. Small queries decode only
the prefix they need. Grow a bigger base with :mod:`primebase.build` (monolithic) or
:mod:`primebase.stream_build` (segmented -- flat RAM, so this is the one for a base
measured in GB), then point at it::

    export PRIMEBASE_GAPS=~/big_base.npy
    export PRIMEBASE_META=~/big_base.meta.json

**Memory.** Every function above is streaming: it holds one 8 MB chunk at a time and
stops as soon as it can, so RAM does not track the base's size. The two exceptions
are :func:`base` and :func:`base_upto`, which by definition hand back an array --
and a decoded prime costs 8 bytes, so an 11.8-billion-prime base is 94 GB. Both
refuse past ``PRIMEBASE_MAX_DECODE`` primes (default 2e8) rather than let the OOM
killer take your session; iterate :func:`stream` instead, which has no ceiling::

    for chunk in pb.stream(limit):     # int64 chunks, flat RAM at any base size
        ...

Note ``pb.mod(n, chunk)`` rather than ``n % chunk`` if you write your own loop: past
9.2e18 (int64) plain numpy raises ``OverflowError``, and a big base's reach is well
past that.
"""
from __future__ import annotations

from ._base import (REACH, base, base_upto, count, count_upto, largest, meta,
                    mod, stream)
from .certify import factor, factorize, is_prime, largest_factor
from .counting import count_range, li, octaves, pi, prime_count
from .neighborhood import neighborhood
from .prove import (build_up, certificate, gcd_certify, is_pocklington_prime,
                    odd_core, pocklington, pock_chain, pock_descend, report,
                    wrap_descend, wrap_scan,
                    quick_prove, roll_chain, self_prove,
                    verify_certificate)
from .query import first, nth, primes, twins
from .rollball import chain_descend, deep_descend, descend, roll_up, rungs

__version__ = "0.1.0"

__all__ = [
    # query
    "first", "nth", "primes", "twins",
    # certify / factor
    "is_prime", "factor", "factorize", "largest_factor",
    # counting
    "pi", "count_range", "octaves", "li", "prime_count",
    # neighborhood (factor anatomy)
    "neighborhood",
    # rollball (Pocklington ladder)
    "descend", "rungs", "roll_up", "chain_descend", "deep_descend",
    # proofs (Pocklington certificates)
    "is_pocklington_prime", "certificate", "verify_certificate", "roll_chain",
    "self_prove", "gcd_certify", "pocklington", "odd_core", "build_up",
    "quick_prove", "pock_descend", "pock_chain", "report", "wrap_descend", "wrap_scan",
    # base + facts
    "base", "base_upto", "stream", "count_upto", "mod",
    "count", "largest", "meta", "REACH",
]
