"""gap_forge.py -- construct prime gaps (the even 'distances') by sieve, not by testing.

Build, don't decode.  The primes and their gaps are CONSTRUCTED by eliminating composites
(Eratosthenes): you never ask "is this prime?" -- what survives elimination *is* prime.
And by the squaring bound, crossing out multiples of the primes up to sqrt(N) builds every
prime and every gap up to N exactly.  Feed it N; extract numpy arrays of primes, gaps
(the even distances), twins, gap counts, jumping champion.

    from sitoolz.gap_forge import GapForge
    f = GapForge(1_000_000)
    f.primes       # array of primes <= N          (constructed, not tested)
    f.gaps         # consecutive distances          (== your pg; gaps[0]=1, rest even)
    f.evens        # the strictly-even distances     (gaps[1:])
    f.twins        # lower member of each twin pair  (where gap == 2)
    f.where_gap(6) # indices where the distance is 6
    f.gap_counts() # {gap value: how many times it occurs}
"""
from __future__ import annotations
import numpy as np


def sieve(N):
    """Boolean prime mask for [0, N].  The constructor: strike out multiples of each prime
    up to sqrt(N) (the squaring bound); whatever survives is prime *by construction*."""
    N = int(N)
    mask = np.ones(N + 1, dtype=bool)
    if N >= 0:
        mask[:min(2, N + 1)] = False
    for p in range(2, int(N ** 0.5) + 1):
        if mask[p]:
            mask[p * p::p] = False
    return mask


def wheel_horizon(wheel_primes=(2, 3, 5, 7, 11, 13)):
    """The exact-build horizon of a wheel of the given primes: (next prime not in it)^2.

    A wheel of primes up to P builds every prime and every gap EXACTLY up to
    next_prime(P)**2 -- the squaring bound.  Feed it the next prime to reach the next square.
    """
    q = max(wheel_primes) + 1
    while any(q % p == 0 for p in range(2, q)):
        q += 1
    return q * q                          # next_prime(max)^2


class GapForge:
    """Construct the primes up to N and expose their gaps (distances / evens) as arrays."""

    def __init__(self, N):
        self.N = int(N)
        self._mask = sieve(self.N)
        self._primes = np.flatnonzero(self._mask)

    # ---- the built objects ----
    @property
    def primes(self):
        """All primes <= N (constructed by elimination, no primality test)."""
        return self._primes

    @property
    def gaps(self):
        """Consecutive prime gaps -- the distances.  gaps[i] = primes[i+1] - primes[i].
        gaps[0] == 1 (the sole odd distance, 2->3); every later distance is even."""
        return np.diff(self._primes)

    @property
    def evens(self):
        """The strictly-even distances (all gaps after the unique 2->3 = 1)."""
        return np.diff(self._primes)[1:]

    @property
    def cumulative(self):
        """Running distance from the first prime (== primes - 2)."""
        return self._primes - self._primes[0]

    # ---- queries ----
    def where_gap(self, g):
        """Indices i where the distance primes[i]->primes[i+1] equals g."""
        return np.flatnonzero(np.diff(self._primes) == int(g))

    @property
    def twins(self):
        """Lower member of each twin pair (distance == 2)."""
        g = np.diff(self._primes)
        return self._primes[:-1][g == 2]

    def gap_counts(self):
        """{gap value: number of occurrences} over all distances up to N."""
        vals, cnts = np.unique(np.diff(self._primes), return_counts=True)
        return dict(zip(vals.tolist(), cnts.tolist()))

    def jumping_champion(self):
        """The most common distance up to N."""
        vals, cnts = np.unique(np.diff(self._primes), return_counts=True)
        return int(vals[cnts.argmax()])

    def coverage(self, M, top_k=10):
        """How well the number M sieves the distances.

        A distance g is 'caught' when g divides M (M % g == 0) -- so M swallows every gap
        whose prime factorization fits inside its own.  Returns catch-rate and the top
        surviving (rough) gaps M cannot reach.  Use it to test any candidate sieve number
        against the real distances instantly.
        """
        M = int(M)                                  # may be a big Python int; keep it out of numpy
        g = np.diff(self._primes)
        vals, cnts = np.unique(g, return_counts=True)
        caught_val = np.array([M % int(v) == 0 for v in vals])   # bignum-safe: mod the few distinct gaps
        caught = int(cnts[caught_val].sum())
        total = int(g.size)
        sv, sc = vals[~caught_val], cnts[~caught_val]
        order = np.argsort(-sc)[:top_k]
        return {
            "M": M,
            "total_gaps": total,
            "caught": caught,
            "survived": total - caught,
            "catch_rate": caught / total,
            "top_survivors": [(int(sv[i]), int(sc[i])) for i in order],
        }

    def __repr__(self):
        return f"GapForge(N={self.N}, {len(self._primes)} primes, {len(self._primes)-1} gaps)"


if __name__ == "__main__":
    f = GapForge(1_000_000)
    print(f)
    print("first 20 primes :", f.primes[:20].tolist())
    print("first 20 gaps   :", f.gaps[:20].tolist())
    print("first 20 evens  :", f.evens[:20].tolist())
    print("twins <= 100    :", f.twins[f.twins <= 100].tolist())
    print("count of gap 6  :", len(f.where_gap(6)))
    print("jumping champion:", f.jumping_champion())
    print("gap counts (small):", {k: v for k, v in f.gap_counts().items() if k <= 14})
    print("wheel(<=13) builds gaps exactly up to:", wheel_horizon())
