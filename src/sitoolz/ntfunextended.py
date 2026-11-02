"""
ntfunextended -- Lucas-Lehmer indexing / orbit toolkit (import as nte).

Everything we built exploring "indexing with LL":
  * the walkable, indexable ladder + step/unstep      ll_ladder, ll_step, ll_unstep
  * orbits on any modulus (zero vs cycle)             ll_orbit, ll_orbit_mode, ll_compare
  * the power-of-2 gem before zero                    ll_zero_row, ll_sqrt2, ll_pow2
  * the clean/shadow popcount law                     ll_clean_landing, ll_clean_scan
  * the binary tree-address of Z/M under LL           ll_preimages, ll_back_tree,
                                                       ll_address, ll_from_address
  * the map variants (XOR/scale twists)               ll_modes, ll_step_mode, LL_MODES
  * NEW: branching-arity = factor counter (algebraic) ll_arity, ll_arity_table

Usage:
    from sitoolz import ntfunextended as nte
    nte.ll_ladder(13)
    nte.ll_arity_table([7,11,13,17,19,23])

The headline law (this module's reason to exist):
    For a PURELY ALGEBRAIC LL map (squaring / powers, no XOR), the functional-graph
    branching arity equals  2**(#distinct prime factors of M)  times a power-map
    correction.  So arity == 2  <=>  M has exactly one prime factor  <=>  M is prime.
    XOR-containing variants BREAK this (they add irregular branching), because XOR is
    non-algebraic and does not respect CRT across the prime factors.
"""
from collections import Counter

from sitoolz import ntfun
# re-export the whole LL toolkit that already lives in ntfun
from sitoolz.ntfun import (                       # noqa: F401
    ll_step, ll_unstep, ll_ladder, ll_zero_row, ll_orbit, ll_compare,
    ll_sqrt2, ll_pow2, ll_clean_landing, ll_clean_scan,
    ll_preimages, ll_back_tree, ll_address, ll_from_address,
    ll_modes, ll_step_mode, ll_orbit_mode, LL_MODES,
)

try:
    import pandas as _pd
except ImportError:
    _pd = None


def ll_arity(M, mode="normal"):
    """Max branching arity of the LL-variant functional graph on Z/M.
    Builds the whole graph (O(M) time+memory) -- fine for small M, the honest
    price of reading global structure.  For algebraic maps this equals
    2**(#prime factors) (times a power-map factor); arity==2 <=> M prime."""
    f = LL_MODES[mode]
    indeg = Counter()
    for x in range(int(M)):
        indeg[f(x, M) % M] += 1
    return max(indeg.values())


def ll_arity_prime_test(M, mode="normal"):
    """Quick-ish primality verdict from branching arity (algebraic modes only):
    True iff arity == 2 (single prime factor).  EXACT for 'normal'/'pow5' etc.,
    but O(M) -- meant for SMALL M.  Not a shortcut for large numbers."""
    return ll_arity(M, mode) == 2


def ll_arity_table(ps, modes=("normal", "pow5", "sq_xor", "xor_insquare")):
    """Branching arity per (Mersenne exponent, mode), with #prime factors, so you
    can watch 2**(#factors) show through the algebraic maps and XOR break it."""
    if _pd is None:
        raise RuntimeError("pandas required")
    import sympy as _sp
    rows = []
    for p in ps:
        M = (1 << p) - 1
        row = {"p": p, "prime": bool(_sp.isprime(M)),
               "n_factors": len(_sp.factorint(M))}
        for m in modes:
            row[m] = ll_arity(M, m)
        rows.append(row)
    return _pd.DataFrame(rows).set_index("p")


def mersenne_arity_falsepos(pmax=20, mode="normal"):
    """QUICK false-positive check: over Mersenne numbers M_p = 2**p - 1 for p up to
    pmax, classify prime-vs-composite by the arity test (arity==2) and compare to
    truth.  Returns (false_positives, false_negatives, detail-list, seconds)."""
    import sympy as _sp
    import time as _t
    t0 = _t.time()
    fp, fn, detail = [], [], []
    for p in range(2, pmax + 1):
        M = (1 << p) - 1
        truth = bool(_sp.isprime(M))
        guess = ll_arity_prime_test(M, mode)
        detail.append((p, M, truth, guess))
        if guess and not truth:
            fp.append(p)          # said prime, actually composite  <-- false positive
        if truth and not guess:
            fn.append(p)          # said composite, actually prime
    return {"false_positives": fp, "false_negatives": fn,
            "detail": detail, "seconds": _t.time() - t0}
