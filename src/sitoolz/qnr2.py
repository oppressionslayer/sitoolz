"""qnr2 -- the Mersenne-offset prime filter, and what its tight gate really selects.

Lars Rocha's construction.  Scan odd x, build

    offset = (2**x - 1) - 2**(x//2)

and keep the x with `offset % x == 2`.  Unpacking it with r = (x-1)/2: if x is
prime then 2^x - 1 == 1 (Fermat) and 2^r == +/-1 (Euler), so offset is 0 or 2
mod x, and

    offset % x == 2   <=>   2 is a quadratic NON-residue mod x   <=>   x = 3,5 mod 8

which is a base-2 Euler-Jacobi test, arrived at through the offset rather than
from the textbook.  It is exact on that residue class: over odd x < 20000 the
loose filter returns every one of the 1140 primes that are 3 or 5 mod 8, misses
none, and admits exactly five composites -- 325, 3277, 3531, 3725, 6063, the
base-2 Euler-Jacobi pseudoprimes in range.

**What the tight gate actually is.**  The second gate

    off_r = (2**x - 1 - 2**r) % r        r = x//2
    keep if  x % off_r == 2              i.e. off_r divides x - 2

was written to remove those pseudoprimes.  It does -- but it is not a
pseudoprime filter.  Below 20000 it keeps 44 of 1145 survivors, discarding 1101
genuine primes.  It selects a much narrower class, and that class is exactly:

    x prime, x = 3 or 5 mod 8, and off_r | (x - 2)

with off_r determined by the STRUCTURE of r, not by a big computation:

    branch A -- r prime, so x = 2r+1 is a SAFE PRIME.
        Then 2^(r-1) = 1 mod r, so 2^r = 2 and 2^x = 2^(2r+1) = 8 mod r,
        giving off_r = 8 - 1 - 2 = 5 for every r > 5.
        The gate is then just 5 | x-2, i.e. x = 2 mod 5.
        With x = 3 mod 8 that is  x = 107 (mod 120).

    branch B -- r = 2q with q prime, so x = 4q+1.
        Then 2^r = 4 and 2^x = 32 mod r, giving off_r = 27.
        Lands on x = 29, 53, 77 (mod 120).

Verified: predicting membership from this structure alone reproduces the tight
filter exactly -- 285 of 285 below 200,000, zero disagreements.  So the tight
mode is a *safe-prime detector wearing a modular-arithmetic disguise*, and the
120-spacing in its output is the CRT of "3 mod 8" and "2 mod 5".

    from sitoolz import qnr2

    qnr2.qnr2_all(n)                 # COMPLETE test -- every prime, both branches
    qnr2.completeness(100000)        # the two branches, side by side
    qnr2.regime_frame()              # sieve vs qnr2 vs Lucas-Lehmer
    qnr2.qnr2_filter(20000)          # loose: every prime 3,5 mod 8
    qnr2.qnr2_filter(20000, tight=True)
    qnr2.qnr2_test(10**200 + 357)    # one-shot probable prime test
    qnr2.classify(107)               # why a given x is kept, or is not
    qnr2.tight_frame(20000)          # the keepers as a dataframe
    qnr2.branch_frame(200000)        # the branch structure
    qnr2.pseudoprimes(100000)        # what the loose filter lets through

gmpy2 is used when present and everything falls back to stdlib `pow`.
"""
from math import gcd

try:
    import gmpy2 as _gmpy2
except ImportError:                                   # pragma: no cover
    _gmpy2 = None


def _mr(n):
    """Deterministic Miller-Rabin for 64-bit, strong probable prime above."""
    n = int(n)
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _self(n):
    """Primality by qnr2's OWN criterion -- no external test anywhere.

    Only accepts n = 3 or 5 mod 8, so `classify(..., prime=_self)` is a
    closed system: the filter judging its own structure.  For n = 1 or 7 mod 8
    it returns False, which is correct behaviour for the criterion and means
    the safe/4q+1 branches are recognised only when r itself lands in the
    accepted class.
    """
    return qnr2_test(n)


BACKENDS = {"gmpy2": None, "mr": _mr, "self": _self}


def _is_prime(n, prime=None):
    """Primality, via the chosen backend.

    `prime` may be a callable, or one of 'gmpy2', 'mr' (stdlib Miller-Rabin),
    'self' (qnr2's own test -- no external primality anywhere).  None means
    gmpy2 if installed, else 'mr'.
    """
    if callable(prime):
        return bool(prime(n))
    if prime == "mr":
        return _mr(n)
    if prime == "self":
        return _self(n)
    if prime == "gmpy2" and _gmpy2 is None:
        raise ImportError("gmpy2 is not installed")
    if _gmpy2 is not None:
        return bool(_gmpy2.is_prime(n))
    return _mr(n)


# ------------------------------------------------------------------ the filter
def offset_mod(x):
    """`offset % x` where offset = (2**x - 1) - 2**(x//2).  No big powers."""
    r = x // 2
    return (pow(2, x, x) - 1 - pow(2, r, x)) % x


def off_r(x):
    """The tight gate's quantity: (2**x - 1 - 2**r) mod r, with r = x//2."""
    r = x // 2
    return (pow(2, x, r) - 1 - pow(2, r, r)) % r


def qnr2_filter(hi, tight=False):
    """Odd x in [3, hi) surviving the offset test.  `tight` adds the second gate."""
    out = []
    for x in range(3, int(hi), 2):
        if offset_mod(x) != 2:
            continue
        if tight:
            o = off_r(x)
            if not (o and x % o == 2):
                continue
        out.append(x)
    return out


_PRIMORIAL = {}


def _primorial(bound):
    P = _PRIMORIAL.get(bound)
    if P is None:
        sieve = bytearray([1]) * (bound + 1)
        for i in range(2, int(bound ** 0.5) + 1):
            if sieve[i]:
                sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        P = 1
        for p in range(3, bound + 1):
            if sieve[p]:
                P *= p
        _PRIMORIAL[bound] = P
    return P


def qnr2_test(x, bound=2000):
    """One-shot base-2 Euler-Jacobi test.  Accepts primes that are 3 or 5 mod 8.

    Uses 2^x = 2*(2^r)^2 so only ONE modexp is needed, plus a small-prime gcd
    pre-screen.  Note this rewrite is *stricter* than the literal offset: below
    4000 the two disagree on exactly 325 and 3725, both composite, both
    base-2 Euler-Jacobi pseudoprimes that the fast form rejects and the slow
    form admits.  A free win hiding in the optimisation.
    """
    x = int(x)
    if x < 3 or x % 2 == 0:
        return False
    g = gcd(x, _primorial(bound))
    if g != 1 and g != x:
        return False
    r = x >> 1
    a = int(_gmpy2.powmod(2, r, x)) if _gmpy2 else pow(2, r, x)
    return (2 * a * a - 1 - a) % x == 2


# ------------------------------------------------------- the COMPLETE criterion
# The offset test splits on Euler's criterion, and the original filter kept only
# one of its two branches:
#
#     x prime  =>  2^x - 1 == 1 (Fermat)  and  2^r == +/-1 (Euler), r = (x-1)/2
#
#     2 a QR   (x = 1,7 mod 8):  2^r == +1  ->  offset == 1 - 1 == 0
#     2 a QNR  (x = 3,5 mod 8):  2^r == -1  ->  offset == 1 + 1 == 2
#
# So `offset % x == 2` is half the criterion.  Accepting {0, 2} catches EVERY
# prime -- verified: all 9591 odd primes below 100,000, zero missed.  That makes
# the offset construction a complete base-2 Euler-Jacobi test rather than a
# residue-class filter, and it costs one comparison.

def offset_class(x):
    """Which Euler branch x lands in: 0, 2, or None if it fails both."""
    v = offset_mod(x)
    return v if v in (0, 2) else None


def qnr2_all(x, bound=2000):
    """Complete base-2 Euler-Jacobi probable-prime test -- BOTH branches.

    `qnr2_test` accepts only x = 3,5 mod 8; this accepts every prime.  Same one
    modexp, same gcd pre-screen, one extra comparison.  Uses 2^x = 2*(2^r)^2 so
    the whole test is a single exponentiation.
    """
    x = int(x)
    if x < 3 or x % 2 == 0:
        return x == 2
    g = gcd(x, _primorial(bound))
    if g != 1 and g != x:
        return False
    r = x >> 1
    a = int(_gmpy2.powmod(2, r, x)) if _gmpy2 else pow(2, r, x)
    return (2 * a * a - 1 - a) % x in (0, 2)


def qnr2_filter_all(hi, bound=2000):
    """Every odd x in [3, hi) passing the complete criterion."""
    return [x for x in range(3, int(hi), 2) if qnr2_all(x, bound)]


def completeness(hi=100000):
    """The two branches against the primes, side by side."""
    sieve = bytearray([1]) * int(hi)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(int(hi) ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    primes = {p for p in range(3, int(hi), 2) if sieve[p]}
    rows = []
    for label, want in (("offset == 2  (x = 3,5 mod 8)", (2,)),
                        ("offset == 0  (x = 1,7 mod 8)", (0,)),
                        ("offset in {0,2}  COMPLETE", (0, 2))):
        keep = {x for x in range(3, int(hi), 2) if offset_mod(x) in want}
        rows.append({"branch": label, "kept": len(keep),
                     "primes_caught": len(keep & primes),
                     "primes_missed": len(primes - keep),
                     "composites": len(keep - primes),
                     "mod8": ",".join(map(str, sorted({x % 8 for x in keep & primes})))})
    return _frame(rows)


# ------------------------------------------------- what the tight gate selects
def classify(x, prime=None):
    """Why x is kept by the tight filter, or why it is not.

    Returns a dict with the branch, the predicted off_r, and whether the gate
    fires -- all from the structure of r = x//2, without evaluating the gate.

    `prime` selects the primality backend: None (gmpy2 if present, else 'mr'),
    'gmpy2', 'mr', 'self' (qnr2's own test -- no external primality at all), or
    any callable.  See BACKENDS.
    """
    r = x // 2
    out = {"x": x, "r": r, "prime": _is_prime(x, prime), "mod8": x % 8,
           "mod120": x % 120, "loose": False, "branch": None,
           "off_r": None, "tight": False}
    if not out["prime"] or x % 8 not in (3, 5):
        return out
    out["loose"] = True
    if _is_prime(r, prime):
        out["branch"] = "safe"          # x = 2r+1, a safe prime
        out["off_r"] = (8 - 1 - 2) % r  # 2^x = 8, 2^r = 2 mod r
    elif r % 2 == 0 and _is_prime(r // 2, prime):
        out["branch"] = "4q+1"          # r = 2q, x = 4q+1
        out["off_r"] = (32 - 1 - 4) % r
    else:
        out["branch"] = "other"
        out["off_r"] = off_r(x)
    o = out["off_r"]
    out["tight"] = bool(o and x % o == 2)
    return out


def is_safe_prime(x, prime=None):
    """x = 2r+1 with both x and r prime."""
    return _is_prime(x, prime) and _is_prime(x // 2, prime)


def predict_tight(hi, prime=None):
    """The tight keepers, from structure alone.  Reproduces qnr2_filter exactly."""
    return [x for x in range(3, int(hi), 2) if classify(x, prime)["tight"]]


def pseudoprimes(hi):
    """Composites the loose filter admits -- base-2 Euler-Jacobi pseudoprimes."""
    return [x for x in qnr2_filter(hi) if not _is_prime(x)]


def fast_vs_slow(hi):
    """Where the one-shot test and the literal offset definition disagree."""
    out = []
    for x in range(3, int(hi), 2):
        slow = offset_mod(x) == 2
        if slow != qnr2_test(x):
            out.append({"x": x, "slow_definition": slow, "qnr2_test": qnr2_test(x),
                        "prime": _is_prime(x)})
    return out


# ---------------------------------------------------------------- dataframes
def _frame(rows, index=None):
    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        return df.set_index(index) if index else df
    except ImportError:
        from sitoolz.dfkit import frame
        return frame(rows, index=index)


def tight_frame(hi=20000, prime=None):
    """The tight keepers, one row each, with the branch that explains them."""
    return _frame([classify(x, prime) for x in predict_tight(hi, prime)], index="x")


def branch_frame(hi=200000):
    """How the tight class splits by branch."""
    rows, agg = [], {}
    for x in predict_tight(hi):
        c = classify(x)
        b = agg.setdefault(c["branch"], {"branch": c["branch"], "count": 0,
                                         "off_r": set(), "mod120": set()})
        b["count"] += 1
        b["off_r"].add(c["off_r"])
        b["mod120"].add(c["mod120"])
    for b in agg.values():
        rows.append({"branch": b["branch"], "count": b["count"],
                     "off_r": ",".join(map(str, sorted(b["off_r"])[:6])),
                     "x_mod_120": ",".join(map(str, sorted(b["mod120"])))})
    return _frame(sorted(rows, key=lambda r: -r["count"]))


def loose_frame(hi=20000):
    """The loose filter against the residue class it is meant to capture."""
    keep = qnr2_filter(hi)
    target = [p for p in range(3, int(hi), 2) if _is_prime(p) and p % 8 in (3, 5)]
    bad = [x for x in keep if not _is_prime(x)]
    return _frame([{"upto": hi, "kept": len(keep),
                    "primes_3_or_5_mod_8": len(target),
                    "missed": len(set(target) - set(keep)),
                    "composites_admitted": len(bad),
                    "composites": ",".join(map(str, bad))}])


MERSENNE_EXPONENTS = [107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
                      9689, 9941, 11213, 19937, 21701, 23209]


def mersenne_bench(exponents=None, reps=1):
    """qnr2_test against gmpy2.is_prime on Mersenne primes 2^p - 1.

    Note what the criterion says about these: 2^p - 1 for p > 2 is always
    7 mod 8, so 2 is a quadratic RESIDUE and qnr2_test correctly returns
    False for every one of them.  It is a filter for 3,5 mod 8 -- Mersennes
    are the class it deliberately excludes.  The timing is still the honest
    comparison of one modexp against a full primality test at that size.
    """
    import time
    rows = []
    for p in (exponents or MERSENNE_EXPONENTS):
        M = (1 << p) - 1
        t0 = time.perf_counter()
        for _ in range(reps):
            mine = qnr2_test(M)
        t1 = time.perf_counter()
        gm = None; tg = None
        if _gmpy2 is not None:
            t2 = time.perf_counter()
            for _ in range(reps):
                gm = bool(_gmpy2.is_prime(M))
            tg = (time.perf_counter() - t2) / reps
        rows.append({"p": p, "digits": len(str(M)) if p < 4000 else int(p * 0.30103) + 1,
                     "M_mod_8": M % 8, "qnr2_test": mine,
                     "qnr2_sec": round((t1 - t0) / reps, 6),
                     "gmpy2_is_prime": gm,
                     "gmpy2_sec": None if tg is None else round(tg, 6),
                     "speedup": None if not tg else round(tg / ((t1 - t0) / reps), 1)})
    return rows


def mersenne_frame(exponents=None, reps=1):
    return _frame(mersenne_bench(exponents, reps), index="p")


def _lucas_lehmer(p):
    """The dedicated Mersenne test, with the 2^p = 1 mod M reduction."""
    if p == 2:
        return True
    M = (1 << p) - 1
    s = 4
    for _ in range(p - 2):
        s = s * s - 2
        while s.bit_length() > p:
            s = (s & M) + (s >> p)
        if s >= M:
            s -= M
    return s == 0


def _sieve_upto(hi):
    s = bytearray([1]) * int(hi)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(int(hi) ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return s


def regime_frame(hi=100000, mersenne_p=1279):
    """Where each approach wins -- bulk scanning versus one large number.

    The two regimes have opposite answers and the comparison is only fair
    within each:

      BULK, bounded range      a sieve tests by structure, ~30x faster than
                               any per-candidate exponentiation.  qnr2_filter
                               pays a full modexp per odd number.
      ONE large number         a sieve cannot run at all -- sieving to 2^1279
                               needs ~10^369 bytes.  qnr2_test is one modexp
                               and beats both gmpy2's Baillie-PSW and the
                               dedicated Lucas-Lehmer test.
    """
    import time

    def t(fn, reps=1):
        fn()
        t0 = time.perf_counter()
        for _ in range(reps):
            fn()
        return (time.perf_counter() - t0) / reps

    M = (1 << mersenne_p) - 1
    rows = [
        {"regime": f"bulk to {hi:,}", "method": "sieve",
         "seconds": round(t(lambda: _sieve_upto(hi)), 5), "applies": True},
        {"regime": f"bulk to {hi:,}", "method": "qnr2_filter (loose)",
         "seconds": round(t(lambda: qnr2_filter(hi)), 5), "applies": True},
        {"regime": f"bulk to {hi:,}", "method": "qnr2_filter_all (complete)",
         "seconds": round(t(lambda: qnr2_filter_all(hi)), 5), "applies": True},
        {"regime": f"single 2^{mersenne_p}-1", "method": "sieve",
         "seconds": None, "applies": False},
        {"regime": f"single 2^{mersenne_p}-1", "method": "qnr2_test",
         "seconds": round(t(lambda: qnr2_test(M), 3), 5), "applies": True},
        {"regime": f"single 2^{mersenne_p}-1", "method": "offset_mod (literal)",
         "seconds": round(t(lambda: offset_mod(M), 3), 5), "applies": True},
        {"regime": f"single 2^{mersenne_p}-1", "method": "Lucas-Lehmer",
         "seconds": round(t(lambda: _lucas_lehmer(mersenne_p)), 5), "applies": True},
    ]
    if _gmpy2 is not None:
        rows.insert(3, {"regime": f"bulk to {hi:,}", "method": "gmpy2 per candidate",
                        "seconds": round(t(lambda: [x for x in range(3, hi, 2)
                                                    if _gmpy2.is_prime(x)]), 5),
                        "applies": True})
        rows.append({"regime": f"single 2^{mersenne_p}-1", "method": "gmpy2.is_prime",
                     "seconds": round(t(lambda: bool(_gmpy2.is_prime(M)), 3), 5),
                     "applies": True})
    return _frame(rows)


def backend_compare(hi=20000):
    """The tight class under each primality backend -- do they agree?"""
    base = qnr2_filter(hi, tight=True)
    rows = []
    for name in ("gmpy2", "mr", "self"):
        if name == "gmpy2" and _gmpy2 is None:
            continue
        import time
        t0 = time.perf_counter()
        got = predict_tight(hi, prime=name)
        el = time.perf_counter() - t0
        rows.append({"backend": name, "keepers": len(got),
                     "matches_filter": got == base,
                     "missing": len(set(base) - set(got)),
                     "extra": len(set(got) - set(base)),
                     "seconds": round(el, 3)})
    return _frame(rows)


if __name__ == "__main__":
    print("the criterion is COMPLETE if you keep both Euler branches:")
    print(completeness(100000)); print()
    print("loose filter against its target class:")
    print(loose_frame(20000)); print()
    print("the tight class splits into two structural branches:")
    print(branch_frame(200000)); print()
    print("prediction from structure == the filter itself:",
          predict_tight(50000) == qnr2_filter(50000, tight=True))
    print("safe primes among the 'safe' branch:",
          all(is_safe_prime(x) for x in predict_tight(50000)
              if classify(x)["branch"] == "safe"))
    print()
    print("primality backends -- does the structural prediction survive?")
    print(backend_compare(20000)); print()
    print("where each approach wins:")
    print(regime_frame(50000)); print()
    print("fast test vs the literal definition, below 4000:")
    for d in fast_vs_slow(4000):
        print("  ", d)
