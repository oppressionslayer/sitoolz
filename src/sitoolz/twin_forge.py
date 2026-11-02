"""twin_forge.py -- hunt provable twin primes of the form  k*2^n +- 1,  in parallel.

Own the fast side. `k*2^n + 1` is a Proth prime: N-1 = k*2^n is factored by construction, so
one Pocklington witness proves it -- fast and rigorous. Every certified +1 is a twin CANDIDATE;
the open half is "is k*2^n - 1 also prime?" (the Riesel / Lucas side, N+1 = k*2^n factored).

This sieves k hard (both sides), screens survivors with BPSW, and on a hit proves BOTH sides
outright (Pocklington for +1, Lucas/BLS for -1), streaming each as a hand-off packet.

Self-contained: gmpy2 optional (fast BPSW screen). The proofs are our own -- no oracle.

    from sitoolz.twin_forge import hunt, estimate
    estimate(256)                       # how big is this search?
    hunt(256, out="twins_256.jsonl")    # find + prove, stream packets
"""
from __future__ import annotations
import math, os, json, time
from math import gcd, isqrt
from multiprocessing import Pool, cpu_count

try:
    import gmpy2
    def _bpsw(x): return bool(gmpy2.is_prime(gmpy2.mpz(x)))
    HAVE_GMPY2 = True
except ImportError:                                     # pragma: no cover
    HAVE_GMPY2 = False
    def _bpsw(n, bases=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)):
        n = int(n)
        if n < 2: return False
        for p in bases:
            if n % p == 0: return n == p
        d, r = n - 1, 0
        while d % 2 == 0: d //= 2; r += 1
        for a in bases:
            x = pow(a, d, n)
            if x in (1, n - 1): continue
            for _ in range(r - 1):
                x = x * x % n
                if x == n - 1: break
            else: return False
        return True

_WIT = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)


# ---------------- sieve ----------------
def _oddprimes(bound):
    bound = int(bound)
    s = bytearray([1]) * (bound + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, isqrt(bound) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(3, bound + 1) if s[i]]     # odd primes only


def _bad_residues(n, bound):
    """k-residues mod each odd prime p<=bound that kill either side: k ≡ ± inv(2^n) (mod p)."""
    out = []
    for p in _oddprimes(bound):
        t = pow(2, n, p)
        if t == 0:
            continue
        r = pow(t, -1, p)
        out.append((p, r % p, (-r) % p))
    return out


def _sieve(n, k_lo, k_hi, bad):
    import numpy as np
    alive = np.ones(k_hi - k_lo, dtype=bool)
    for p, r1, r2 in bad:
        alive[(r1 - k_lo) % p::p] = False
        alive[(r2 - k_lo) % p::p] = False
    return (k_lo + np.flatnonzero(alive)).tolist()


# ---------------- proofs (our own) ----------------
def _factor_small(k):
    fp, d = {2}, 3
    while d * d <= k:
        while k % d == 0: fp.add(d); k //= d
        d += 2
    if k > 1: fp.add(k)
    return sorted(fp)


def _pocklington(N, n, fps):                            # +1 side: N-1 = k*2^n
    for p in fps:
        for a in _WIT:
            if pow(a, N - 1, N) != 1:
                return None
            if gcd(pow(a, (N - 1) // p, N) - 1, N) == 1:
                break
        else:
            return None
    return True                                         # N-1 = k*2^n fully factored (> sqrt N)


def _jacobi(a, m):
    a %= m; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if m % 8 in (3, 5): r = -r
        a, m = m, a
        if a % 4 == 3 and m % 4 == 3: r = -r
        a %= m
    return r if m == 1 else 0


def _lucas_uv(m, P, Q, mod):
    D = P * P - 4 * Q; inv2 = pow(2, -1, mod)
    U, V, Qk = 0, 2, 1
    for b in bin(m)[2:]:
        U, V = (U * V) % mod, (V * V - 2 * Qk) % mod; Qk = (Qk * Qk) % mod
        if b == "1":
            U, V = ((P * U + V) * inv2) % mod, ((D * U + P * V) * inv2) % mod; Qk = (Qk * Q) % mod
    return U, V


def _lucas(N, fps):                                     # -1 side: N+1 = k*2^n
    if N < 3 or N % 2 == 0 or isqrt(N) ** 2 == N:
        return False
    Ds, D = [], 5
    while len(Ds) < 60:
        j = _jacobi(D, N)
        if j == 0 and abs(D) % N != 0: return False
        if j == -1: Ds.append(D)
        D = -D + 2 if D < 0 else -D - 2
    for D in Ds:
        for P in range(1, 40, 2):
            Q = (P * P - D) // 4
            if (P * P - D) % 4 or gcd(N, 2 * Q * D) != 1: continue
            if _lucas_uv(N + 1, P, Q, N)[0] != 0: continue
            if all(gcd(_lucas_uv((N + 1) // q, P, Q, N)[0], N) == 1 for q in fps):
                return True
    return False


# ---------------- worker: test one candidate (the fast side first) ----------------
def _test_one(args):
    n, k = args
    base = 1 << n
    if not _bpsw(k * base + 1):                         # +1 Proth side first (fast)
        return (k, False, False)
    return (k, True, _bpsw(k * base - 1))               # (k, +1 prime, -1 prime?)


# ---------------- estimate ----------------
def estimate(n, sieve_bound=100_000):
    lnN = n * math.log(2)
    digits = int(n * math.log10(2)) + 1
    C2 = 0.6601618
    k_per_plus = lnN / 2                                 # ~k tried per +1 Proth prime
    k_per_twin = lnN * lnN / (2 * C2)                    # ~k tried per twin
    print(f"n={n}  ->  {digits}-digit numbers,  ln N = {lnN:.0f}")
    print(f"  ~1 in {k_per_plus:,.0f} sieved k gives a +1 Proth prime")
    print(f"  ~1 in {k_per_twin:,.0f} sieved k gives a TWIN (both sides prime)")
    return dict(digits=digits, k_per_plus=k_per_plus, k_per_twin=k_per_twin)


# ---------------- driver ----------------
def hunt(n, k_start=1, workers=None, sieve_bound=100_000, span=2_000_000,
         want_twins=10**9, out=None, checkpoint=None, verbose=True):
    """Search k for k*2^n +- 1 twins, streaming a line per candidate.

    Sieves a span of k (fast), then tests survivors' +1 side across workers; on a +1 Proth
    prime it proves it (Pocklington) and writes a hand-off packet to `out`, and if the -1
    side is prime too it proves that (Lucas) -- a full twin. Checkpointed and resumable.
    Runs until `want_twins` full twins are found (default: effectively forever, so it keeps
    streaming +1 candidates). Returns the list of proven twins.
    """
    workers = workers or max(1, cpu_count() - 1)
    base = 1 << n
    bad = _bad_residues(n, sieve_bound)
    k = k_start
    if checkpoint and os.path.exists(checkpoint):
        k = json.load(open(checkpoint)).get("next_k", k_start)
    found, plus, tested, t0 = [], 0, 0, time.time()
    print(f"[hunt] n={n} ({int(n*math.log10(2))+1} digits), {workers} workers, "
          f"start k={k:,}", flush=True)
    with Pool(workers) as pool:
        while len(found) < want_twins:
            survivors = _sieve(n, k, k + span, bad)
            print(f"[sieve] k in [{k:,}, {k+span:,}) -> {len(survivors):,} survivors to test",
                  flush=True)
            for kk, plus_prime, minus_prime in pool.imap_unordered(
                    _test_one, [(n, s) for s in survivors]):
                tested += 1
                if plus_prime:
                    plus += 1
                    hi, lo = kk * base + 1, kk * base - 1
                    fps = _factor_small(kk)
                    up = _pocklington(hi, n, fps) is not None
                    packet = {"k": kk, "n": n, "digits": len(str(hi)),
                              "upper": str(hi), "upper_proven_pocklington": up,
                              "lower": str(lo),
                              "lower_status": "probable_prime" if minus_prime else "composite"}
                    if minus_prime:
                        packet["lower_proven_lucas"] = _lucas(lo, fps)
                        found.append(packet)
                    if out:
                        with open(out, "a") as f:
                            f.write(json.dumps(packet) + "\n")
                    tag = "*** TWIN ***" if minus_prime else "+1 Proth prime (twin candidate)"
                    print(f"  [HIT] k={kk}  {packet['digits']} digits  {tag}  +1proven={up}",
                          flush=True)
                if verbose and tested % 4 == 0:
                    dt = (time.time() - t0) / 60
                    print(f"[hunt] tested {tested:,} | +1 primes {plus} | twins {len(found)} "
                          f"| {dt:.1f} min | {tested/max(dt,1e-9):.1f}/min", flush=True)
            k += span
            if checkpoint:
                json.dump({"next_k": k}, open(checkpoint, "w"))
    return found
