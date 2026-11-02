"""Roll primes up and down the Pocklington ladder ``N = 2*R*q + 1``.

Descending a prime means factoring ``N-1`` to recover the prime ``q`` it was built
from -- and factoring ``N-1`` is exactly what the base does up to ``REACH``. Rolling
up builds a larger proven prime from a seed; ``prove.build_up`` does the same by construction.
"""
from __future__ import annotations

from math import gcd, isqrt

from ._base import REACH, mod, stream
from .certify import factorize, is_prime


def _odd_part(x: int) -> tuple:
    """Strip all factors of 2: return ``(odd_part, v)`` with ``x = 2**v * odd_part``."""
    x, v = int(x), 0
    while x % 2 == 0:
        x >>= 1
        v += 1
    return x, v


try:                                   # gmpy2's is_prime is BPSW in C -- same test as
    from gmpy2 import is_prime as _gmpy_is_prime   # ours below, ~20-30x faster.
    _HAVE_GMPY2 = True
except Exception:
    _HAVE_GMPY2 = False


def _jacobi(a: int, n: int) -> int:
    """Jacobi symbol ``(a/n)`` for odd ``n > 0``. Returns -1, 0, or 1."""
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):        # (2/n) = -1 when n = 3,5 (mod 8)
                result = -result
        a, n = n, a                     # reciprocity flip
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def _lucas_strong_prp(n: int) -> bool:
    """Strong Lucas probable-prime test, Selfridge parameters -- BPSW's second half.

    Picks the first ``D`` in ``5, -7, 9, -11, ...`` with ``(D/n) = -1``, sets
    ``P=1, Q=(1-D)/4``, then checks ``U_d = 0`` or ``V_{d*2^r} = 0`` for some
    ``r < s``, where ``n+1 = d*2^s``. Pairing this with a base-2 Miller-Rabin is
    Baillie-PSW: no counterexample is known at any size.
    """
    n = int(n)
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    if isqrt(n) ** 2 == n:              # BPSW needs this: squares break D-selection
        return False

    D = 5                                # Selfridge Method A
    while True:
        j = _jacobi(D % n, n)
        if j == -1:
            break
        if j == 0 and abs(D) % n != 0:
            return False                 # D shares a factor with n -> composite
        D = -(D + 2) if D > 0 else -(D - 2)

    P, Q = 1, (1 - D) // 4
    d, s = n + 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    U, V, Qk = 1, P, Q % n               # U_1, V_1, Q^1
    for bit in bin(d)[3:]:               # ladder over d's bits (skip leading 1)
        U, V = (U * V) % n, (V * V - 2 * Qk) % n      # double: U_2k, V_2k
        Qk = (Qk * Qk) % n
        if bit == "1":                                # step: U_2k+1, V_2k+1
            U, V = P * U + V, D * U + P * V
            if U & 1:
                U += n
            if V & 1:
                V += n
            U, V = (U >> 1) % n, (V >> 1) % n
            Qk = (Qk * Q) % n

    if U == 0 or V == 0:
        return True
    for _ in range(1, s):                # V_{d*2^r} == 0 for some r < s ?
        V = (V * V - 2 * Qk) % n
        if V == 0:
            return True
        Qk = (Qk * Qk) % n
    return False


def _probable_prime(n: int) -> bool:
    """Strong probable-prime test -- the gate that stops pseudoprimes.

    Uses gmpy2's **BPSW** when available (one Miller-Rabin base 2 plus one strong
    Lucas test, in C). BPSW has no known counterexample at any size and costs ~2
    tests, versus the 12 Miller-Rabin rounds of the pure-python fallback -- which
    only buys determinism below 3.317e24, a bound that is irrelevant for the
    thousand-digit numbers this library actually chases.

    Without gmpy2 it runs the same BPSW in pure python: :func:`_miller_rabin` base 2
    plus :func:`_lucas_strong_prp`. Two tests instead of twelve, and strictly stronger
    -- 341, 2047, the Carmichaels and every composite Mersenne die on one half or the
    other.
    """
    n = int(n)
    if n < 2:
        return False
    if _HAVE_GMPY2:
        return bool(_gmpy_is_prime(n))
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):   # cheap trial division
        if n % p == 0:
            return n == p
    return _miller_rabin(n, 2) and _lucas_strong_prp(n)      # BPSW


def _miller_rabin(n: int, a: int) -> bool:
    """One Miller-Rabin round, base ``a`` -- the strong-probable-prime test."""
    n = int(n)
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    x = pow(a, d, n)
    if x in (1, n - 1):
        return True
    for _ in range(r - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False


def _is_prime_any(n: int) -> bool:
    """Proven prime up to REACH; Miller-Rabin (deterministic < 3.3e24) beyond."""
    n = int(n)
    if n < 2:
        return False
    return is_prime(n) if n <= REACH else _probable_prime(n)


def descend(N: int) -> list:
    """Roll a prime down: at each rung take the largest prime ``q | (N-1)``.

    Continues while ``q > sqrt(N)`` (the Pocklington condition), returning the chain
    ``[N, q1, q2, ..., seed]``. Each rung ``a`` satisfies ``a = 2*R*q + 1``.

    Factor-based -- needs ``N <= REACH``. For the ``(x-1)/2^v`` halving walk that
    works at any size, use :func:`chain_descend` instead (they are deliberately
    kept separate).
    """
    N = int(N)
    if N > REACH:
        raise ValueError(
            f"N exceeds provable reach {REACH:,}. Use chain_descend(N) for the "
            f"safe-prime halving walk, or factor N-1 with a bigger base / Pollard rho."
        )
    chain, n = [N], N
    while n > 3:
        q = max(factorize(n - 1))
        if q * q <= n:                    # no single large prime factor -> stop
            break
        chain.append(q)
        n = q
    return chain


def rungs(N: int) -> list:
    """The descent as explicit ``(a, R, q)`` triples with ``a = 2*R*q + 1``."""
    chain = descend(N)
    out = []
    for a, q in zip(chain, chain[1:]):
        out.append((a, (a - 1) // (2 * q), q))
    return out


def roll_up(q: int, count: int = 6, rmax: int = 1_000_000) -> list:
    """From prime ``q``, the first ``count`` proven primes ``N = 2*R*q + 1 <= REACH``.

    ``strong`` in each ``(R, N, strong)`` triple means ``2*R < q`` (i.e. ``q > sqrt N``),
    the case a single Pocklington witness proves outright.
    """
    q = int(q)
    if not is_prime(q):
        raise ValueError(f"seed {q} is not prime")
    out = []
    for R in range(1, rmax):
        N = 2 * R * q + 1
        if N > REACH:
            break
        if is_prime(N):
            out.append((R, N, 2 * R < q))
            if len(out) >= count:
                break
    return out


def chain_descend(N: int, target: int = 3) -> list:
    """Walk ``x -> oddpart(x-1)`` down to ``target`` -- the safe-prime / Cunningham
    descent -- verifying each rung as a genuine Pocklington step.

    Every link writes ``x = 2**v * m + 1`` with ``m = oddpart(x-1)``. That link is a
    *real* chain rung only if ``m`` is prime **and** ``m > sqrt(x)`` (``rung=True``);
    otherwise the odd part is composite/small and the value is just a reduction, not
    a member of any prime chain. Works at any size -- primes are proven up to REACH
    and Miller-Rabin (deterministic < 3.3e24) beyond.

    Returns a list of link records ``{x, v, m, rung, in_reach}``. This is distinct
    from :func:`descend`, which factors ``N-1`` for the true largest-prime ladder.
    """
    links, x = [], int(N)
    while x > target:
        m, v = _odd_part(x - 1)
        links.append({
            "x": x,
            "v": v,
            "m": m,
            "rung": (m > 1) and (m * m > x) and _is_prime_any(m),
            "in_reach": x <= REACH,
        })
        if m < 2:
            break
        x = m
    return links


def _pollard_rho(n: int) -> int:
    """A nontrivial factor of composite ``n`` -- last-resort fallback only."""
    n = int(n)
    if n % 2 == 0:
        return 2
    for c in range(1, 40):
        x = y = 2
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n:
            return d
    return n


def _largest_prime_factor_big(m: int) -> int:
    """Largest prime factor of ``m`` at any size.

    Divides ``m`` by our base primes (the smooth part ``R`` of a Pocklington
    ``N-1`` is base-smooth), returning the instant the leftover goes prime -- that
    leftover is the next prime ``q`` in the chain. Pollard rho fires only if a
    genuinely hard cofactor survives (multiple primes above the base), which a
    well-formed Pocklington ``N-1`` never has.

    The base is streamed and only *actual* divisors reach python, so the cost is a
    vectorized ``mod`` per chunk rather than a python step per prime -- the difference
    between seconds and days once the base holds billions.
    """
    m = int(m)
    if m < 2:
        return 1
    big = 1
    for c in stream(isqrt(m)):
        if int(c[0]) ** 2 > m:                    # every remaining prime > sqrt(m)
            break
        for p in c[mod(m, c) == 0]:               # vectorized: only true divisors
            p = int(p)
            while m % p == 0:
                m //= p
            big = p
            if m == 1 or _is_prime_any(m):        # leftover is q -> done, no rho
                return max(big, m)
    if m == 1:
        return big
    if _is_prime_any(m):
        return max(big, m)
    stack, best = [m], big                         # genuine hard cofactor: rho fallback
    while stack:
        t = stack.pop()
        if t == 1:
            continue
        if _is_prime_any(t):
            best = max(best, t)
            continue
        d = _pollard_rho(t)
        stack.append(d)
        stack.append(t // d)
    return best


def deep_descend(N: int) -> list:
    """Descend a prime of **any** size to a seed the base can prove.

    Above REACH each rung is the largest prime of ``x-1`` (peel + Pollard rho,
    primality by Miller-Rabin); at/below REACH the proven :func:`descend` finishes.
    Stops early if some ``x-1`` is too smooth (no prime factor ``> sqrt(x)``).

    Feasible only while ``x-1`` is factorable -- immediate for primes *built* by
    construction (their ``N-1`` is known), generally not for arbitrary huge primes.
    """
    chain, x = [int(N)], int(N)
    while x > REACH:
        q = _largest_prime_factor_big(x - 1)
        if q * q <= x:                    # N-1 too smooth: no Pocklington rung
            break
        chain.append(q)
        x = q
    if x <= REACH:
        chain += descend(x)[1:]           # proven descent the rest of the way
    return chain
