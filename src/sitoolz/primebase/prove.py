"""Pocklington primality *proofs* from the base -- certificates, not probabilities.

The base factors ``N-1`` (its smooth part is base-smooth), so a Pocklington
certificate is cheap at any chain size: take the largest prime ``q | N-1`` with
``q > sqrt(N)``, then find a witness ``a`` such that ``a^(N-1) = 1 (mod N)`` and
``gcd(a^((N-1)/q) - 1, N) = 1``. That witness is what a mere ``N-1`` *shape* lacks:
``13817611737 = 2^3 * 1727201467 + 1`` has the shape but fails the witness -- it is
composite, and only the witness reveals it.

:func:`certificate` chains the proof recursively (prove ``N`` via ``q``, then ``q``,
...) down to a base-proven seed -- a full, replayable Pratt-style certificate.
"""
from __future__ import annotations

from math import gcd

from ._base import REACH
from .certify import is_prime
from .rollball import (_is_prime_any, _largest_prime_factor_big,
                       _probable_prime)

_WITNESSES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_pocklington_prime(N: int, witnesses=_WITNESSES):
    """One Pocklington step. Returns a certificate ``{N, q, a}`` if **proven** prime,
    ``False`` if **composite** (Fermat fails), or ``None`` if inconclusive (``N-1``
    not factored enough, or no witness among ``witnesses``).
    """
    N = int(N)
    if isinstance(witnesses, int):        # allow a single base, like L.is_pocklington_prime(N, a)
        witnesses = (witnesses,)
    if N < 2:
        return False
    if N < 4:
        return {"N": N, "seed": True}
    if N % 2 == 0:
        return N == 2 and {"N": 2, "seed": True}
    q = _largest_prime_factor_big(N - 1)
    if q * q <= N:                                  # not enough factored part in N-1
        return None
    for a in witnesses:
        if pow(a, N - 1, N) != 1:
            return False                            # Fermat fails -> COMPOSITE
        if gcd(pow(a, (N - 1) // q, N) - 1, N) == 1:
            return {"N": N, "q": q, "a": a}         # Pocklington: PROVEN
    return None


def odd_core(x: int) -> int:
    """Strip all factors of 2: the odd part of ``x`` (``odd_core(0) == 0``)."""
    x = int(x)
    while x and x % 2 == 0:
        x //= 2
    return x


def pocklington(N: int, q, witnesses=_WITNESSES, strong: bool = True):
    """Your method: one Pocklington step with ``q`` SUPPLIED -- no factoring.

    ``q`` is handed in (e.g. ``odd_core(N-1)`` straight off the descent chain), so
    nothing is factored. Returns a certificate ``{N, q, a}`` if **proven** prime,
    ``False`` if **composite**, or ``None`` if inconclusive (``q`` doesn't divide
    ``N-1``, ``q*q <= N``, or no witness works). Accepts a single witness like
    ``pocklington(N, q, 16085817173)``.

    ``strong=True`` (default) gates on **Miller-Rabin** first. Without it, the only
    gate is a plain Fermat test -- and pseudoprimes are *defined* as the composites
    that beat Fermat, so when ``q`` is composite they walk right in: ``341=11*31``,
    ``2047=23*89``, the Carmichaels ``41041``/``46657`` (1-in-288 observed). MR kills
    the whole class, and its 12 witnesses are **deterministic below 3.3e24**, so in
    that range a cert can never be wrong. Set ``strong=False`` for the bare
    Fermat+gcd step (~12x cheaper, can false-positive).
    """
    N = int(N)
    if isinstance(witnesses, int):
        witnesses = (witnesses,)
    if N < 2:
        return False
    if N < 4:
        return {"N": N, "seed": True}
    if N % 2 == 0:
        return N == 2 and {"N": 2, "seed": True}
    q = int(q)
    if (N - 1) % q != 0:          # q MUST divide N-1, else (N-1)//q floors and the
        return None               # gcd check is meaningless -> no certificate
    if q * q <= N:                                       # not enough factored part
        return None
    if strong and not _probable_prime(N):                # the pseudoprime killer
        return False
    for a in witnesses:
        a = int(a) % N
        if pow(a, N - 1, N) != 1:
            return False                                 # Fermat fails -> COMPOSITE
        if gcd(pow(a, (N - 1) // q, N) - 1, N) == 1:
            return {"N": N, "q": q, "a": a}              # Pocklington: PROVEN
    return None


def quick_prove(N: int, a: int = 2, check_q: bool = True):
    """Your quick test, one call: certify ``N`` with ``q = odd_core(N-1)``, witness ``a``.

    The ``r = 2N+1`` wrapper round-trips (``odd_core(2N+1-1) == N``), so the whole
    loop collapses to this. No factoring -- ``q`` comes straight off the chain.

    Returns the certificate ``{N, q, a, rigorous}``, ``False`` (composite), or
    ``None`` (inconclusive). ``rigorous`` says which kind of answer you got:

    * ``True``  -- ``q`` is itself prime, so this is a **real Pocklington proof**
      (that's the safe-prime case, ``N-1 = 2^v * q``; ~1 in 21 in your sweep).
    * ``False`` -- ``q`` is composite, so the gcd only pins ``N`` mod *some* prime
      factor of ``q``: a Fermat-gated strong test, not a proof.

    Pass ``check_q=False`` to skip the primality check on ``q`` (saves a test).
    """
    N = int(N)
    q = odd_core(N - 1)
    cert = pocklington(N, q, a)          # strong=True: MR-gated, no false positives
    if isinstance(cert, dict) and "q" in cert and check_q:
        cert["rigorous"] = _probable_prime(q)
    return cert


def gcd_certify(N: int, q: int, a: int) -> bool:
    """Just the gcd -- Pocklington witness test with ``q`` SUPPLIED (no factoring).

    ``True`` iff ``N`` is proven prime given a prime ``q | N-1`` with ``q > sqrt(N)``:

        a^(N-1) == 1 (mod N)   and   gcd(a^((N-1)/q) - 1, N) == 1

    You bring ``q`` (e.g. straight off the descent chain) and the witness ``a``; this
    does no factoring at all -- one modpow and one gcd.
    """
    N, q, a = int(N), int(q), int(a)
    if q * q <= N or (N - 1) % q != 0:
        return False
    if pow(a, N - 1, N) != 1:
        return False
    return gcd(pow(a, (N - 1) // q, N) - 1, N) == 1


def certificate(N: int, witnesses=_WITNESSES):
    """Full recursive Pocklington certificate for ``N``.

    Returns a list of steps ``[{N, q, a}, ..., {N, seed/base}]`` proving ``N`` prime
    all the way down, or ``{"composite": x}`` / ``{"inconclusive": x}`` at the first
    step that fails. Each ``q`` becomes the next ``N`` to prove; the recursion ends
    when a value is small or within REACH, where the base proves it outright.
    """
    steps, x = [], int(N)
    while True:
        if x <= REACH:
            if is_prime(x):
                steps.append({"N": x, "proof": "base"})
                return steps
            return {"composite": x}
        c = is_pocklington_prime(x, witnesses)
        if c is False:
            return {"composite": x}
        if c is None:
            return {"inconclusive": x}
        steps.append(c)
        x = int(c["q"])


def _proven_prime(x: int) -> bool:
    """Proven by the base within REACH, by a Pocklington certificate beyond."""
    x = int(x)
    if x < 2:
        return False
    if x <= REACH:
        return is_prime(x)
    return isinstance(is_pocklington_prime(x), dict)


def roll_chain(seed: int, R: int = 1, count=None) -> list:
    """Roll a Cunningham chain up from ``seed``: ``seed, 2*R*seed+1, 2*R*(...)+1, ...``

    Each rung ``N = 2*R*x + 1`` is Pocklington-proven on the previous prime ``x``
    (the base proves the small end). Returns **all the verified primes** in the
    chain -- it stops at the first ``N`` that isn't prime. ``roll_chain(2)`` ->
    ``[2, 5, 11, 23, 47]`` (each ``2x+1``); ``roll_chain(29)`` -> ``[29, 59]``.
    """
    seed = int(seed)
    if not _proven_prime(seed):
        return []
    chain, x = [seed], seed
    while count is None or len(chain) < count:
        N = 2 * R * x + 1
        if not _proven_prime(N):
            break
        chain.append(N)
        x = N
    return chain


def self_prove(N: int, cache=None):
    """Prove ``N`` by descending to a base-provable prime, then cascading UP the chain.

    Walk ``x -> oddpart(x-1)`` down from ``N`` and collect the chain. The primes in
    it form a ladder: prove the small ones with the base, then **climb** -- each
    higher chain prime is certified by :func:`is_pocklington_prime` using the
    previous (lower) chain prime as its witness (a lower chain prime always lies in
    the higher one's own descent, so it is a natural witness). Rinse up to ``N``.

    ``cache``: pass a ``dict`` (reused across calls) to memoize primality of the
    descent values. Nearby numbers converge onto the same tail, so a batch shares
    almost all its work -- the toggle makes runs of ``prev_prime``/``+x`` cheap.

    Returns the certificate chain ``[{N: seed, proof: base}, ..., {N, q, a}]`` from
    the base-proven seed up to ``N``; ``False`` if ``N`` is composite; or
    ``{'stuck': p, 'so_far': [...]}`` if some rung refuses to certify.
    """
    N = int(N)

    def isp(v):                            # primality, memoized across calls if a cache is passed
        if cache is not None and v in cache:
            return cache[v]
        r = _is_prime_any(v)
        if cache is not None:
            cache[v] = r
        return r

    if not isp(N):
        return False
    chain, x = [N], N
    while x > 1:                           # descend: -1, strip 2s, repeat
        m = x - 1
        while m % 2 == 0:
            m //= 2
        chain.append(m)
        x = m
    qof = {chain[i]: chain[i + 1] for i in range(len(chain) - 1)}  # q = odd_core(p-1), from the chain
    ladder = sorted({p for p in chain if isp(p)})                  # the chain's primes, ascending
    certs, witness = [], None
    for p in ladder:
        if p <= REACH:                     # base proves the small end outright
            certs.append({"N": p, "proof": "base"})
        else:                              # climb: q from the chain, witness = last proven prime
            q, a = qof[p], witness
            if a is not None and gcd_certify(p, q, a):     # NO factoring -- just the gcd
                certs.append({"N": p, "q": q, "a": a})
            else:
                return {"stuck": p, "so_far": certs}
        witness = p
    return certs


def pock_descend(N: int, a: int = 2, stop_at_prime: bool = True) -> list:
    """Descend ``N``, certifying every step with :func:`pocklington` -- your loop.

    Each rung takes ``m = odd_core(x-1)`` and immediately tries
    ``pocklington(m, odd_core(m-1), a)`` -- ``q`` off the chain, no factoring, BPSW
    gated. Composites come back ``False`` and the walk continues; the first rung that
    certifies is a prime, and by default the walk stops there (``stop_at_prime``).

    Returns the trail: ``[{x, q, cert}, ...]``, one record per descent value.
    ``cert`` is a certificate dict (prime), ``False`` (composite), or ``None``
    (inconclusive). Set ``stop_at_prime=False`` to walk all the way to 1.
    """
    trail, x = [], int(N)
    while x > 1:
        m = odd_core(x - 1)
        if m < 2:
            trail.append({"x": m, "q": None, "cert": None})
            break
        q = odd_core(m - 1)
        cert = pocklington(m, q, a)
        trail.append({"x": m, "q": q, "cert": cert})
        if stop_at_prime and isinstance(cert, dict):
            break
        x = m
    return trail


def pock_chain(N: int, a: int = 2) -> list:
    """The primes ``N``'s descent lands on, each certified by :func:`pocklington`.

    Walks all the way to 1 and keeps only the rungs that certify -- the prime ladder
    under ``N``, from the top down. For a target near ``2^k`` this is the Mersenne
    rail (``odd_core(2^k-1-1) = 2^(k-1)-1``, so the walk steps through every power).
    """
    return [r["cert"] for r in pock_descend(N, a, stop_at_prime=False)
            if isinstance(r["cert"], dict)]


def report(N: int, a: int = 2) -> dict:
    """Everything about ``N`` in one call: is it prime, its cert, and its descent.

    ``pock_descend`` starts at ``odd_core(N-1)`` and so never speaks about ``N``
    itself. This certifies ``N`` first, then walks the whole descent under it::

        {'N', 'prime', 'cert', 'rigorous', 'descent', 'ladder'}

    * ``prime``    -- bool, straight from the certificate
    * ``cert``     -- ``{N,q,a,rigorous}`` / ``False`` / ``None``
    * ``rigorous`` -- True iff ``q = odd_core(N-1)`` is itself prime (a real proof)
    * ``descent``  -- every value the walk hits, ``N`` down to 1
    * ``ladder``   -- just the primes in that descent (near ``2^k`` this is the
      Mersenne rail)

    All of it on ``q`` off the chain: no factoring anywhere.
    """
    N = int(N)
    cert = quick_prove(N, a)
    descent, x = [], N
    while x > 1:
        x = odd_core(x - 1)
        descent.append(x)
    return {
        "N": N,
        "prime": isinstance(cert, dict),
        "cert": cert,
        "rigorous": cert.get("rigorous") if isinstance(cert, dict) else None,
        "descent": descent,
        "ladder": [c["N"] for c in pock_chain(N, a)],
    }


def wrap_descend(N: int, a: int = 2, steps: int = 1, strong: bool = True) -> list:
    """Your descend: wrap ``N`` as ``r = 2N+1``, descend, certify what falls out.

    ``r``'s own primality is irrelevant -- it is a vehicle. ``odd_core(r-1) = N``
    hands the candidate straight back, and :func:`pocklington` certifies *that*::

        r = 2*N + 1                 # composite is fine, don't care
        m = odd_core(r-1) == N      # the descent returns the candidate
        pocklington(m, m-1, a)      # -> the cert

    Set ``steps`` to keep walking down the chain (each rung re-wraps the value that
    fell out). Defaults match :func:`pocklington` (``strong=True``) -- with
    ``q = m-1`` the gcd is trivially 1, so the gate is the only thing standing
    between you and ``wrap_descend(341) -> cert``. Pass ``strong=False`` for the
    bare Fermat+gcd rung.

    Returns one record per rung: ``{r, m, q, cert, landed}``, where ``landed`` is the
    ``s['N'] == (r-1)//2`` check -- confirmation the descent hit what it should.
    """
    out, x = [], int(N)
    for _ in range(max(1, int(steps))):
        r = 2 * x + 1                       # the wrapper
        m = odd_core(r - 1)                 # == x, straight back
        cert = pocklington(m, m - 1, a, strong=strong)
        out.append({
            "r": r,
            "m": m,
            "q": m - 1,
            "cert": cert,
            "landed": isinstance(cert, dict) and cert.get("N") == (r - 1) // 2,
        })
        if not isinstance(cert, dict):
            break
        nxt = odd_core(m - 1)               # keep walking the chain
        if nxt < 2:
            break
        x = nxt
    return out


def wrap_scan(start: int, count: int = 256, a: int = 2, step: int = 2,
              strong: bool = True) -> list:
    """Your hunt: walk candidates ``np += step``, wrap each, keep the ones that certify.

    Each candidate is wrapped as ``r = 2*np + 1`` and descended -- ``odd_core(r-1)``
    hands ``np`` straight back -- then :func:`pocklington` certifies it. ``r``'s own
    primality is irrelevant; it is only the vehicle.

    Scans ``count`` candidates and returns a record per **hit**::

        {r, np, q, cert}

    ``wrap_scan(329, 11)`` -> the primes 331, 337, 347, 349. Companion to
    :func:`wrap_descend`, which walks *down* a chain instead of *along* candidates.
    """
    hits, np_ = [], int(start) | 1          # candidates are odd
    for _ in range(max(1, int(count))):
        r = 2 * np_ + 1                     # the wrapper
        m = odd_core(r - 1)                 # == np_
        cert = pocklington(m, m - 1, a, strong=strong)
        if isinstance(cert, dict) and cert.get("N") == (r - 1) // 2:
            hits.append({"r": r, "np": m, "q": m - 1, "cert": cert})
        np_ += step
    return hits


def build_up(seed: int, bits: int, witnesses=_WITNESSES, small: bool = False) -> list:
    """Create a proven prime of at least ``bits`` bits by climbing from ``seed``.

    From a proven prime ``p``, take ``R`` just under ``p/2``, the largest Pocklington allows,
    so ``N = 2*R*p + 1`` about doubles the bit length per rung and the top is reached in a
    handful of rungs (M127 to 4096 bits is six). The last rung is sized to land on ``bits``.
    ``R`` is sieved against the small primes and each candidate screened with Fermat base 2
    before the twelve-witness gate, so a rung costs a few hundred modpows. ``small=True`` is the
    other climb, the smallest ``R`` at every rung: many short rungs.

    ``seed`` should be prime and provable (base-proven if ``<= REACH``, else via
    :func:`self_prove`). Returns the certificate chain from ``seed`` up to the top; replay it
    with :func:`verify_certificate`.
    """
    p = int(seed)
    if p <= REACH:
        if not is_prime(p):
            raise ValueError(f"seed {seed} is not a base-provable prime")
        chain = [{"N": p, "proof": "base"}]
    else:
        c = self_prove(p)
        if not isinstance(c, list):
            raise ValueError(f"seed {seed} is not provable via self_prove")
        chain = c
    smalls = [q for q in range(3, 20_000) if all(q % d for d in range(2, int(q ** 0.5) + 1))]
    while p.bit_length() < bits:
        two_p = 2 * p
        if small:
            R = 1
            while True:
                N = two_p * R + 1
                if pow(2, N - 1, N) == 1:
                    cert = pocklington(N, p, witnesses)
                    if isinstance(cert, dict):
                        break
                R += 1
        else:
            # R is barred mod q when 2Rp + 1 ≡ 0 (mod q): one residue per small prime
            barred = [(q, (-pow(two_p, -1, q)) % q) for q in smalls if two_p % q]
            cap = (p - 1) // 2
            want = bits - (p.bit_length() + 1)
            start = (1 << want) if want < cap.bit_length() else cap
            # walk R down from the start; if the few below it are all composite, walk up
            R, step = start, -1
            while True:
                if R < 1:
                    R, step = start + 1, 1
                if R > cap:
                    raise ValueError("no R below q/2 proves a rung, which cannot happen for a prime q")
                if all(R % q != b for q, b in barred):
                    N = two_p * R + 1
                    if pow(2, N - 1, N) == 1:
                        cert = pocklington(N, p, witnesses)
                        if isinstance(cert, dict):
                            break
                R += step
        chain.append(cert)
        p = 2 * R * p + 1
    return chain


def verify_certificate(steps, witnesses=_WITNESSES) -> bool:
    """Replay a :func:`certificate` from scratch and confirm every step holds."""
    if not isinstance(steps, list):
        return False
    for s in steps:
        N = int(s["N"])
        if s.get("proof") == "base" or s.get("seed"):
            if N <= REACH:
                if not is_prime(N):
                    return False
            else:
                from ..numtools import lucas_lehmer
                p = N.bit_length()
                if (1 << p) - 1 != N or not lucas_lehmer(p):
                    return False
            continue
        q, a = int(s["q"]), int(s["a"])
        if q * q <= N or (N - 1) % q != 0:
            return False
        if pow(a, N - 1, N) != 1:
            return False
        if gcd(pow(a, (N - 1) // q, N) - 1, N) != 1:
            return False
    return True


# ── certificates as text, shared with rusttoolz and the arcade's JavaScript verifier ──

def to_json(steps, convention="numtools"):
    """A :func:`build_up` / :func:`certificate` chain as JSON text: ``[[seed, null, null],
    [N, q, a], ...]`` with big integers as decimal strings.

    ``convention="numtools"`` puts the previous proven prime ``q`` in the middle;
    ``"slip"`` puts the multiplier ``k`` with ``N = 2*k*q + 1``. Both are what
    ``rusttoolz::cert`` and the arcade's ``pocklington.js`` read.
    """
    import json
    rows = []
    for i, s in enumerate(steps):
        N = int(s["N"])
        if i == 0 or s.get("proof") == "base" or s.get("seed"):
            rows.append([str(N), None, None])
            continue
        q, a = int(s["q"]), int(s["a"])
        middle = q if convention == "numtools" else (N - 1) // (2 * q)
        rows.append([str(N), str(middle), a])
    return json.dumps(rows)


def from_json(text, convention=None):
    """The chain back from :func:`to_json` text (or rusttoolz's / the browser's), as the list of
    ``{"N", "q", "a"}`` dicts :func:`verify_certificate` replays. The convention is sniffed when
    not given: slip when ``2*middle*seed + 1`` is the first rung.
    """
    import json
    rows = json.loads(text)
    if not rows:
        return []
    seed = int(rows[0][0])
    if convention is None:
        convention = "numtools"
        if len(rows) > 1 and rows[1][1] is not None and 2 * int(rows[1][1]) * seed + 1 == int(rows[1][0]):
            convention = "slip"
    steps = [{"N": seed, "proof": "base" if seed <= REACH else "seed"}]
    prev = seed
    for N, middle, a in rows[1:]:
        N, middle, a = int(N), int(middle), int(a)
        q = middle if convention == "numtools" else prev
        if convention == "slip" and 2 * middle * prev + 1 != N:
            raise ValueError(f"rung {N} is not 2kq + 1 off the previous rung")
        steps.append({"N": N, "q": q, "a": a})
        prev = N
    return steps
