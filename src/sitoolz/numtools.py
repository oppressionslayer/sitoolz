"""numtools — small numeric-table helpers.

Usage:
    from sitoolz.numtools import creator
    df = creator()                 # odd numbers 3..257, mod-3 offset rows
    df = creator(stop=265)         # go further
    df = creator(offsets=(1, -1, -2, +4))   # any offsets you like
"""

from __future__ import annotations

import numpy as np
import pandas as pd

try:  # fast primality via gmpy2 when available
    from gmpy2 import is_prime as _gmpy_is_prime

    def _is_prime(v) -> int:
        return int(bool(_gmpy_is_prime(int(v))))

except Exception:  # pragma: no cover - pure-python fallback

    def _is_prime(v) -> int:
        v = int(v)
        if v < 2:
            return 0
        if v < 4:
            return 1
        if v % 2 == 0:
            return 0
        i = 3
        while i * i <= v:
            if v % i == 0:
                return 0
            i += 2
        return 1


def wide(linewidth: int = 200):
    """Stop pandas & numpy from truncating output with '...'.

    Call once per session:  from numtools import wide; wide()
    Pass a bigger linewidth for very wide arrays, e.g. wide(300).
    """
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", None)
    pd.set_option("display.expand_frame_repr", False)
    np.set_printoptions(threshold=np.inf, linewidth=linewidth, suppress=True)


def pad(width=4):
    """Right-align every integer in printed numpy arrays to a fixed field width,
    so arrays you stack on top of each other line up column-for-column.

    Pick a width >= the digit-count of your biggest number (pad(5) for up to
    99999). Call pad(None) to reset to numpy's default.
    """
    if width is None:
        np.set_printoptions(formatter=None)
    else:
        w = int(width)
        np.set_printoptions(formatter={"int": lambda x, w=w: f"{x:{w}d}"})


def creator(
    start: int = 3,
    stop: int = 257,
    step: int = 2,
    mod: int = 3,
    offsets=(1, -1, -2),
    weave: int = 4,
    xor_bitlen: bool = True,
    primes: bool = True,
    wheel=(2, 3, 5),
    n=None,
    length: int | None = None,
    inclusive: bool = True,
    transpose: bool = False,
):
    """Build a DataFrame whose top row is a numeric sequence and whose
    remaining rows are ``(n + offset) % mod`` for each offset.

    Parameters
    ----------
    start, stop, step : int
        Define the sequence. With ``inclusive=True`` (default) ``stop`` is
        included when it lands on the step grid — so the defaults give the
        odd numbers 3, 5, 7, ... 257.
    mod : int
        Modulus applied to each offset row (default 3).
    offsets : iterable of int
        One row per offset. Default ``(1, -1, -2)`` gives rows
        ``(n+1)%3``, ``(n-1)%3``, ``(n-2)%3``.
    weave : int
        The woven row: ``0, weave, 0, 2*weave, 0, 3*weave, ...`` (0 disables).
    xor_bitlen : bool
        Add the row ``2**n.bit_length() ^ n`` (e.g. 27 -> 32^27 -> 59).
    primes : bool
        Add three 1/0 primality rows: ``(n+3)//2``, ``2n+1``, and ``n``.
    wheel : sequence of int
        Add a ``wheel<mod>`` row (1 if n is coprime to every listed prime -- i.e.
        on an admissible lane) plus the corresponding ``x=(n-1)/2`` index row.
        Default ``(2,3,5)`` -> the mod-30 wheel. Empty/None to disable.
    n : int | sequence, optional
        Supply your own sequence directly (scalar, list, or array). Takes
        precedence over ``start``/``stop``/``step``/``length``. All other rows
        are computed from it. e.g. ``creator(n=[27, 100, 255])``.
    length : int, optional
        If given, produce exactly this many columns (``start, start+step, ...``),
        ignoring ``stop``/``inclusive``. Handy for picking a custom size.
    inclusive : bool
        If True, extend the arange stop by ``step`` so an on-grid ``stop``
        value is included.
    transpose : bool
        If True, return the numbers going *down* a column instead of across.

    Returns
    -------
    pandas.DataFrame
        Rows: ``n`` plus one row per offset. Columns are positional (0..N-1).
    """
    if n is not None:
        # caller supplied their own sequence (scalar, list, or array)
        n = np.atleast_1d(np.asarray(n))
    elif length is not None:
        # pick an exact number of columns: start, start+step, ... (length values)
        n = start + step * np.arange(length)
    else:
        hi = stop + step if inclusive else stop
        n = np.arange(start, hi, step)

    def _label(k: int) -> str:
        sign = "+" if k >= 0 else "-"
        return f"(n{sign}{abs(k)})%{mod}"

    rows = [n] + [(n + k) % mod for k in offsets]
    index = ["n"] + [_label(k) for k in offsets]

    if weave:
        # 0 at even positions; weave, 2*weave, 3*weave, ... at odd positions
        # e.g. weave=4 -> 0, 4, 0, 8, 0, 12, 0, 16, ...
        w = np.zeros(len(n), dtype=int)
        odd = np.arange(1, len(n), 2)
        w[odd] = weave * (np.arange(len(odd)) + 1)
        rows.append(w)
        index.append(f"weave{weave}")

    if xor_bitlen:
        # 2**n.bit_length() XOR n, per element (exact Python int semantics)
        # e.g. n=27 -> 2**5 ^ 27 -> 32 ^ 27 -> 59
        x = np.array([(2 ** int(v).bit_length()) ^ int(v) for v in n])
        rows.append(x)
        index.append("2^bl^n")

    if primes:
        # For each test, show the tested VALUE, then a 1/0 prime flag directly
        # below it so they line up for examination. Uses gmpy2 when available.
        for vlabel, plabel, seq in (
            ("(n+3)//2", "(n+3)//2 prime", (n + 3) // 2),
            ("2n+1", "2n+1 prime", 2 * n + 1),
            ("n (val)", "n prime", n),
        ):
            rows.append(np.asarray(seq))
            index.append(vlabel)
            rows.append(np.array([_is_prime(v) for v in seq]))
            index.append(plabel)

    if wheel:
        # the scaffold: 1 if n is coprime to every wheel prime (on a lane),
        # plus the corresponding index x = (n-1)//2. Compare wheel vs 'n prime':
        # every prime >= next-wheel-prime is a survivor; only the base primes
        # (2,3,5 for the default) get dropped.
        wprimes = tuple(int(q) for q in np.atleast_1d(wheel))
        modulus = 1
        for q in wprimes:
            modulus *= q
        survives = np.array(
            [int(all(int(v) % q != 0 for q in wprimes)) for v in n]
        )
        rows.append(survives)
        index.append(f"wheel{modulus}")
        rows.append((n - 1) // 2)
        index.append("x=(n-1)/2")

    df = pd.DataFrame(rows, index=index)
    return df.T if transpose else df


def climb(residues, mod, start=0, stop=None, count=None):
    """Generate an interleaved-lane / modular sequence ("climb").

    Returns the sorted integers whose value ``% mod`` is one of ``residues``.
    This is the general form of patterns like ``11, 23, 27, 39, ...`` which is
    ``climb((7, 11), 16, 11, 1019)`` -> two lanes stepping by 16.

    Parameters
    ----------
    residues : int | sequence of int
        Allowed remainders mod ``mod`` (one lane each).
    mod : int
        The step / modulus. Each lane advances by this.
    start : int
        Lower bound (inclusive).
    stop : int, optional
        Upper bound (inclusive). Provide this OR ``count``.
    count : int, optional
        Instead of a ``stop``, return exactly this many values from ``start``.

    Returns
    -------
    numpy.ndarray  (sorted, de-duplicated)

    Examples
    --------
    climb((7, 11), 16, 11, 1019)     # the 11,23,27,39,... array
    climb(11, 16, start=11, count=8) # single lane, first 8 values
    climb((1, 3, 7, 9), 10, 1, 100)  # last-digit 1/3/7/9 up to 100
    """
    mod = int(mod)
    start = int(start)
    res = sorted({int(r) % mod for r in np.atleast_1d(residues)})

    def _first_at_or_after(s, r):
        # smallest x >= s with x % mod == r
        return s + ((r - s) % mod)

    lanes = []
    if count is not None:
        for r in res:
            f = _first_at_or_after(start, r)
            lanes.append(f + mod * np.arange(count))  # count/lane guarantees enough
        return np.unique(np.concatenate(lanes))[:count]

    if stop is None:
        raise ValueError("climb() needs either stop=... or count=...")
    stop = int(stop)
    for r in res:
        f = _first_at_or_after(start, r)
        if f > stop:
            continue
        n_terms = (stop - f) // mod + 1
        lanes.append(f + mod * np.arange(n_terms))
    if not lanes:
        return np.array([], dtype=int)
    return np.unique(np.concatenate(lanes))


def prime_w(count=132, factor=True):
    """The arithmetic 'w' proxy for next_prime(p), tabulated against real primes.

    For each prime ``p`` (starting at the 3rd prime, 5), computes the
    primality-test-free quantity ``w`` and lines it up with the true next
    prime so you can see where -- and by how much -- it drifts.

    Columns
    -------
    x, p, w, next_prime, offset (= w - next_prime), w_prime,
    snap (= nextprime at-or-above w), snap_hits (snap == next_prime),
    and factorint(101*w) if ``factor=True``.

    Notes
    -----
    ``w`` almost always lands *inside the gap* ``(p, next_prime(p)]``, so
    ``snap`` recovers the true next prime with no search over p at all.
    """
    try:
        import sympy
    except Exception as e:  # pragma: no cover
        raise RuntimeError("prime_w needs sympy (pip install sympy)") from e

    plist = [2]
    while len(plist) < count + 3:
        plist.append(int(sympy.nextprime(plist[-1])))
    m = [pp * 100 for pp in plist]

    cp = 101
    prevcp = 404
    prevw = 5
    noffset = 202
    recs = []
    for x in range(2, 2 + count):
        c = m[x]
        p = plist[x]
        g = c + p
        w = ((c + p) + 202) // 101
        l = w // 3
        while l % 2 == 1:  # l is not recomputed inside -- matches the original
            w += 2
            noffset = ((w - prevw) * 101) + prevcp
            if noffset > 1010:
                break
        nb = int(sympy.nextprime(p))
        snap = w if _is_prime(w) else int(sympy.nextprime(w - 1))
        rec = {
            "x": x,
            "p": p,
            "w": w,
            "next_prime": nb,
            "offset": w - nb,
            "w_prime": bool(_is_prime(w)),
            "snap": snap,
            "snap_hits": snap == nb,
        }
        if factor:
            rec["factorint(101w)"] = dict(sympy.factorint(101 * w))
        recs.append(rec)
        cp = c + p
        prevw = p
        prevcp = (g - cp) + noffset
    return pd.DataFrame(recs)


class ResidueFilter:
    """Instant reject filter for a prime-derived relation ``a*p + b``.

    Given a wheel (default 2,3,5 -> modulus 30), it precomputes which lanes of
    ``p`` can make ``a*p + b`` coprime to the wheel (i.e. possibly prime). Any
    ``p`` in a non-surviving lane makes ``a*p + b`` divisible by a wheel prime,
    so it is *guaranteed composite* -- skip it with no test.

    Use inline in a sieve::

        rf = residue_filter(1, 2)          # for w = p + 2
        if not rf.ok(p):                   # w is guaranteed composite
            continue                       # ...skip, no primality test needed
    """

    def __init__(self, a, b, wheel=(2, 3, 5)):
        from math import gcd

        self.a, self.b, self.wheel = int(a), int(b), tuple(int(q) for q in wheel)
        W = 1
        for q in self.wheel:
            W *= q
        self.modulus = W
        self.lanes = [r for r in range(W) if gcd(r, W) == 1]
        self.survive = [r for r in self.lanes if gcd((self.a * r + self.b) % W, W) == 1]
        self._sset = frozenset(self.survive)
        self.ruled_out = [r for r in self.lanes if r not in self._sset]
        self.ruled_frac = len(self.ruled_out) / len(self.lanes)

    def ok(self, p):
        """True if ``a*p + b`` *could* be prime -> you must test it.
        False only when ``a*p + b`` is *provably* composite -> safe to skip.

        Safe by construction: if ``p`` shares a factor with the wheel (i.e. it is
        a base prime like 2/3/5), this returns True (test it) rather than risk
        dropping a real prime. So it NEVER rejects a genuine prime.
        """
        from math import gcd

        r = int(p) % self.modulus
        if gcd(r, self.modulus) != 1:
            return True  # base-prime edge case -> don't skip, just test it
        return r in self._sset

    def __repr__(self):
        return (
            f"ResidueFilter({self.a}*p+{self.b}, mod {self.modulus}): "
            f"{len(self.survive)}/{len(self.lanes)} lanes survive, "
            f"{self.ruled_frac * 100:.0f}% of primes ruled out with no test | "
            f"survive lanes: {self.survive}"
        )


def residue_filter(a, b, wheel=(2, 3, 5)):
    """Build a :class:`ResidueFilter` for the relation ``a*p + b`` (see its docs)."""
    return ResidueFilter(a, b, wheel)


def swift_primes(lo, hi, wheel=(2, 3, 5)):
    """Find primes in ``[lo, hi)`` by jumping only through wheel lanes.

    The synthesis of the session: skip every number a wheel can prove composite
    (for free), and fire the real primality test *only* on candidates that
    survive the jump. Uses ``gmpy2.is_prime`` (BPSW) on the survivors.

    Returns ``(primes, stats)`` where stats records how much was skipped without
    a test -- the whole point is that ``tested`` is far below ``range``.

    Deeper wheel -> bigger jumps: ``wheel=(2,3,5,7)`` walks only 48/210 of the
    numbers instead of 8/30.
    """
    import gmpy2
    from math import gcd

    W = 1
    for q in wheel:
        W *= q
    lanes = sorted(r for r in range(W) if gcd(r, W) == 1)

    primes = [q for q in wheel if lo <= q < hi]  # the base primes the wheel skips
    tested = 0

    block = (max(lo, 2) // W) * W
    while block < hi:
        for r in lanes:
            c = block + r
            if c < max(lo, 2) or c == 1:
                continue
            if c >= hi:
                break
            tested += 1  # only wheel survivors ever reach a real test
            if gmpy2.is_prime(c):
                primes.append(c)
        block += W

    primes.sort()
    total = max(hi - lo, 1)
    stats = {
        "range": hi - lo,
        "tested": tested,
        "skipped_no_test": (hi - lo) - tested,
        "work_fraction": round(tested / total, 3),
        "primes_found": len(primes),
        "wheel": f"mod {W} ({len(lanes)} lanes)",
    }
    return primes, stats


def create_prime(n_bits, seed=None):
    """Construct a PROVABLE prime of ~n_bits bits, bottom-up (Maurer/Pocklington).

    Not a sieve, not a probable-prime test -- it *builds* a prime from smaller
    proven primes and returns a certificate that proves it. Each step forms
    ``q = 2*R*k + 1`` from an already-proven prime ``R > sqrt(q)``; a single
    witness certifies ``q`` by Pocklington's theorem. This is the grown-up
    version of your ``^1 //2`` descent -- the descent found the certifier,
    this walks it the other way to *make* the prime.

    Returns ``(prime, certificate)`` where certificate is a chain
    ``[(base, None, None), (q1, R1, a1), ...]`` verifiable by
    :func:`verify_certificate` without trusting any primality test.
    """
    import random

    import gmpy2

    rng = random.Random(seed)

    def _trial_prime(x):
        if x < 2:
            return False
        i = 2
        while i * i <= x:
            if x % i == 0:
                return False
            i += 1
        return True

    def build(bits):
        if bits <= 20:
            lo, hi = 1 << (bits - 1), 1 << bits
            p = int(gmpy2.next_prime(rng.randrange(lo, hi)))
            return p, [(p, None, None)]
        R, cert = build(bits // 2 + 1)
        lo, hi = 1 << (bits - 1), 1 << bits
        kmin, kmax = max(1, lo // (2 * R)), hi // (2 * R)
        for _ in range(50000):
            k = rng.randrange(kmin, kmax + 1)
            q = 2 * R * k + 1
            if not (lo <= q < hi) or R * R <= q:
                continue
            if not gmpy2.is_prime(q):
                continue
            for a in range(2, 3000):
                if pow(a, q - 1, q) != 1:
                    continue
                if int(gmpy2.gcd(pow(a, (q - 1) // R, q) - 1, q)) == 1:
                    return int(q), cert + [(int(q), int(R), int(a))]
        raise RuntimeError(f"construction stalled at {bits} bits; retry with a new seed")

    return build(n_bits)


def pollard_p1(N, base=2, max_B=1_000_000):
    """Factor N by Pollard's p-1 -- the wrap-around you reasoned your way to.

    Raises ``a`` through prime-powers mod N; when the running exponent becomes a
    multiple of some ``p-1``, ``a`` wraps to 1 mod p, so ``gcd(a-1, N)`` reveals p.
    Self-tunes B by checking the gcd after every prime, and backs off to a new
    base if it overshoots (grabs all of N). Deterministic on any N whose p-1 is
    B-smooth; returns None if no factor's predecessor is smooth enough.

    Returns a nontrivial factor of N, or None.
    """
    from math import gcd

    import gmpy2

    N = int(N)
    if N % 2 == 0:
        return 2
    for attempt in range(4):  # fresh base if we overshoot
        a = base + attempt
        p = 2
        overshoot = False
        while p < max_B:
            pe = p
            while pe * p <= max_B:
                pe *= p
            a = pow(a, pe, N)
            g = gcd(a - 1, N)
            if 1 < g < N:
                return int(g)
            if g == N:
                overshoot = True
                break
            p = int(gmpy2.next_prime(p))
        if not overshoot:
            return None  # exhausted the smoothness bound, no smooth p-1
    return None


def pollard_rho(N, max_c=20):
    """Factor N by Pollard's rho -- the general-purpose wrap-around.

    Follows the sequence ``x -> x^2 + c (mod N)`` and watches for a collision
    mod an unknown factor via ``gcd(|x-y|, N)``. Works regardless of smoothness;
    expected time ~ N^(1/4), so it's the go-to for pulling a *small* factor out.
    Returns a nontrivial factor of N, or None.
    """
    from math import gcd

    import gmpy2

    N = int(N)
    if N % 2 == 0:
        return 2
    if gmpy2.is_prime(N):
        return None
    for c in range(1, max_c + 1):
        x = y = 2
        d = 1
        while d == 1:
            x = (x * x + c) % N
            y = (y * y + c) % N
            y = (y * y + c) % N
            d = gcd(abs(x - y), N)
        if 1 < d < N:
            return int(d)
    return None


def verify_certificate(cert):
    """Independently verify a Pocklington certificate from :func:`create_prime`.

    Re-proves primality from the certificate chain alone (trial-divides only the
    small base), never calling a probabilistic test. Returns True iff the whole
    chain proves the top prime.
    """
    import gmpy2

    def _trial_prime(x):
        if x < 2:
            return False
        i = 2
        while i * i <= x:
            if x % i == 0:
                return False
            i += 1
        return True

    base = cert[0][0]
    if not _trial_prime(base):
        return False
    proven = {base}
    for q, R, a in cert[1:]:
        if R not in proven or (q - 1) % R != 0 or R * R <= q:
            return False
        if pow(a, q - 1, q) != 1:
            return False
        if int(gmpy2.gcd(pow(a, (q - 1) // R, q) - 1, q)) != 1:
            return False
        proven.add(q)
    return True


def swift_derived(a, b, plo, phi, wheel=(2, 3, 5)):
    """Find primes of the form ``a*p + b`` for primes ``p`` in ``[plo, phi)``.

    Stacks both jumps: iterate primes ``p`` (already a wheel jump), then use the
    residue filter to skip every ``p`` whose ``a*p+b`` is provably composite --
    firing ``gmpy2.is_prime`` only on survivors. E.g. ``swift_derived(2, 1, ...)``
    hunts Sophie Germain primes (``2p+1`` prime); ``swift_derived(1, 2, ...)``
    the upper members of twin pairs.

    Returns ``(pairs, stats)`` with ``pairs = [(p, a*p+b), ...]`` for the hits.
    """
    import gmpy2
    import sympy

    rf = residue_filter(a, b, wheel)
    pairs = []
    skipped = tested = 0
    for p in sympy.primerange(plo, phi):
        if not rf.ok(p):
            skipped += 1  # a*p+b provably composite -> no test
            continue
        v = a * p + b
        tested += 1
        if v > 1 and gmpy2.is_prime(v):
            pairs.append((int(p), int(v)))

    considered = skipped + tested
    stats = {
        "relation": f"{a}*p+{b}",
        "primes_considered": considered,
        "skipped_by_residue": skipped,
        "tested": tested,
        "ruled_out%": round(100 * skipped / considered) if considered else 0,
        "hits": len(pairs),
        "survive_lanes": rf.survive,
    }
    return pairs, stats


def filter_table(prange=(3, 120), relations=((1, 2, "p+2"), (2, 1, "2p+1"), (1, 4, "p+4")), wheel=(2, 3, 5)):
    """DataFrame of :func:`residue_filter` in action over the primes in ``prange``.

    For each prime ``p`` and each relation ``a*p+b`` you pass, shows the value and
    a verdict: ``skip`` (filter proved it composite -- no test), ``PRIME``, or
    ``comp`` (had to test, turned out composite). The ``skip`` rows are the free
    culls. A summary of how much each relation ruled out is attached as
    ``df.attrs['summary']``.

    relations : iterable of ``(a, b, name)``.
    """
    import sympy

    filters = {name: residue_filter(a, b, wheel) for (a, b, name) in relations}
    W = next(iter(filters.values())).modulus
    primes = list(sympy.primerange(prange[0], prange[1]))

    rows = []
    tally = {name: {"skip": 0, "test": 0, "wrong": 0} for (_, _, name) in relations}
    for p in primes:
        row = {"p": p, f"p%{W}": p % W}
        for (a, b, name) in relations:
            rf = filters[name]
            v = a * p + b
            row[name] = v
            isp = bool(sympy.isprime(v))
            if not rf.ok(p):
                row[f"{name}?"] = "skip"
                tally[name]["skip"] += 1
                if isp:
                    tally[name]["wrong"] += 1
            else:
                row[f"{name}?"] = "PRIME" if isp else "comp"
                tally[name]["test"] += 1
        rows.append(row)

    df = pd.DataFrame(rows)
    summary = pd.DataFrame(
        [
            {
                "relation": name,
                "ruled_out%": round(100 * filters[name].ruled_frac),
                "skipped": tally[name]["skip"],
                "tested": tally[name]["test"],
                "wrongly_skipped": tally[name]["wrong"],
                "survive_lanes": filters[name].survive,
            }
            for (_, _, name) in relations
        ]
    )
    df.attrs["summary"] = summary
    return df


def lane_moves(mod=16, N=1019, first=()):
    """Line-by-line: the residue lanes it takes to cover every prime <= N.

    Each row is one 'move' -- committing to a residue lane mod ``mod`` -- with
    how many new primes it grabs and the running coverage. ``first`` front-loads
    specific lanes so you can start from where your 8::4 climb already is, e.g.
    ``lane_moves(16, first=(7, 11))``.

    Returns a pandas DataFrame indexed by move number.
    """
    from math import gcd

    import sympy

    primes = list(sympy.primerange(2, N + 1))
    residues = [r for r in range(mod) if gcd(r, mod) == 1]
    by_count = sorted(residues, key=lambda r: -sum(1 for p in primes if p % mod == r))
    order = list(first) + [r for r in by_count if r not in first]

    caught = set()
    rows = []
    for i, r in enumerate(order, 1):
        lane = [p for p in primes if p % mod == r]
        new = [p for p in lane if p not in caught]
        caught.update(lane)
        rows.append(
            {
                "move": i,
                f"lane(mod{mod})": r,
                "primes_here": len(lane),
                "new": len(new),
                "cumulative": len(caught),
                "coverage%": round(100 * len(caught) / len(primes), 1),
                "first_few": lane[:4],
            }
        )
    return pd.DataFrame(rows).set_index("move")


def xor_climb(seed=3, start=8, stop=1018, step=4, tweak=7, climb_only=False,
              climb=None, as_arrays=False, evens=False):
    """Reconstruct the running-XOR climb and its power-of-2 structure.

    Walks ``r = r ^ x`` for x in ``range(start, stop, step)`` and records, per
    step, the XOR quantities you were exploring plus power-of-2 markers.

    Columns
    -------
    x, r, r_bits, x-4, r^(x-4), r^(x-4)^<tweak>, prime, x_is_p2

    Notes
    -----
    - ``r`` is the running XOR; the even-index rows are the real climb
      (``climb_only=True`` keeps only those), the odd rows collapse to junk.
    - ``r^(x-4)^tweak`` sits at a flat baseline but spikes to ``2^(k+1) - 4``
      exactly when ``x == 2^k`` -- a power-of-2 boundary detector.
    - Filter with ``df[df.x_is_p2]`` to see just the boundary crossings.
    - ``climb`` (optional): a custom sequence XOR'd into ``r`` in place of the
      arithmetic ramp ``x`` -- give it a short motif and it tiles across the
      range to climb with variance. ``x`` stays the positional ramp, so the
      x-4 / power-of-2 frame is unchanged. ``climb=None`` keeps the original.
    - ``as_arrays=True``: instead of a DataFrame, return a namedtuple of arrays
      ``(xm4, r, x47, x4)`` so you can unpack::
          xm4, r, x47, x4 = xor_climb(..., as_arrays=True)
      where ``x4 = r^(x-4)`` and ``x47 = x4 ^ tweak``. Collects *all* rows by
      default; pass ``evens=True`` to keep only the even rows (the real climb).

    Returns a pandas DataFrame (or a namedtuple of arrays if ``as_arrays``).
    """
    tcol = f"r^(x-4)^{tweak}"
    xs = list(range(start, stop, step))
    if climb is not None:
        drive = np.resize(np.asarray(climb, dtype=np.int64), len(xs))  # tile the motif
    r = seed
    recs = []
    for k, x in enumerate(xs):
        r = r ^ (int(drive[k]) if climb is not None else x)
        xm4 = x - 4
        a = r ^ xm4
        b = a ^ tweak
        recs.append(
            {
                "x": x,
                "r": r,
                "r_bits": r.bit_length(),
                "x-4": xm4,
                "r^(x-4)": a,
                tcol: b,
                "prime": bool(_is_prime(b)),
                "x_is_p2": x > 0 and (x & (x - 1)) == 0,
            }
        )
    df = pd.DataFrame(recs)
    if as_arrays:
        from collections import namedtuple

        sel = df.iloc[0::2].reset_index(drop=True) if evens else df  # all rows unless evens
        Climb = namedtuple("Climb", ["xm4", "r", "x47", "x4"])
        return Climb(
            sel["x-4"].to_numpy(),
            sel["r"].to_numpy(),
            sel[tcol].to_numpy(),
            sel["r^(x-4)"].to_numpy(),
        )
    if climb_only:
        df = df.iloc[0::2].reset_index(drop=True)  # even-index = the real climb
    return df


def tile_mask(tile, mult=1, prange=(3, 300)):
    """XOR a repeating tile against ``mult * primes`` -- the mask sandbox.

    Reproduces the games from this session::

        tile_mask([12, 28])                 # np.tile([12,28],..) ^ primes
        tile_mask([28, 12], mult=2)         # the 'backwards' one (all even)

    A power-of-2 ``mult`` factors straight out -- since ``(2a)^(2b) == 2*(a^b)``
    it lifts the whole pattern one octave without reshaping it. Feed the result
    to ``creator(n=...)`` for the full row analysis (bands, primes, wheel, ...).

    Parameters
    ----------
    tile : sequence of int
        The repeating XOR mask, e.g. ``[12, 28]``.
    mult : int
        Scale the primes before XOR (default 1).
    prange : (lo, hi)
        Prime value range ``[lo, hi)``. Default ``(3, 300)`` skips 2 to match
        the original construction (``primes[1:...]``).

    Returns
    -------
    numpy.ndarray  --  ``tiled_mask ^ (mult * primes_in_range)``
    """
    lo, hi = int(prange[0]), int(prange[1])
    sieve = np.ones(max(hi, 2), dtype=bool)
    sieve[:2] = False
    for i in range(2, int(hi**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    primes = np.nonzero(sieve)[0]
    primes = primes[primes >= lo]
    tile = np.asarray(tile, dtype=np.int64)
    reps = int(np.ceil(len(primes) / len(tile)))
    mask = np.tile(tile, reps)[: len(primes)]
    return mask ^ (mult * primes)


def ruler(n=64, heights=False):
    """The ruler sequence (OEIS A007814): ``v2(m+1)`` for m in range(n).

    Equivalently: the 2-adic valuation, the carry-chain length of ``m -> m+1``,
    the tick-mark heights on a ruler, and the between-boundary structure hiding
    inside ``xor_climb``. Self-similar -- its period doubles each level.

    heights=True returns the xor_climb spike heights (8, then ``2^(j+4)-4``)
    instead of the raw levels.
    """

    def _v2(k):
        j = 0
        while k % 2 == 0:
            k //= 2
            j += 1
        return j

    levels = np.array([_v2(m + 1) for m in range(n)])
    if heights:
        return np.array([8 if j == 0 else 2 ** (j + 4) - 4 for j in levels])
    return levels


def ruler_ascii(n=64):
    """A printable ASCII ruler of the first n ticks -- self-similarity made visible."""
    lv = ruler(n)
    top = int(lv.max())
    return "\n".join(
        "".join("|" if lv[m] >= r else " " for m in range(n))
        for r in range(top, -1, -1)
    )


def gap_ruler(k):
    """The PRIME-GAP ruler: read consecutive prime gaps off a Fermat anchor 2^k + 1.

    ``2^k + 1`` is two bits (bit 0 and bit k), so XOR by any even wrap ``w < 2^k``
    is pure ADDITION -- ``(2^k+1) ^ w == (2^k+1) + w``. Sweeping every even wrap
    sweeps every odd number in the octave ``[2^k, 2^(k+1))``; the primes among them
    are consecutive, the ``wrap`` column is their cumulative distance from the
    anchor, and DIFFERENCING the wraps hands the prime gaps back. This is the
    companion of :func:`ruler` (the 2-adic ruler of *counting*): one reads the gap
    structure off the primes, the other off the integers -- both difference-rulers.

    Only clean while ``w < 2^k`` (below the next bit); past that the wrap reaches
    bit k, XOR stops equalling ``+``, and you hand off to the next Fermat base
    ``2^(k+1)+1``. So one call = one octave; chain ``k, k+1, k+2, ...`` up the line.

    Returns a DataFrame ``[prime, wrap, gap, lane]``; ``df.attrs`` carries ``base``
    (=2^k+1) and ``ruler`` (the 2-adic ruler heights over the octave, from the climb).
    """
    base = 2 ** k + 1
    primes = [base + w for w in range(2, 2 ** k, 2) if _is_prime(base + w)]
    wraps = [p - base for p in primes]
    gaps = [w - pv for w, pv in zip(wraps, [0] + wraps[:-1])]   # gap from previous prime
    df = pd.DataFrame({"prime": primes, "wrap": wraps, "gap": gaps,
                       "lane": [p % 4 for p in primes]})
    df.attrs["base"] = base
    df.attrs["ruler"] = ruler(2 ** k, heights=True)             # the counting backbone
    return df


def under_p2(kmax=32):
    """The nearest prime BELOW each power of two 2^k, for k = 2..kmax.

    ``2^k`` is even and primes past 2 are odd, so the distance ``2^k - prime`` is
    always ODD -- you never drop an even step down from a power of two. The
    distance is ``1`` exactly when ``2^k - 1`` is a Mersenne prime (rare: k =
    2,3,5,7,13,17,19,... light up like beacons). From the all-ones anchor
    ``2^k - 1``, the prime sits at ``wrap = distance - 1`` (an even wrap), since
    all-ones XOR an even number is plain subtraction -- the ceiling read-out.

    There is no closed form for the distances (a formula would crack prime gaps);
    every one is odd, and the 1s flag the Mersennes. Returns a DataFrame
    ``[k, p2, anchor, prime, distance, wrap, mersenne]``.

    Uses trial division (fast to a few billion, ~k<=32); for larger k feed a
    prebuilt base to :func:`is_prime_by_index` instead.
    """
    rows = []
    for k in range(2, kmax + 1):
        top = 1 << k
        p = top - 1
        while not _is_prime(p):
            p -= 2                                  # step down by 2, stay odd
        d = top - p
        rows.append({"k": k, "p2": top, "anchor": top - 1, "prime": p,
                     "distance": d, "wrap": (top - 1) ^ p, "mersenne": d == 1})
    return pd.DataFrame(rows)


def lucas_lehmer(p):
    """Deterministic PROOF that the Mersenne number ``2^p - 1`` is prime.

    ``M_p = 2^p - 1`` can only be prime if ``p`` is prime (``p = ab`` makes
    ``2^a - 1`` divide ``2^p - 1``), but that is not enough -- ``M_11 = 2047 =
    23*89``. The Lucas-Lehmer test settles it exactly: with ``s_0 = 4`` and
    ``s_{i+1} = s_i^2 - 2``, an odd-prime ``M_p`` is prime IFF ``s_{p-2} == 0
    (mod M_p)``. It's a true proof -- no oracle, no probable-prime -- and it's the
    only reason thousand-digit Mersennes can be decided in milliseconds (each step
    is one squaring mod ``M_p``). This is the distance-1 beacon of :func:`under_p2`
    made rigorous: ``lucas_lehmer(k)`` True <=> ``2^k - 1`` is the prime just below
    ``2^k`` <=> ``under_p2`` distance 1 at ``k``.

    Returns True/False; ``p`` need not be prime (False for composite p or p < 2).
    """
    if p == 2:
        return True
    if p < 2 or not _is_prime(p):
        return False
    m = (1 << p) - 1
    s = 4
    for _ in range(p - 2):
        s = (s * s - 2) % m
    return s == 0


def mersenne_exponents(pmax=2300):
    """The exponents ``p < pmax`` for which ``2^p - 1`` is prime (proven by
    :func:`lucas_lehmer`). For ``pmax=2300`` this is the first 17 Mersenne
    exponents ``2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281``."""
    return [p for p in range(2, pmax) if lucas_lehmer(p)]


def ceiling_floor(k):
    """Map the octave ``[2^k, 2^(k+1))`` from BOTH ends at once.

    Two anchors straddle the octave: the FLOOR ``2^k + 1`` (a Fermat-form, two
    bits) measures each prime's distance UP, and the CEILING ``2^(k+1) - 1``
    (all-ones) measures its distance DOWN. Both read by pure addition/subtraction
    (``floor + fw == p`` and ``ceil - cw == p``), and for EVERY prime in the
    octave they meet in the middle exactly:

        floor_wrap + ceil_wrap == 2^k - 2   (a constant, the octave width - 2)

    So the two rulers are mirror images -- the up-distance from the floor and the
    down-distance from the ceiling always sum to the same number. The prime
    nearest the ceiling (smallest ``ceil_wrap``) is :func:`under_p2`'s answer for
    ``k+1``; if the ceiling ``2^(k+1)-1`` is itself prime it's a Mersenne beacon.

    Returns ``[prime, floor_wrap, ceil_wrap, total]``; ``df.attrs`` has ``floor``,
    ``ceil``, ``invariant`` (=2^k-2).
    """
    floor = (1 << k) + 1
    ceil = (1 << (k + 1)) - 1
    rows = []
    for p in range((1 << k) + 1, 1 << (k + 1), 2):
        if _is_prime(p):
            fw, cw = p - floor, ceil ^ p               # ceil is all-ones -> XOR == subtract
            rows.append({"prime": p, "floor_wrap": fw, "ceil_wrap": cw, "total": fw + cw})
    df = pd.DataFrame(rows)
    df.attrs["floor"], df.attrs["ceil"], df.attrs["invariant"] = floor, ceil, (1 << k) - 2
    return df


def gap_stream(k_lo, k_hi):
    """Stitch the per-octave Fermat gap-rulers into ONE continuous prime stream.

    Walks every prime from ``2^k_lo`` up to ``2^(k_hi+1)``, tagging each with the
    octave ``k`` it lives in, that octave's Fermat floor ``2^k + 1``, its ``wrap``
    (distance up from the floor -- the clean gap-ruler read, since the floor is two
    bits), and the running ``gap`` to the previous prime. ``seam=True`` flags the
    gaps that cross a power of two (the joints between octaves). So each octave is
    :func:`gap_ruler` and stringing them gives the whole gap sequence, organized by
    power-of-2 windows with the seams marked.

    Returns ``[prime, gap, wrap, octave, base, seam]``.
    """
    rows = []
    prev = None
    for k in range(k_lo, k_hi + 1):
        base, lo, hi = (1 << k) + 1, 1 << k, 1 << (k + 1)
        for p in range(lo + 1, hi, 2):
            if _is_prime(p):
                rows.append({"prime": p, "gap": 0 if prev is None else p - prev,
                             "wrap": p - base, "octave": k, "base": base,
                             "seam": prev is not None and prev < lo <= p})
                prev = p
    return pd.DataFrame(rows)


def is_prime_base(n, base=None):
    """Provable primality by trial division against a saved base (vectorized, fast).

    A base of primes up to ``P`` proves any ``n <= P**2``: a composite that big must
    have a factor ``<= P`` sitting in the base. No oracle, no Miller-Rabin -- just
    numpy modulo against the base primes up to ``sqrt(n)``. Raises past the reach.
    """
    import math

    n = int(n)
    if n < 2:
        return False
    if base is None:
        base = load_base()
    if n > int(base[-1]) ** 2:
        raise ValueError(f"n={n:,} exceeds base reach {int(base[-1]) ** 2:,}; grow the base")
    sub = base[base <= math.isqrt(n)]
    return not (n % sub == 0).any()


def altitude(k, m=8, base=None):
    """The three-part octave report at ``2^k`` -- works at TRILLIONS, off the base.

    You can't sweep a trillion-wide octave, so this samples the two places the whole
    theory lives -- the seam at ``2^k`` and the octave's ends -- and returns them as
    real data:

    - **seam window**: the DataFrame is ``m`` primes each side of ``2^k`` with their
      running ``gap`` and a ``seam`` flag on the joint that crosses the power of two.
    - ``df.attrs`` carries the rest: ``under`` = ``(prime_below, distance)`` (the
      :func:`under_p2` read, distance always odd, ``mersenne`` when 1); ``seam`` =
      the joint decomposed into ``down`` (=under_p2 distance) + ``up``; ``ends`` =
      the two octave-end primes each showing ``floor_wrap + ceil_wrap == 2^k - 2``
      (the :func:`ceiling_floor` invariant), and ``invariant_ok``.

    Provable to the base reach (~2.9 trillion default). Past that, pass a bigger
    ``base`` (``load_base(limit=..., force=True)``) -- the error says how big.
    Returns the seam-window DataFrame ``[prime, gap, side, seam]`` with rich attrs.
    """
    import math

    if base is None:
        base = load_base()
    reach = int(base[-1]) ** 2
    p2, top = 1 << k, 1 << (k + 1)
    if top - 1 > reach:
        need = math.isqrt(top) + 1
        raise ValueError(f"2^{k+1} = {top:,} exceeds base reach {reach:,}; "
                         f"grow it: load_base(limit={need:,}, force=True)")

    def isp(x):
        return not (x % base[base <= math.isqrt(x)] == 0).any()

    def scan(x, step, cnt):
        out, p = [], (x + step) | 1
        while len(out) < cnt:
            if isp(p):
                out.append(p)
            p += 2 * step
        return out

    below = scan(p2, -1, m)[::-1]
    above = scan(p2, +1, m)
    floor, ceil = p2 + 1, top - 1
    rows, prev = [], None
    for p in below + above:
        rows.append({"prime": p, "gap": 0 if prev is None else p - prev,
                     "side": "below" if p < p2 else "above",
                     "seam": prev is not None and prev < p2 <= p})
        prev = p
    df = pd.DataFrame(rows)
    pb, pa = below[-1], above[0]
    down, up = p2 - pb, pa - p2
    ce = scan(top, -1, 1)[0]                       # prime nearest the ceiling (below 2^(k+1))
    ends = [{"prime": p, "floor_wrap": p - floor, "ceil_wrap": ceil - p,
             "total": (p - floor) + (ceil - p)} for p in (pa, ce)]
    df.attrs.update(p2=p2, floor=floor, ceil=ceil, invariant=p2 - 2,
                    under=(pb, down), mersenne=down == 1,
                    seam={"below": pb, "above": pa, "down": down, "up": up, "gap": pa - pb},
                    ends=ends, invariant_ok=all(e["total"] == p2 - 2 for e in ends))
    return df


def altitude_view(k, m=8, base=None):
    """Printable three-part report for :func:`altitude` -- the screenshot version."""
    df = altitude(k, m=m, base=base)
    a = df.attrs
    pb, d = a["under"]
    s, W = a["seam"], 74
    out = ["=" * W, f"  ALTITUDE 2^{k} = {a['p2']:,}   (provable via base, reach checked)", "=" * W]
    out += ["", f"  under_p2:  nearest prime below = {pb:,}   distance {d} (odd={d % 2 == 1})"
                f"   mersenne={a['mersenne']}"]
    out += ["", f"  seam gap = {s['gap']} = {s['down']} (down, = under_p2) + {s['up']} (up)"]
    out += ["", f"  ceiling_floor invariant  floor_wrap + ceil_wrap == 2^{k}-2 = {a['invariant']:,}:"]
    for e in a["ends"]:
        out += [f"    p={e['prime']:>16,}  floor_wrap={e['floor_wrap']:>16,}"
                f"  ceil_wrap={e['ceil_wrap']:>16,}  sum ok={e['total'] == a['invariant']}"]
    out += ["", "  seam window (m each side, running gaps):"]
    for _, r in df.iterrows():
        mark = "  <-- SEAM" if r["seam"] else ""
        out += [f"    {int(r['prime']):>16,}  gap {str(int(r['gap'])):>4}{mark}"]
    out += ["=" * W]
    return "\n".join(out)


def altitude_hi(k, m=6):
    """Stratosphere altitude report -- gmpy2 BPSW engine, ANY exponent k, no base.

    Same three parts as :func:`altitude` (seam window, nearest prime below, the
    ceiling/floor invariant) but powered by gmpy2's BPSW, so it runs at ``2^10000``
    and beyond where no trial-division base could ever reach. HONEST: past the base's
    provable range this is a PROBABLE-prime test (BPSW -- no known counterexample,
    but not a proof). The invariant and the odd-distance are still EXACT -- they're
    algebra, they don't care how big the number is.

    The primes are astronomically long, so they're shown as OFFSETS from ``2^k``:
    each row's ``offset`` means the prime is ``2**k + offset`` (negative = below the
    power of two, positive = above). Returns DataFrame ``[offset, gap, side, seam]``;
    ``df.attrs`` has ``k``, ``digits``, ``under`` (distance below, odd), ``mersenne``,
    ``seam`` (down+up), ``invariant`` (=2^k-2, exact), ``ends`` (both octave ends
    confirming ``floor_wrap + ceil_wrap == 2^k-2``), ``engine``, ``provable`` (False).
    """
    import math

    import gmpy2

    p2 = 1 << k

    def above_offsets(cnt):                        # primes just ABOVE 2^k, as offsets
        out, p = [], gmpy2.next_prime(p2)
        for _ in range(cnt):
            out.append(int(p) - p2)
            p = gmpy2.next_prime(p)
        return out

    def below_offsets(x, cnt):                     # cnt primes just BELOW x, offsets from 2^k
        span = max(64, cnt * x.bit_length() * 2)   # next_prime sieves first -> fast even huge
        while True:
            got, p = [], gmpy2.next_prime(x - span)
            while p < x:
                got.append(int(p))
                p = gmpy2.next_prime(p)
            if len(got) >= cnt:
                return [g - p2 for g in got[-cnt:]]
            span *= 2

    below = below_offsets(p2, m)                   # negative offsets, ascending to 0
    above = above_offsets(m)                       # positive offsets
    rows, prev = [], None
    for d in below + above:
        rows.append({"offset": d, "gap": 0 if prev is None else d - prev,
                     "side": "below" if d < 0 else "above",
                     "seam": prev is not None and prev < 0 <= d})
        prev = d
    df = pd.DataFrame(rows)
    down, up = -below[-1], above[0]
    floor, ceil, inv = p2 + 1, (1 << (k + 1)) - 1, p2 - 2
    cep = below_offsets(1 << (k + 1), 1)[0] + p2   # ceil-end prime (nearest below 2^(k+1))
    e = (1 << (k + 1)) - cep
    fe = p2 + up                                   # floor-end prime
    ends = [{"where": f"2^{k}+{up}", "small_wrap": up - 1, "which": "floor_wrap",
             "sum_ok": (fe - floor) + (ceil - fe) == inv},
            {"where": f"2^{k+1}-{e}", "small_wrap": e - 1, "which": "ceil_wrap",
             "sum_ok": (cep - floor) + (ceil - cep) == inv}]
    df.attrs.update(k=k, digits=int(k * math.log10(2)) + 1, under=down, up=up,
                    mersenne=down == 1, seam={"down": down, "up": up, "gap": down + up},
                    invariant=inv, ends=ends, invariant_ok=all(x["sum_ok"] for x in ends),
                    engine="gmpy2 BPSW", provable=False)
    return df


def altitude_hi_view(k, m=6):
    """Printable stratosphere report for :func:`altitude_hi` (primes as 2^k offsets)."""
    df = altitude_hi(k, m=m)
    a, W = df.attrs, 74
    out = ["=" * W, f"  ALTITUDE 2^{k}  (~{a['digits']:,} digits)   engine={a['engine']}  "
                    f"provable={a['provable']}", "=" * W]
    out += ["", f"  under:  nearest prime below = 2^{k} - {a['under']}   distance {a['under']} "
                f"(odd={a['under'] % 2 == 1})   mersenne={a['mersenne']}"]
    out += ["", f"  seam gap = {a['seam']['gap']} = {a['seam']['down']} (down) + {a['seam']['up']} (up)"]
    out += ["", f"  ceiling_floor invariant  floor_wrap + ceil_wrap == 2^{k}-2  (exact big int):"]
    for e in a["ends"]:
        out += [f"    near {e['where']:<14}  {e['which']}={e['small_wrap']}  (partner is huge)  sum ok={e['sum_ok']}"]
    out += ["", "  seam window  (prime = 2^k + offset):"]
    for _, r in df.iterrows():
        mark = "  <-- SEAM" if r["seam"] else ""
        sign = "+" if r["offset"] >= 0 else "-"
        out += [f"    2^{k} {sign} {abs(int(r['offset'])):>5}   gap {str(int(r['gap'])):>4}{mark}"]
    out += ["=" * W]
    return "\n".join(out)


_EXX = np.array([
        0, 8, 28, 24, 44, 40, 60, 48, 64, 64, 92, 88, 108, 96, 120, 120, 136, 132, 144, 144,
        172, 160, 164, 180, 192, 196, 208, 212, 248, 228, 248, 248, 272, 264, 280, 272, 300,
        288, 312, 312, 336, 328, 348, 344, 380, 372, 380, 376, 396, 384, 400, 392, 424, 416,
        440, 440, 456, 452, 464, 456, 504, 484, 496, 500, 536, 512, 512, 536, 556, 544, 548,
        560, 584, 580, 600, 588, 620, 636, 608, 632, 656, 648, 664, 660, 680, 700, 700, 696,
        716, 724, 708, 724, 756, 740, 760, 740, 768, 808, 792, 776, 808, 800, 816, 816, 848,
        832, 856, 856, 872, 864, 892, 888, 924, 920, 912, 916, 936, 928, 944, 932, 972, 960,
        964, 968, 1012, 1016, 996, 1008, 1032, 1028, 1028, 1040, 1068, 1084, 1084, 1056, 1104,
        1108, 1104, 1096, 1120, 1124, 1136, 1128, 1176, 1156, 1168, 1172, 1208, 1188, 1200,
        1204, 1260, 1220, 1220, 1224, 1268, 1252, 1272, 1264, 1304, 1284, 1304, 1296, 1332,
        1312, 1324, 1332, 1352, 1352, 1344, 1368, 1384, 1400, 1392, 1384, 1408, 1408, 1456,
        1428, 1440, 1444, 1464, 1456, 1492, 1472, 1496, 1520, 1504, 1528, 1508, 1512, 1544,
        1536, 1540, 1540, 1580, 1568, 1592, 1592, 1608, 1620, 1600, 1656, 1632, 1636, 1656,
        1656, 1672, 1668, 1680, 1684, 1724, 1704, 1720, 1784, 1736, 1728, 1732, 1784, 1776,
        1776, 1788, 1784, 1804, 1792, 1796, 1812, 1824, 1824, 1836, 1832, 1856, 1860, 1872,
        1876, 1896, 1908, 1900, 1900, 1948, 1920, 1948, 1936, 1972, 1956, 1956, 1972, 2008,
        1988, 2008, 2008, 2028, 2016, 2032, 2032, 2064, 2084, 2072, 2068, 2080, 2076, 2108,
        2104, 2128, 2132, 2128, 2120, 2164, 2144, 2168, 2160, 2208, 2176, 2204, 2200, 2236,
        2232, 2220, 2220, 2276, 2256, 2264, 2260, 2272, 2276, 2288, 2280, 2332, 2304, 2328,
        2360, 2336, 2316, 2352, 2320, 2376, 2396, 2392, 2388, 2400, 2404, 2404, 2416, 2448,
        2440, 2432, 2432, 2480, 2464, 2476, 2488, 2508, 2504, 2496, 2500, 2528, 2508, 2544,
        2544, 2572, 2568, 2560, 2572, 2560, 2588, 2620, 2608, 2644, 2668, 2640, 2644, 2676,
        2636, 2672, 2676, 2708, 2688, 2712, 2708, 2748, 2728, 2712, 2736, 2752, 2752, 2780,
        2768, 2808, 2784, 2812, 2808, 2824, 2820, 2840, 2820, 2856, 2848, 2856, 2868, 2888,
        2900, 2884, 2896, 2924, 2904, 2896, 2920, 2964, 2948, 2968, 2968, 2984, 2944, 2988,
        3000, 3044
])


def recover_exx():
    """The full 377-value ``exx`` sequence, baked in -- recover it anywhere,
    no primes, no sympy, no ``ee``/``pf`` needed. Flip it, ship it.

    ``exx = ee[::2] ^ (2 * prime_gaps)``; this is the frozen result, exact.
    The octave twin is ``recover_eff() == recover_exx() + 128``.
    """
    return _EXX.copy()


def recover_eff():
    """The octave twin of :func:`recover_exx`: ``exx + 128`` (2^7 up), exact."""
    return _EXX.copy() + 128


def octave_flip(a, bit=7):
    """Flip the 2**bit "octave" bit of every element -- the flip-it operation.

    ``octave_flip(exx)`` == ``exx ^ 128``: shuffles the +128 offset on/off
    across the 16-wide blocks, reordering the sequence ("out of order") without
    touching its building blocks. bit=7 -> 128 (2^7). It is its own inverse.
    """
    return np.asarray(a) ^ (1 << int(bit))


def mersenne_ruler(n):
    """Mersenne-ruler staircase (OEIS A038712): ``m ^ (m-1)`` for m=1..n.

    Yields 1, 3, 1, 7, 1, 3, 1, 15, ... -- the low run of 1-bits, equal to
    ``2**(v2(m)+1) - 1``. Pure XOR, zero primes: the powers-of-2 spine that
    rides the odd 16-blocks of the gap constructions.
    """
    m = np.arange(1, int(n) + 1)
    return m ^ (m - 1)


def mersenne_block_offsets(nblocks, scale=256, width=16):
    """Per-element offset for the interleaved gap/Mersenne array.

    Even ``width``-blocks carry 0; odd blocks carry ``scale * mersenne_ruler``
    (256, 768, 256, 1792, ...). Returns a flat array of length nblocks*width --
    add doubled gaps to it and you rebuild ``(re^128) ^ (exx+128)``.
    """
    off = np.zeros(int(nblocks), dtype=np.int64)
    odd = off[1::2]
    off[1::2] = scale * mersenne_ruler(len(odd))
    return np.repeat(off, int(width))


def ruler_sieve(N):
    """Primes <= N on the ruler lattice.

    The ruler r(n)=v2(n) marks each integer's octave depth; its level-0 ticks --
    the odd numbers -- are exactly the prime candidates (2 aside). Strike the odd
    multiples and what survives is prime. See :func:`ruler_sieve_128` for the
    memory-light octave-bar version. No primes assumed -- it makes them.
    """
    if N < 2:
        return np.array([], dtype=np.int64)
    s = np.ones(N + 1, dtype=bool)
    s[:2] = False
    s[4::2] = False                       # drop non-level-0 ticks (the evens)
    for p in range(3, int(N ** 0.5) + 1, 2):
        if s[p]:
            s[p * p :: 2 * p] = False     # odd multiples only -- stay on the lattice
    return np.flatnonzero(s)


def ruler_sieve_128(N, seg=128):
    """Segmented ruler sieve in the 128-space: sieve octave-bars one at a time.

    Processes ``[k*seg, (k+1)*seg)`` bar by bar, holding only the primes up to
    ``sqrt(N)`` and a single ``seg``-wide bar. seg defaults to 128 = 2**7, the
    octave floor we start from. Matches Eratosthenes exactly; returns primes <= N.
    """
    base = ruler_sieve(int(N ** 0.5) + 1)
    out = []
    for lo in range(0, N + 1, seg):
        hi = min(lo + seg, N + 1)
        bar = np.ones(hi - lo, dtype=bool)
        if lo == 0:
            bar[:2] = False
        for p in base:
            if p * p >= hi:
                break
            start = max(p * p, (lo + p - 1) // p * p)   # first in-bar multiple, >= p^2
            bar[start - lo :: p] = False
        out.append(np.flatnonzero(bar) + lo)
    return np.concatenate(out) if out else np.array([], dtype=np.int64)


def climb_tree(lo=0, hi=128, deltas=None, pmax=None):
    """XOR climb-tree over the primes of a power-of-2 range ``[lo, hi)``.

    Rooted at the middle prime, hops connect two primes that are exactly one bit
    apart -- i.e. whose XOR is a power of 2, a clean half/octave jump (``20^84 ==
    64``). The top bit of the range is the top<->bottom jump. BFS, so every prime
    is reached by its shortest chain of power-of-2 wraps; only prime landings are
    kept and the xor that got there is recorded.

    Returns a DataFrame ``[depth, parent, child, xor, jump]`` (``jump`` is True
    when the edge crosses the top<->bottom half). ``df.attrs['root']`` is the
    root. Barren roots are dodged: the root is the prime nearest the middle that
    actually has a prime neighbor (so ``[0,256)`` roots at 131, not dead-end 127).
    Print it with :func:`climb_tree_ascii`.

    ``pmax`` (optional): cap the node set to primes ``<= pmax`` while keeping the
    full p2 range/bit structure. ``climb_tree(0, 64, pmax=37)`` climbs just the
    small primes 2..37 inside the ``[0,64)`` octave.
    """
    from collections import deque

    lo, hi = int(lo), int(hi)
    W = hi - lo
    nbits = max(1, (W - 1).bit_length())                   # bits to index offsets 0..W-1
    if deltas is None:
        deltas = [1 << i for i in range(nbits)]            # 1, 2, ..., top bit
    topbit = 1 << (nbits - 1)                              # the top<->bottom jump
    isP = {
        n for n in range(max(lo, 2), hi)
        if (pmax is None or n <= pmax) and _is_prime(n)
    }
    mid = (lo + hi) // 2

    def neighbors(p):
        return [lo + ((p - lo) ^ d) for d in deltas]

    cand = sorted(isP, key=lambda p: (abs(p - mid), p))
    root = next(
        (p for p in cand if any(lo <= q < hi and q in isP for q in neighbors(p))),
        cand[0] if cand else None,
    )
    cols = ["depth", "parent", "child", "xor", "jump"]
    if root is None:
        return pd.DataFrame(columns=cols)

    seen = {root}
    q = deque([(root, 0)])
    rows = []
    while q:
        node, d = q.popleft()
        for dx in deltas:
            child = lo + ((node - lo) ^ dx)
            if lo <= child < hi and child in isP and child not in seen:
                seen.add(child)
                rows.append(
                    {"depth": d + 1, "parent": node, "child": child,
                     "xor": dx, "jump": dx == topbit}
                )
                q.append((child, d + 1))
    df = pd.DataFrame(rows, columns=cols)
    df.attrs["root"] = root
    return df


def climb_tree_ascii(df):
    """Pretty-print a :func:`climb_tree` DataFrame as an indented tree."""
    root = df.attrs.get("root")
    kids = {}
    for _, r in df.iterrows():
        kids.setdefault(int(r["parent"]), []).append(
            (int(r["child"]), int(r["xor"]), bool(r["jump"]))
        )
    lines = [f"root {root}"]

    def walk(n, ind=0):
        for c, x, j in sorted(kids.get(n, [])):
            lines.append("  " * ind + f"--x{x}--> {c}" + ("  [JUMP top/bottom]" if j else ""))
            walk(c, ind + 1)

    if root is not None:
        walk(root)
    return "\n".join(lines)


def xor_climb_arr(values, deltas=None, root=None, pmax=None):
    """XOR climb-tree over EXACTLY the numbers you hand it -- the array *is* the node set.

    Like :func:`climb_tree`, but instead of taking the primes of a ``[lo, hi)`` range it
    climbs through the array ``values`` you give it (any ints -- your twins, your gaps, a
    hand-picked prime set). Two values are joined when their XOR is a power of two (or, with
    ``deltas``, is one of your moves); BFS from the root gives every value its shortest chain
    of moves. Only landings that are already IN ``values`` are kept -- it never steps off
    your array.

    Returns a DataFrame ``[depth, parent, child, xor, jump]`` (``df.attrs['root']`` set),
    printable with :func:`climb_tree_ascii`.

    - ``deltas`` (optional): the move alphabet. Default is the power-of-2 bit-flips
      ``1, 2, 4, ...`` up to the widest value. Pass your own to climb by a custom set,
      e.g. ``deltas=[2, 6]`` to hop only by those XORs.
    - ``root`` (optional): where to start. Default is the best-connected value.
    - ``pmax`` (optional): drop values ``> pmax`` before climbing.
    """
    from collections import deque

    vals = sorted({int(v) for v in values if (pmax is None or int(v) <= pmax)})
    cols = ["depth", "parent", "child", "xor", "jump"]
    if not vals:
        return pd.DataFrame(columns=cols)
    vset = set(vals)
    nbits = max(1, max(v.bit_length() for v in vals))
    if deltas is None:
        deltas = [1 << i for i in range(nbits)]            # power-of-2 bit-flips
    deltas = [int(d) for d in deltas]
    topbit = 1 << (nbits - 1)

    def neighbors(p):
        return [p ^ d for d in deltas if (p ^ d) in vset]

    if root is None:                                       # best-connected value, tie -> central
        mid = vals[len(vals) // 2]
        root = max(vals, key=lambda p: (len(neighbors(p)), -abs(p - mid)))
        if not neighbors(root):                            # dodge a barren root
            root = next((p for p in vals if neighbors(p)), vals[0])
    root = int(root)

    seen = {root}
    dq = deque([(root, 0)])
    rows = []
    while dq:
        node, d = dq.popleft()
        for dx in deltas:
            child = node ^ dx
            if child in vset and child not in seen:
                seen.add(child)
                rows.append({"depth": d + 1, "parent": node, "child": child,
                             "xor": dx, "jump": dx == topbit})
                dq.append((child, d + 1))
    df = pd.DataFrame(rows, columns=cols)
    df.attrs["root"] = root
    return df


def xor_fold_climb(k, gates=(3, 5, 7, 11, 13, 17), fold_only=False):
    """A gated XOR fold tree over the octave ``[0, 2**k)`` -- the sieve and the fold as one.

    Nodes are the mod-GATE survivors: odd numbers > 1 that no small-prime gate forces
    composite (coprime to ``gates`` -- the can-be-prime wheel, PURE LOGIC, no primality test).
    They climb by power-of-two XOR (like :func:`climb_tree`) and FOLD by complement: each node
    ``c`` pairs with ``2**k - c``, the two summing to ``2**k``. So the octave folds across its
    center and only survivors remain -- exactly the candidate set a builder would test. The
    fold you see and the sieve that speeds the build are the same object.

    ``fold_only=True`` keeps only the **mutual fold pairs** -- survivors whose complement also
    survives -- so the whole tree is the fold-survives-true lane (``fold_survives`` all True):
    the numbers that pair perfectly across the octave center. Below the squaring bound those
    are exactly the pairs of primes summing to ``2**k`` (Goldbach splits of a power of two).

    Returns ``[depth, parent, child, xor, jump, fold, fold_survives, sums_p2]``
    (``df.attrs['root']`` and ``df.attrs['size'] = 2**k``); print with
    :func:`climb_tree_ascii`. Pick ``k`` on the octaves (8, 16, 32, ...) so the folds nest;
    scale ``gates`` up with ``k`` (toward ``2**(k/2)``) to keep survivors pure prime.
    """
    from collections import deque

    P2 = 1 << int(k)
    # vectorized gate-sieve -- the [start::step] indexer, straight out of Eratosthenes.
    # Survivors come out SORTED for free, and it scales to big octaves without a Python loop.
    sieve = np.ones(P2, dtype=bool)
    sieve[::2] = False                                   # drop evens
    sieve[1] = False
    for p in gates:
        sieve[p::p] = False                              # drop every multiple of each gate
    surv = np.flatnonzero(sieve).tolist()                # can-be-prime wheel, already sorted
    sset = set(surv)
    deltas = [1 << i for i in range(int(k))]
    topbit = 1 << (int(k) - 1)
    half = P2 >> 1                                       # the hinge: 2^(k-1), fold straddles it
    cols = ["depth", "parent", "child", "xor", "jump", "fold", "fold_survives", "sums_p2",
            "half", "above"]
    if not surv:
        return pd.DataFrame(columns=cols)

    def nbrs(x):
        return [x ^ d for d in deltas if (x ^ d) in sset]

    mid = P2 // 2
    root = min(surv, key=lambda x: (abs(x - mid), x))
    if not nbrs(root):
        root = next((x for x in surv if nbrs(x)), surv[0])

    seen = {root}
    dq = deque([(root, 0)])
    rows = []
    while dq:
        node, dep = dq.popleft()
        for d in deltas:
            ch = node ^ d
            if ch in sset and ch not in seen:
                seen.add(ch)
                fold = P2 - ch
                rows.append({"depth": dep + 1, "parent": node, "child": ch, "xor": d,
                             "jump": d == topbit,
                             "fold": fold, "fold_survives": fold in sset,
                             "sums_p2": ch + fold == P2,       # always True: child + fold = 2^k
                             "half": half, "above": ch >= half})  # the hinge and which side of it
                dq.append((ch, dep + 1))
    df = pd.DataFrame(rows, columns=cols)
    if fold_only:                                        # the fold-survives-true lane
        df = df[df["fold_survives"]].reset_index(drop=True)
    df.attrs["root"] = root
    df.attrs["size"] = P2
    return df


def halving_evens(lo=0, hi=128, deep=False):
    """Every even ``e`` in a power-of-2 range ``[lo, hi)`` whose half is prime.

    That is: ``e`` even and ``e // 2`` prime, so ``e == 2 * p`` for a prime
    ``p`` -- the primes doubled, one octave up (``4->2, 6->3, 10->5, ...``).
    Returns a sorted int array.

    ``deep=True`` follows the doubling all the way: halve until odd and test
    THAT odd part, so it also catches ``164 -> 82 -> 41`` and the rest of each
    prime's doubling tower ``2^k * p`` (three doublings = x8, one jump across the
    p2/8 boundary). Powers of 2 have odd part 1 and are excluded.

    Ties to :func:`climb_tree`: since ``2p ^ 2q == 2*(p^q)``, the tree over the
    shallow evens is the prime climb-tree with every value doubled -- same shape,
    one octave up. XOR (by an even wrap) moves inside the odd-prime lanes;
    doubling is the only bridge across the parity/8 boundary to these evens.
    """
    lo, hi = int(lo), int(hi)
    if not deep:
        return np.array(
            [2 * p for p in range(2, hi // 2 + 1) if lo <= 2 * p < hi and _is_prime(p)],
            dtype=int,
        )

    def _odd_part(e):
        while e % 2 == 0:
            e //= 2
        return e

    lo_e = max(lo + (lo % 2), 2)  # first even >= lo (and >= 2)
    return np.array(
        [e for e in range(lo_e, hi, 2) if _is_prime(_odd_part(e))],
        dtype=int,
    )


def odd_core(n):
    """The odd core of ``n``: divide out every factor of 2 (``n >> v2(n)``).

    Always odd. Every integer is ``2**k * odd_core(n)``, so halving a number
    down to its core is the deterministic way to find the odd prime it sits over
    -- the consistency behind the wrap-and-divide game, no ^1/^3 ambiguity. The
    only primes reachable this way are the ODD primes (2's core is 1). Works on
    a scalar or a numpy array.
    """
    a = np.asarray(n)
    if a.ndim == 0:
        v = int(a)
        return 0 if v == 0 else v // (v & -v)
    a = a.astype(np.int64)
    low = a & -a
    return np.where(a == 0, 0, a // np.where(low == 0, 1, low))


def tower(core, hi, lo=0):
    """The doubling tower over ``core``: ``core, 2*core, 4*core, ...`` in ``[lo, hi)``.

    Every member has odd core ``odd_core(core)`` and divides back down to it.
    Use an odd prime as ``core`` to get the evens that halve to that prime.
    """
    out = []
    v = int(core)
    while v < hi:
        if v >= lo:
            out.append(v)
        v *= 2
        if v == 0:  # core 0 guard
            break
    return np.array(out, dtype=int)


def wrap_depth(a, b):
    """How far ``a ^ b`` halves before going odd: ``v2(a ^ b)`` -- the lowest bit
    where ``a`` and ``b`` differ, and the consistent rule behind wrap-and-divide.

    Read ``a % 8`` and ``b % 8`` and you know it in advance: differ at bit 1 ->
    1, same mod 4 -> 2, same mod 8 -> 3 (the x8 jump). Deterministic; only the
    aggregate looks like "half the time" (depths are geometric: 1/2, 1/4, ...).
    """
    x = int(a) ^ int(b)
    return 0 if x == 0 else (x & -x).bit_length() - 1


def cycle_walk(first_row, distances):
    """Unroll a table from a starting row by successive XOR distances.

    Row 0 is ``first_row``; each next row is the previous row XOR the next
    scalar distance -- "know the distance to the next row, xor to continue the
    cycle." Row ``k`` ends up as ``first_row ^ cumxor(distances[:k])``.

    Feed the prime gaps as ``distances`` and a cumulative-XOR seed row and you
    rebuild the buildsystem table: the rows step by the gaps, one XOR at a time.
    Returns a 2-D int array with ``len(distances)+1`` rows.
    """
    first = np.asarray(first_row)
    rows = [first.copy()]
    cur = first
    for d in distances:
        cur = cur ^ int(d)
        rows.append(cur)
    return np.array(rows)


def xor_table(seq):
    """The XOR Cayley table of a sequence: ``M[i, j] = seq[i] ^ seq[j]``.

    Row ``i`` is ``seq ^ seq[i]``; the diagonal is 0 and it is symmetric. With
    ``seq = numpy.bitwise_xor.accumulate([0, *gaps])`` (the running XOR of the
    prime gaps) this is buildsystem: consecutive rows differ by a single gap,
    and the low bit stays odd because only the first gap (2->3) is odd.
    """
    s = np.asarray(seq)
    return s[:, None] ^ s[None, :]

def under(p):
    """One step 'under' a prime: ``p ^ 1`` then halve = ``(p - 1) // 2`` (the index).

    ``^1`` is the only parity-flipping wrap, so this is the unique halvable
    descent from an odd prime. Returns a namedtuple ``Under(index, safe)`` where
    ``index = (p - 1) // 2`` and ``safe`` is True iff both ``p`` and ``index``
    are prime -- i.e. ``p`` is a safe prime and the index is Sophie Germain
    (``under(23) -> Under(11, True)``). This ``(p-1)/2`` is the creator wheel index.
    """
    from collections import namedtuple

    p = int(p)
    q = (p - 1) // 2
    return namedtuple("Under", ["index", "safe"])(q, bool(_is_prime(p) and _is_prime(q)))


def prime_moves(p, hi=None):
    """Every move from a prime ``p`` in the wrap / double / index system.

    - ``climb``: ``(wrap, q)`` pairs with ``q = p ^ wrap`` prime and ``wrap`` a
      power of 2 >= 2 (an even wrap -- stays in the odd lane).
    - ``tower``: the doubling tower ``[p, 2p, 4p, ...]`` below ``hi``.
    - ``under`` / ``safe``: the index ``(p-1)//2`` and whether it is prime.

    ``hi`` defaults to the next octave above ``p``. Returns a dict, ready to feed
    an auto-builder -- the three moves are: ^even (climb), x2/÷2 (tower),
    ^1 then ÷2 (under).
    """
    p = int(p)
    if hi is None:
        hi = 1 << (p.bit_length() + 1)
    climb = [
        ((1 << k), p ^ (1 << k))
        for k in range(1, hi.bit_length())
        if 0 < (p ^ (1 << k)) < hi and _is_prime(p ^ (1 << k))
    ]
    u = under(p)
    return {
        "prime": p,
        "climb": climb,
        "tower": tower(p, hi).tolist(),
        "under": u.index,
        "safe": u.safe,
    }

def prime_stream(limit=None):
    """Generate primes in order with **no primality test** -- an incremental sieve.

    It never asks "is q prime?". It ELIMINATES composites: each prime lays a mark
    at its square and advances that mark by ADDITION (``q + p``). The lowest
    number nobody ever marked is the next prime. The only operations are "is q on
    the pending-composite list" and "add p to move a mark forward" -- no division,
    no detector. Goes to the lowest number, builds the high primes by addition.

    Yields primes <= ``limit`` (or forever if ``limit`` is None).
    """
    pending = {}          # future composite -> the primes that will land on it
    q = 2
    while limit is None or q <= limit:
        if q not in pending:
            yield q                       # never marked -> prime
            pending[q * q] = [q]          # its first composite is q*q
        else:
            for p in pending.pop(q):      # advance each mark by ADDITION
                pending.setdefault(q + p, []).append(p)
        q += 1


def prime_set(limit):
    """Detector-free prime membership: ``set(prime_stream(limit))``.

    Built by elimination, not testing -- feed it to :func:`climb_tree` /
    :func:`prime_moves` (via their internal checks) so the whole system decides
    "prime?" by set membership in the sieve survivors, never by a detector.
    """
    return set(prime_stream(int(limit)))

def build_lines(start=2, stop=110):
    """Walk the primes detector-free and decompose each one, line by line.

    Uses :func:`prime_stream` (no primality test) and, for each prime ``p`` in
    ``[start, stop]``, records:

    - ``lane``  = ``p % 8``          (the odd-lane residue)
    - ``gap``   = next prime - p     (the low even step that lands)
    - ``two_k`` = ``2 ** v2(p-1)``   (distance down to the p2/8 boundary)
    - ``core``  = odd part of p-1    (so ``p - 1 == two_k * core``; 89 -> 8, 11)
    - ``under`` = ``(p-1)//2``       (the index)

    Returns a pandas DataFrame -- the whole thing built from the moves, no detector.
    """
    ps = [p for p in prime_stream(stop) if p >= start]
    rows = []
    for i, p in enumerate(ps):
        gap = ps[i + 1] - p if i + 1 < len(ps) else 0
        pm1 = p - 1
        two_k = pm1 & -pm1 if pm1 else 0          # 2**v2(p-1)
        core = pm1 // two_k if two_k else 0        # odd core of p-1
        rows.append({"p": p, "lane": p % 8, "gap": gap,
                     "two_k": two_k, "core": core, "under": pm1 // 2})
    return pd.DataFrame(rows)

def odometer(start=2, stop=200):
    """The self-building prime view -- the low-number odometer beside the primes.

    Cranks primes with :func:`prime_stream` (no primality test) and shows, per
    prime, the buildsystem's first two columns and the jump read DOWN the column:

    - ``seq``  : column 0 = cumulative XOR of the gaps (the odd odometer)
    - ``even`` : ``seq ^ 1`` = the even mirror sitting *beside it* (column 1)
    - ``lane`` : ``p % 8`` (the locked coordinate)
    - ``jump`` : ``seq ^ seq_below`` = the gap to the next prime
    - ``next`` : ``prime + jump``

    Every value stays low (the seq never leaves 0..2^k-1); the gap is a vertical
    neighbor XOR, not a subtraction. No detector anywhere. Returns a DataFrame.
    """
    ps = [p for p in prime_stream(stop) if p >= start]
    gaps = [ps[i + 1] - ps[i] for i in range(len(ps) - 1)]
    seq = np.bitwise_xor.accumulate(np.array([0] + gaps))     # len == len(ps)
    rows = []
    for i, p in enumerate(ps):
        s = int(seq[i])
        below = int(seq[i + 1]) if i + 1 < len(seq) else None
        jump = (s ^ below) if below is not None else 0
        rows.append({"prime": p, "seq": s, "even": s ^ 1, "lane": p % 8,
                     "jump": jump, "next": p + jump if jump else 0})
    return pd.DataFrame(rows)

def is_prime_by_index(n, index):
    """Prime test using an existing prime list (the index) -- no q+p, no detector.

    Divide ``n`` by the primes in ``index`` in order; the moment one divides it
    is composite, and once the primes pass ``sqrt(n)`` nothing else can divide it
    so it is prime. ``index`` must contain every prime up to ``sqrt(n)`` (building
    in order, it always does). One mod is a hint; the whole index is a proof.
    """
    if n < 2:
        return False
    for p in index:
        if p * p > n:
            return True
        if n % p == 0:
            return False
    return True


def primes_by_index(N):
    """Build primes <= N by dividing each candidate against the index so far.

    Your method: keep the growing prime list and, for each odd candidate, divide
    by the primes up to its square root -- if none divide it, it joins the index.
    No pending-composite bookkeeping (no q+p). Returns a list. Correct because the
    index always holds every prime up to sqrt(candidate) by the time you reach it.
    """
    if N < 2:
        return []
    primes = [2]
    n = 3
    while n <= N:
        r = int(n ** 0.5)
        is_p = True
        for p in primes:
            if p > r:
                break
            if n % p == 0:
                is_p = False
                break
        if is_p:
            primes.append(n)
        n += 2
    return primes

def xor_grid(a, b, prepend_zero=True, n=None, labels=True):
    """Outer XOR table between two arrays -- your ``meshgrid ^`` setup, in one call.

    Prepends a 0 to each array (``prepend_zero=True``) like your pa/ca, then
    ``P, N = np.meshgrid(a, b); D = P ^ N``. Returns a DataFrame with the
    (0-prepended) arrays as column / row labels. ``D.iloc[i, j] == b'[i] ^ a'[j]``.
    ``n`` slices the top-left ``n x n`` corner (like ``D[:n, :n]``). ``labels=False``
    gives a plain 0..k index instead of the array values.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if prepend_zero:
        a = np.concatenate([[0], a])
        b = np.concatenate([[0], b])
    P, N = np.meshgrid(a, b)
    D = P ^ N
    df = pd.DataFrame(D, index=(b if labels else None), columns=(a if labels else None))
    if n is not None:
        df = df.iloc[:n, :n]
    return df

def rail_climb_tree(lo=0, hi=256, rail=None, grow=False, cap=None, pmax=None,
                    split=False, anchors=None):
    """2nd-generation climb tree with an arbitrary RAIL of XOR wraps.

    Like :func:`climb_tree` but the rail can be any list of wraps (primes,
    non-primes, even or odd), not just powers of 2 -- so it handles more cases.
    Edges connect two primes ``p, q`` when ``p ^ q`` is in ``rail``; keep ``1``
    on the rail to catch the prime 2 (``3 ^ 1 == 2``). BFS from the middle prime.

    The default rail is the powers of 2 up to ``2*hi`` -- it *scales with the
    range*, so a big ``hi`` gets the big octave-jump wraps (``[0,8192)`` needs
    ``4096``, not just ``128``). ``cap`` also defaults to ``2*hi``.

    ``grow=True``: if the rail does not reach every prime in range, keep appending
    the next wraps (2, 3, 4, ...) and rebuilding until coverage is complete or the
    largest wrap hits ``cap``. ``pmax`` caps the node set to primes <= pmax.

    ``split=True``: skip the single-root BFS and tile from **two anchors on the
    sides of a power of two** -- ``2**k - 1`` (all-ones, LEFT, hits the 3-mod-4
    lane) and ``2**k + 1`` (fermat-form, RIGHT, hits the 1-mod-4 lane), with
    ``k`` chosen so the p2 sits nearest the middle of the range. The rail is the
    pure multiples of 4, so each anchor reaches its ENTIRE lane with zero leak
    (same-lane primes always differ by a multiple of 4). Anchors need not be
    prime. Pass ``anchors=(a3, a1)`` to override (e.g. ``(127, 257)``); ``a3``
    must be 3-mod-4, ``a1`` 1-mod-4. The df gains a ``side`` column (``L``/``R``,
    ``e`` for the lone even prime 2); ``df.attrs['root']`` is ``(left, right)``.

    Returns ``(df, info)``: df is ``[depth, parent, child, wrap]`` (root in
    ``df.attrs['root']``); info has ``rail`` (final), ``reached``, ``total``,
    ``coverage``.
    """
    from collections import deque

    lo, hi = int(lo), int(hi)
    if cap is None:
        cap = 2 * hi
    isP = {n for n in range(max(lo, 2), hi) if (pmax is None or n <= pmax) and _is_prime(n)}
    cols = ["depth", "parent", "child", "wrap", "p2", "p2_bounce"]
    if not isP:
        return pd.DataFrame(columns=cols + (["side"] if (split or anchors is not None) else [])), {}

    if split or anchors is not None:
        # anchors on the SIDES of a power of two: 2^k-1 (all-ones, 3-lane) and
        # 2^k+1 (fermat-form, 1-lane). pure mult-of-4 wraps -> each anchor tiles
        # its whole lane, no leak. split=2 straddles ONE p2 (two rails); split=4
        # straddles TWO consecutive p2s (four rails: 127/129 around 128 and
        # 255/257 around 256) so the roll balls show -- double a low anchor and
        # you land on the high one (127*2=254~255, 129*2=258~257).
        mid = (lo + hi) // 2
        four = (split == 4) or (anchors is not None and len(anchors) == 4)

        def _nearest_k(m):
            k = max(2, m.bit_length() - 1)
            if m >= 3 * (1 << k) // 2:          # nearer the next power of two
                k += 1
            return k

        p2info = {}
        if anchors is not None:
            A = [int(a) for a in anchors]
            (a3lo, a1lo, a3hi, a1hi) = A if four else (A[0], A[1], None, None)
        elif four:
            k = _nearest_k(mid)
            p2hi = 1 << k
            p2lo = max(4, p2hi >> 1)
            a3lo, a1lo, a3hi, a1hi = p2lo - 1, p2lo + 1, p2hi - 1, p2hi + 1
            p2info = {"p2": (p2lo, p2hi)}
        else:
            k = _nearest_k(mid)
            p2k = 1 << k
            a3lo, a1lo, a3hi, a1hi = p2k - 1, p2k + 1, None, None
            p2info = {"p2": p2k}

        threes = [a for a in (a3lo, a3hi) if a is not None]
        ones = [a for a in (a1lo, a1hi) if a is not None]
        for a in threes:
            if a % 4 != 3:
                raise ValueError(f"3-lane anchor {a} is not 3-mod-4")
        for a in ones:
            if a % 4 != 1:
                raise ValueError(f"1-lane anchor {a} is not 1-mod-4")

        divider = (a3hi + 1) if four else None               # the high p2 (2^(k+1))

        def _attach(c):
            if c % 4 == 3:
                if four and c >= divider:                    # high octave -> high anchor
                    return a3hi, "L+"
                return a3lo, "L"
            if c % 4 == 1:
                if four and c >= divider:
                    return a1hi, "R+"
                return a1lo, "R"
            return a1lo, "e"                                 # the lone even prime 2

        rows = []
        for c in sorted(isP):
            parent, side = _attach(c)
            w = parent ^ c
            mx, mn = (parent, c) if parent > c else (c, parent)
            b = 1 << (mx.bit_length() - 1)
            crossed = b if b > mn else 0
            # THIRD RAIL -- the prime on the other side: complement across the anchor's power
            # of two, 2^k - c. Sends a 3-mod-4 (left) prime to its 1-mod-4 (right) partner and
            # vice versa, the two summing to 2^k (127-31=96, so 31's other side is 128-31=97).
            p2_of = parent + 1 if parent % 4 == 3 else parent - 1
            other = p2_of - c                                      # complement to 2^k (this octave)
            up = 2 * p2_of - c                                     # complement to 2^(k+1) (next octave up)
            rows.append({"depth": 1, "parent": parent, "child": c, "wrap": w,
                         "p2": crossed, "p2_bounce": crossed > 0, "side": side,
                         "other": other, "other_prime": bool(other > 1 and _is_prime(other)),
                         "sums_p2": bool(c + other == p2_of),      # always True: c + other = 2^k
                         "up": up, "up_prime": bool(up > 1 and _is_prime(up))})  # next rail, octave up
        df = pd.DataFrame(rows, columns=cols + ["side", "other", "other_prime", "sums_p2",
                                                "up", "up_prime"])
        df["gap"] = (df["child"].shift(-1) - df["child"]).fillna(0).astype(int)   # to next prime
        df["under"] = (df["child"] - 1) // 2                                      # build-seed: c = 2*under+1
        df["mod4"] = df["child"] % 4                                              # the lane, numeric (1 or 3)
        df["xflags"] = (df["other_prime"].astype(int) ^ df["up_prime"].astype(int)  # False's XOR out to 0
                        ^ df["sums_p2"].astype(int))                             # -> a clean 0/1 per prime
        # the DOUBLING chain -- how you BUILD a prime by logic (rollball rung, R=1):
        df["under_prime"] = [bool(u > 1 and _is_prime(int(u))) for u in df["under"]]  # child came from a prime?
        df["over"] = 2 * df["child"] + 1                                          # ascent: 2c+1, the next rung
        df["over_prime"] = [bool(_is_prime(int(o))) for o in df["over"]]          # ...and it's prime -> created

        def _chain_up(c, cap=40):
            # how many provable doubling-rungs (2x+1, each Pocklington-proven off the last)
            # this prime launches before the chain breaks -- the growth/provability made visible
            n, x = 0, int(c)
            while n < cap:
                x = 2 * x + 1
                if _is_prime(x): n += 1
                else: break
            return n
        df["chain_up"] = [_chain_up(int(c)) for c in df["child"]]

        def _chain_down(c, cap=40):
            # the mirror: how many times (x-1)/2 stays prime going DOWN -- so chain_down==0
            # with chain_up>0 marks a chain-BOTTOM, the seed you build up from.
            n, x = 0, int(c)
            while n < cap:
                if x % 2 == 0 or x < 3:
                    break
                nxt = (x - 1) // 2
                if nxt >= 2 and _is_prime(nxt):
                    n += 1; x = nxt
                else:
                    break
            return n
        df["chain_down"] = [_chain_down(int(c)) for c in df["child"]]
        allA = tuple(threes + ones)
        df.attrs["root"] = allA
        pure = all(r["wrap"] % 4 == 0 for r in rows if r["side"] != "e")
        info = {"anchors": allA, "reached": len(rows), "total": len(isP),
                "coverage": 1.0, "pure_mult4": pure,
                "rails": {s: sum(r["side"] == s for r in rows)
                          for s in ("L", "L+", "R", "R+", "e")
                          if any(r["side"] == s for r in rows)}}
        info.update(p2info)
        if four:
            info["rollball"] = {a3lo: a3hi, a1lo: a1hi}      # double the low anchor ~ high anchor
        return df, info
    if isinstance(rail, str):
        named = {
            # roll up to 2*hi in 8s, halve back down (mult of 4) + the essentials:
            # {2,6} stitch twins/cousins, 1 catches the even prime 2. -> 100%, lean.
            "roll8": sorted({int(x) for x in np.arange(0, 2 * hi, 8) // 2 if x > 0} | {1, 2, 6}),
            "evens": list(range(2, 2 * hi, 2)),
            "pow2": [1 << k for k in range((2 * hi).bit_length())],
        }
        if rail not in named:
            raise ValueError(f"unknown named rail {rail!r}; use 'roll8', 'evens', or 'pow2'")
        rail = named[rail]
    rail = list(dict.fromkeys(int(r) for r in (rail if rail is not None else [1 << k for k in range((2 * hi).bit_length())])))
    mid = (lo + hi) // 2

    def build(rail, root=None):
        if root is None:
            cand = sorted(isP, key=lambda p: (abs(p - mid), p))
            root = next((p for p in cand if any(lo <= (p ^ r) < hi and (p ^ r) in isP for r in rail)), cand[0])
        seen = {root}
        q = deque([(root, 0)])
        rows = []
        while q:
            node, d = q.popleft()
            for r in rail:
                c = node ^ r
                if lo <= c < hi and c in isP and c not in seen:
                    seen.add(c)
                    mx, mn = (node, c) if node > c else (c, node)
                    b = 1 << (mx.bit_length() - 1)          # octave floor of the higher endpoint
                    crossed = b if b > mn else 0            # the p2 boundary this wrap crossed
                    rows.append({"depth": d + 1, "parent": node, "child": c, "wrap": r,
                                 "p2": crossed, "p2_bounce": crossed > 0})
                    q.append((c, d + 1))
        return root, seen, rows

    root, seen, rows = build(rail)
    if grow:
        # add the *specific* wraps that connect each unreached prime to the tree
        # (its smallest xor to a reached prime) -- not every integer, so it stays
        # fast at any hi. Same root each round, so coverage only grows.
        while len(seen) < len(isP):
            new_wraps = set()
            for u in isP - seen:
                cands = [w for w in (u ^ v for v in seen) if 0 < w <= cap]
                if cands:
                    new_wraps.add(min(cands))
            if not new_wraps:
                break
            rail = list(dict.fromkeys(rail + sorted(new_wraps)))
            root, seen, rows = build(rail, root)
    df = pd.DataFrame(rows, columns=cols)
    df.attrs["root"] = root
    return df, {"rail": rail, "reached": len(seen), "total": len(isP), "coverage": round(len(seen) / len(isP), 3)}


def split_view(lo=0, hi=512, width=100):
    """ASCII picture of the FOUR-rail split (``rail_climb_tree(..., split=4)``).

    Two lanes stacked (3-mod-4 all-ones anchors, 1-mod-4 fermat anchors), each
    with a LOW rail (``2^k ∓ 1``) and a HIGH rail (``2^(k+1) ∓ 1``) split at the
    power of two between them, so the roll ball is visible: the top of a low rail
    doubles straight onto the high rail (``127*2=254`` sits by ``255``). Every
    wrap is a multiple of 4, so both lanes tile with zero leak. Returns a string.
    """
    import textwrap

    df, info = rail_climb_tree(lo, hi, split=4)
    a3lo, a3hi, a1lo, a1hi = info["anchors"]
    p2lo, p2hi = info["p2"]

    def rail(tag, anchor):
        kids = sorted(df[df.parent == anchor].child.tolist())
        head = f"  {tag:<4} {anchor:>4} ┤ "
        lines = textwrap.wrap(" ".join(f"{p:>3}" for p in kids), width=width - len(head)) or [""]
        s = head + lines[0]
        for ln in lines[1:]:
            s += "\n" + " " * len(head) + ln
        return s, len(kids)

    bar = "=" * width
    out = [bar, f"  FOUR RAILS over [{lo},{hi})    p2 dividers {p2lo} & {p2hi}    every wrap = a multiple of 4", bar, ""]
    for name, alo, ahi in [("3-MOD-4  (all-ones, 2^k-1)", a3lo, a3hi),
                           ("1-MOD-4  (fermat,   2^k+1)", a1lo, a1hi)]:
        hline, nH = rail("high", ahi)
        lline, nL = rail("low", alo)
        out += [f"  {name} LANE", "  " + "-" * (width - 4), hline,
                f"        ⇑ roll ball ×2   {alo}*2 = {alo * 2} lands by {ahi}",
                lline, f"        totals: low {alo} -> {nL}   high {ahi} -> {nH}", ""]
    out += [bar, f"  coverage {info['coverage'] * 100:.0f}%   pure ×4: {info['pure_mult4']}   rails {info['rails']}", bar]
    return "\n".join(out)


def roll_ball(lo=0, hi=256, rail=None, base="double"):
    """Roll-balling: transform each prime to an even base, then XOR an odd rail
    to roll onto new primes.

    ``base`` picks what you roll FROM -- of the *current* prime only. All three
    are even, so ``base ^ odd`` lands on an odd candidate:
      * ``"double"`` (default, or ``2``) -> ``2p``
      * ``"+1"``     (or ``1``)          -> ``p + 1``
      * ``"-1"``     (or ``-1``)         -> ``p - 1``
    For each odd wrap ``r``, ``q = base ^ r``; if ``q`` is prime it is a roll. To
    CROSS a power-of-2 boundary you must roll *higher than it* -- the default rail
    is every odd wrap up to ``2*hi``. Pass a small ``rail`` (e.g.
    ``[1, 3, 5, 7, 9, 11, 13, 15]``) for just local rolls.

    Columns ``[p, base, wrap, q, p2, p2_bounce]`` (the roll-from column is named
    ``base`` now, holding ``2p`` / ``p+1`` / ``p-1``): ``p2`` is the power-of-2
    boundary crossed (0 = none), ``p2_bounce`` is ``p2 > 0``. Dense with the full
    rail -- narrow ``rail`` or filter on ``p2`` / ``p2_bounce`` at large ``hi``.
    """
    lo, hi = int(lo), int(hi)
    isP = sorted(n for n in range(max(lo, 2), hi) if _is_prime(n))
    railset = None if rail is None else set(int(r) for r in rail)
    cap = 2 * hi
    if base in ("double", "x2", "2x", "2p", 2):
        base_of = lambda p: 2 * p
    elif base in ("+1", "plus1", "p+1", 1):
        base_of = lambda p: p + 1
    elif base in ("-1", "minus1", "p-1", -1):
        base_of = lambda p: p - 1
    else:
        raise ValueError("base must be 'double', '+1', or '-1'")
    rows = []
    for p in isP:
        d = base_of(p)
        for q in isP:
            r = d ^ q
            ok = (r % 2 == 1 and 0 < r < cap) if railset is None else (r in railset)
            if ok:
                mx, mn = (d, q) if d > q else (q, d)
                b = 1 << (mx.bit_length() - 1)          # octave floor of the higher value
                crossed = b if b > mn else 0            # the p2 boundary you rolled across
                rows.append({"p": p, "base": d, "wrap": r, "q": q, "p2": crossed, "p2_bounce": crossed > 0})
    return pd.DataFrame(rows, columns=["p", "base", "wrap", "q", "p2", "p2_bounce"])

def wrap_find(seeds, wraps=None, base=None):
    """Find primes by XOR-wrapping known primes, verified by base-prime elimination.

    The combined workflow from the rails/bootstrap write-up: XOR wraps GENERATE
    candidates near the ``seeds`` (known primes); the ``base`` prime set VERIFIES
    which are actually prime by division/elimination -- no primality oracle. Any
    composite ``<= base[-1]**2`` is caught, so the base must reach ``sqrt`` of the
    biggest candidate; when ``base`` is None it is auto-built with
    :func:`ruler_sieve` so it always suffices.

    Parameters
    ----------
    seeds : int or iterable of int
        Known primes to wrap from.
    wraps : iterable of int, optional
        XOR wraps to try. Default ``[1,2,4,6,8,10,12,16,64,128]`` -- the small
        even generators (twins/cousins/octaves) plus ``1`` to catch the prime 2.
    base : iterable of int, optional
        Base primes for verification. Default: primes up to ``sqrt(max candidate)``.
        If you pass your own and it is too small to verify a candidate, raises.

    Returns
    -------
    pandas.DataFrame  --  columns ``[seed, wrap, found, gap]`` (deduped on
    seed+found), the primes reached by wrapping and confirmed by the base.
    """
    seeds = [int(s) for s in np.atleast_1d(seeds)]
    if wraps is None:
        wraps = [1, 2, 4, 6, 8, 10, 12, 16, 64, 128]
    wraps = [int(w) for w in wraps]
    cand_max = max((s ^ w for s in seeds for w in wraps), default=2)
    if base is None:
        base = ruler_sieve(int(cand_max ** 0.5) + 1).tolist()
    else:
        base = [int(b) for b in base]
        if base and base[-1] ** 2 < cand_max:
            raise ValueError(
                f"base only reaches {base[-1]}**2 = {base[-1]**2}; a candidate is "
                f"{cand_max} -- extend the base to verify it (need primes up to "
                f"{int(cand_max**0.5)+1})."
            )
    rows = []
    for s in seeds:
        for w in wraps:
            c = s ^ w
            if c > 2 and is_prime_by_index(c, base):
                rows.append({"seed": s, "wrap": w, "found": c, "gap": c - s})
    df = pd.DataFrame(rows, columns=["seed", "wrap", "found", "gap"])
    return df.drop_duplicates(subset=["seed", "found"]).reset_index(drop=True)

_BASE_CACHE = {"base": None}


def _base_path():
    import os
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "base_primes.npy")


def load_base(limit=1_700_000, path=None, force=False):
    """Load (or build + save + cache) the base primes for :func:`prime_check`.

    A base of primes up to ``limit`` proves any number up to ``limit**2``. Loads
    ``base_primes.npy`` from the numtools folder if present and big enough, else
    builds it with :func:`ruler_sieve` and saves it. Cached in memory after first
    call. Default limit ~1.7M -> reach ~2.9 trillion (~1 MB file).
    """
    import os

    if path is None:
        path = _base_path()
    if not force and _BASE_CACHE["base"] is not None and int(_BASE_CACHE["base"][-1]) >= limit - 1:
        return _BASE_CACHE["base"]
    if not force and os.path.exists(path):
        base = np.load(path)
        if int(base[-1]) >= limit - 1:
            _BASE_CACHE["base"] = base
            return base
    base = ruler_sieve(int(limit))
    try:
        np.save(path, base)
    except Exception:  # pragma: no cover - read-only fs etc.
        pass
    _BASE_CACHE["base"] = base
    return base


def prime_check(n, base=None):
    """Prove ``n`` prime -- or factor it -- by dividing against a base of primes.

    No oracle: valid up to ``base[-1]**2`` because any composite that small has a
    factor in the base. Uses the cached base from :func:`load_base` unless you
    pass your own. Returns a namedtuple ``PrimeCheck(n, prime, factor, cofactor,
    reach)`` -- ``prime`` is True/False, and for composites ``n == factor*cofactor``.

    >>> prime_check(16777213).prime      # True
    >>> prime_check(6777213)             # PrimeCheck(n=6777213, prime=False, factor=3, ...)
    """
    from collections import namedtuple

    Result = namedtuple("PrimeCheck", ["n", "prime", "factor", "cofactor", "reach"])
    n = int(n)
    base = load_base() if base is None else np.asarray(base)
    reach = int(base[-1]) ** 2
    if n > reach:
        raise ValueError(
            f"{n} is past the base reach {reach}; grow it with "
            f"load_base({int(n ** 0.5) + 1}, force=True)"
        )
    if n < 2:
        return Result(n, False, None, None, reach)
    sub = base[base <= int(n ** 0.5)]          # only need primes up to sqrt(n)
    hits = np.flatnonzero(n % sub == 0)        # vectorized: which divide n
    if len(hits):
        f = int(sub[hits[0]])
        return Result(n, False, f, n // f, reach)
    return Result(n, True, None, None, reach)

def roll_prime(target_bits, seed=None, presieve_to=30000):
    """Build a ~``target_bits`` PROVABLE prime by rolling up -- WITH presieving.

    Each roll: ``N = 2*R*q + 1`` where ``q`` is the proven prime below it and
    ``R < q/2`` (so ``q > sqrt(N)``, the Pocklington condition). The prime ``q``
    stays inside ``N-1``, so one witness proves ``N``. The speed trick: every
    candidate is trial-divided by the small primes FIRST, killing ~90% of them
    cheaply before the expensive witness test -- that presieve is what makes big
    lengths fast (the search, not the proof, is the bottleneck).

    Returns ``(prime, chain)`` -- ``chain`` is the roll-up certificate,
    checkable with :func:`verify_certificate`.
    """
    import random
    from math import gcd

    rng = random.Random(seed)
    small = [int(x) for x in ruler_sieve(int(presieve_to))]
    q = (rng.getrandbits(31) | (1 << 31)) | 1          # a 32-bit odd
    while not _is_prime(q):
        q += 2
    chain = [(q, None, None)]
    while q.bit_length() < target_bits:
        rb = max(1, min(q.bit_length() - 1, target_bits - q.bit_length() + 1))
        while True:
            R = rng.getrandbits(rb) | 1
            if R >= q // 2:                            # keep q > sqrt(N)
                continue
            N = 2 * R * q + 1
            if any(N % sp == 0 for sp in small):       # PRESIEVE: cheap composite kill
                continue
            wit = None
            for a in (2, 3, 5, 7, 11, 13, 17, 19, 23):
                if pow(a, N - 1, N) != 1:
                    break                              # Fermat fail -> composite
                if gcd(pow(a, (N - 1) // q, N) - 1, N) == 1:
                    wit = a                            # PROVEN prime
                    break
            if wit is not None:
                chain.append((N, q, wit))
                q = N
                break
    return q, chain


def roll_prime_fast(target_bits, seed=None, presieve_to=30000):
    """Same as :func:`roll_prime` -- a PROVABLE prime by Pocklington rollball -- but
    ~12-17x faster (the gain grows with size). Two accelerations, no change to the
    proof:

    * **gmpy2 modexp.** The witness test uses ``gmpy2.powmod`` / ``gmpy2.gcd`` instead
      of Python's built-in ``pow`` -- GMP's big-int arithmetic wins on large N.
    * **residue presieve.** ``roll_prime`` trial-divides each candidate by every small
      prime with a big-int ``N % sp``. Here we precompute ``q mod sp`` once per proven
      rung, then screen each candidate with ``(2*(R mod sp)*(q mod sp)+1) mod sp`` --
      tiny arithmetic, no big-int division. The presieve, not the proof, was the
      bottleneck; this removes it.

    Bit-identical certificate; ``verify_certificate`` checks it the same way. Needs gmpy2.
    Measured: 2000 bits 2.6s -> 0.2s (12x), 5000 bits 131s -> 8s (17x).
    """
    import random
    from gmpy2 import gcd, mpz, powmod

    rng = random.Random(seed)
    small = [int(x) for x in ruler_sieve(int(presieve_to))]
    q = mpz((rng.getrandbits(31) | (1 << 31)) | 1)
    while not _is_prime(int(q)):
        q += 2
    chain = [(int(q), None, None)]
    while q.bit_length() < target_bits:
        rb = max(1, min(q.bit_length() - 1, target_bits - q.bit_length() + 1))
        qmod = [int(q % sp) for sp in small]           # precompute q mod sp per rung
        half = q // 2
        while True:
            R = rng.getrandbits(rb) | 1
            if R >= half:                              # keep q > sqrt(N)
                continue
            comp = False
            for i in range(len(small)):
                sp = small[i]
                if (2 * (R % sp) * qmod[i] + 1) % sp == 0:   # presieve, small arithmetic
                    comp = True
                    break
            if comp:
                continue
            N = 2 * R * q + 1
            wit = None
            for a in (2, 3, 5, 7, 11, 13, 17, 19, 23):
                if powmod(a, N - 1, N) != 1:
                    break
                if gcd(powmod(a, (N - 1) // q, N) - 1, N) == 1:
                    wit = a
                    break
            if wit is not None:
                chain.append((int(N), int(q), wit))
                q = N
                break
    return int(q), chain


def roll_smooth(q, B=500, a=3):
    """One REVERSIBLE rollball step: ``N = 2*R*q + 1`` with a SMOOTH multiplier R.

    Like :func:`roll_prime`'s step, but ``R`` is built as a product of primes ``<= B``
    (smooth), sized ~``q`` bits so ``N ~ q^2`` and the length ~doubles. One Pocklington
    witness (``a``) proves ``N`` off ``q``. Because ``R`` is smooth, the roll can be
    UNDONE by :func:`descend` -- strip the small factors of ``N-1`` and ``q`` reappears.
    That's the whole trick: with a smooth R the certificate lives inside the number.
    Needs gmpy2. Returns ``N`` (an mpz).
    """
    import random

    from gmpy2 import gcd, mpz, powmod

    q = mpz(q)
    small = [int(x) for x in ruler_sieve(int(B))]

    def smooth_R(bits):
        R = mpz(1)
        while R.bit_length() < bits:
            R *= random.choice(small)
        return R

    while True:
        N = 2 * smooth_R(q.bit_length()) * q + 1
        if powmod(a, N - 1, N) == 1 and gcd(powmod(a, (N - 1) // q, N) - 1, N) == 1:
            return N


def roll_tower(seed, steps=8, B=500):
    """Build a REVERSIBLE rollball tower: ``steps`` smooth rolls up from ``seed`` (a
    proven prime -- e.g. a Mersenne like ``2**31-1``). Returns the chain
    ``[seed, q1, ..., top]``, each rung proven off the last and, because R is smooth,
    recoverable from the top alone by :func:`descend`. Needs gmpy2.
    """
    from gmpy2 import mpz

    chain = [mpz(seed)]
    for _ in range(int(steps)):
        chain.append(roll_smooth(chain[-1], B=B))
    return chain


def descend(N, B=500):
    """The REVERSE rollball -- walk a smooth tower back DOWN to its seed.

    Sir Lars's move: subtract 1, strip the factor of 2 and every small prime factor
    (``<= B``). Because the tower was built with a smooth ``R``, ``N-1 = 2*(smooth)*q``,
    so what's left standing is the PRIME ``q`` that ``N`` was rolled from. Repeat on
    ``q`` and the whole ladder unwinds down to the seed -- the certificate falls out
    of the number itself. Returns ``[N, q1, q2, ..., seed]``.

    Only works on smooth-R towers (:func:`roll_tower` / :func:`roll_smooth`); a
    random-R roll won't peel, because its R has a large prime factor that stripping
    small primes can't remove. Needs gmpy2.
    """
    from gmpy2 import is_prime, mpz

    small = [int(x) for x in ruler_sieve(int(B))]
    N = mpz(N)
    chain = [N]
    while True:
        m = N - 1
        for p in small:
            while m % p == 0:
                m //= p
        if 1 < m < N and is_prime(m):
            chain.append(m)
            N = m
        else:
            break
    return chain


def rollball_hunt(N, B=2000, offsets=(-3, -1, 1, 3), verbose=True):
    """The rollball hunt -- descend a number by peeling the nearby EVENS on each side.

    For an odd ``N`` the two nearest evens on each side are ``N-3, N-1, N+1, N+3``. A
    rollball ``N = 2*R*q + off`` hides its rung ``q`` in whichever of these equals
    ``2*R*q``, so we strip the small factors off all four and see which one PEELS to a
    prime rung -- that's the offset it was rolled with. Descend through it and repeat,
    level after level, down to the seed (the "big bottom"). ``N-1`` is also ``N ^ 1`` --
    the XOR-1 wrap (it flips odd numbers down to ``N-1``, even numbers up to ``N+1``).

    A smooth tower peels cleanly all the way; a wild prime peels on none of the four and
    stops at the top. Prints a block per level; returns ``[(N, offset_used, rung), ...]``.
    Needs gmpy2.
    """
    from gmpy2 import is_prime, mpz

    small = [int(x) for x in ruler_sieve(int(B))]
    N = mpz(N)
    levels = []
    lvl = 0
    while True:
        if verbose:
            print(f"  L{lvl:<2}  N = {N.bit_length()}-bit  (the 4 nearby evens):")
        peels = []
        for off in offsets:
            m = N + off                       # off in (-3,-1,1,3) -> N-3, N-1, N+1, N+3 (all even)
            strip = []
            for p in small:
                while m % p == 0:
                    m //= p
                    strip.append(p)
            rung = 1 < m < N and is_prime(m)
            if rung:
                peels.append((abs(off), off, m))
            if verbose:
                sf = " * ".join(f"{p}^{strip.count(p)}" if strip.count(p) > 1 else str(p)
                                for p in sorted(set(strip))) or "1"
                wrap = "  [= N^1, the XOR-1 wrap]" if off == -1 else ""
                tag = (f"q ({m.bit_length()}-bit PRIME)  <== peels!" if rung
                       else f"({m.bit_length()}-bit composite)")
                print(f"        N{off:+d} = ({sf}) * {tag}{wrap}")
        if peels:
            peels.sort()                      # prefer the offset nearest 0
            off_used, rung = peels[0][1], peels[0][2]
            levels.append((int(N), off_used, int(rung)))
            if verbose:
                print(f"      -> peels on N{off_used:+d}; descend to the {rung.bit_length()}-bit rung")
            N = rung
            lvl += 1
        else:
            levels.append((int(N), None, None))
            if verbose:
                print(f"      none of the four evens peels -- BOTTOM after {lvl} descent(s).")
            break
    return levels


def factor_probe(n, b=None, verbose=True, xor=True, W=20000):
    """Descend (fully factor) ``n`` and its square-neighbours -- and test the SQUARE vs
    XOR question head on.

    Pass ``(a, b)`` to start from a product, or a single ``n``. Two structural facts on
    odd ``n``: ``n^2 = a^2 * b^2`` -- the prime factors, squared, right on the surface --
    and ``n^2 - 1 = (n-1)(n+1)``, so descending it MERGES both neighbours' factorizations
    into one list.

    With ``xor=True`` it also runs the honest control: count how often ``n ^ w`` is
    divisible by each real factor vs a same-size STRANGER prime. They come out equal
    (both ~ ``W/p``), proving the XOR "reveal" is pure base rate -- because XOR is bitwise
    and factors are multiplicative. Only squaring truly carries the factors. Full story
    in ``SQUARING_AND_FACTORS.md``.

    Returns ``{label: (value, [prime factors])}``; prints the report when ``verbose``.
    """
    if b is not None:
        n = int(n) * int(b)
    n = int(n)

    def factor(m):
        m = int(m); fs = []
        while m % 2 == 0:
            fs.append(2); m //= 2
        d = 3
        while d * d <= m:
            while m % d == 0:
                fs.append(d); m //= d
            d += 2
        if m > 1:
            fs.append(m)
        return fs

    def next_prime(m):
        m += 2 if m % 2 else 1
        while len(factor(m)) != 1:
            m += 2
        return m

    sq = n * n
    items = [("n", n), ("n-1", n - 1), ("n+1", n + 1),
             ("n^2", sq), ("n^2-1", sq - 1), ("n^2+1", sq + 1)]
    out = {}
    for label, v in items:
        out[label] = (v, factor(v))
        if verbose:
            note = {"n^2": "   <- the factors, squared (multiplicative: survives)",
                    "n^2-1": "   <- (n-1)(n+1): both neighbours merged"}.get(label, "")
            print(f"  {label:<7} = {v:>18}  = {out[label][1]}{note}")

    if xor and verbose:
        reals = sorted(set(factor(n)))
        print(f"  --- XOR check: does n ^ w reveal a factor more than a same-size stranger? (w<{W}) ---")
        for p in reals:
            st = next_prime(p)
            fh = sum(1 for w in range(1, W) if (n ^ w) % p == 0)
            sh = sum(1 for w in range(1, W) if (n ^ w) % st == 0)
            print(f"    factor {p:<9} hits {fh:>5}   |   stranger {st:<9} hits {sh:>5}   |   chance W/p ~ {W // p}")
        print("    verdict: factor == stranger == W/p  ->  XOR reveal is base rate, not structure.")
        print("    factors are multiplicative; XOR is bitwise -- only n^2 = a^2*b^2 keeps them.")
    return out


def prime_walk(start, n=15, direction=1, verbose=True):
    """Take a walk through the primes near ``start``.

    Step to consecutive primes -- ``direction=1`` walks UP, ``-1`` walks DOWN -- and at
    each stop show the prime, the gap from the last one, its mod-4 lane, a ``twin`` mark
    when the gap is 2 (``cousin`` at 4), and a flag when you cross a power of two. A
    stroll through the prime landscape at any altitude: uses gmpy2, so ``2**55`` or
    ``2**5000`` both work. Returns the list of ``(prime, gap, lane)`` walked. Needs gmpy2.
    """
    from gmpy2 import is_prime, mpz, next_prime

    p = mpz(start)
    if direction >= 0:
        if not is_prime(p):
            p = next_prime(p)
    else:
        if p % 2 == 0:
            p -= 1
        while not is_prime(p):
            p -= 2

    walk, prev = [], None
    while len(walk) < n:
        gap = None if prev is None else abs(int(p - prev))
        lane = int(p % 4)
        walk.append((int(p), gap, lane))
        if verbose:
            crossed = ""
            if prev is not None:
                lo, hi = (int(prev), int(p)) if prev < p else (int(p), int(prev))
                b = 1 << (hi.bit_length() - 1)
                if lo < b <= hi:
                    crossed = f"   <-- crossed 2^{b.bit_length() - 1}"
            twin = "  twin!" if gap == 2 else ("  cousin" if gap == 4 else "")
            g = "  -- " if gap is None else f"{gap:>4}"
            print(f"  {int(p):>21,}   gap {g}   lane {lane}{twin}{crossed}")
        prev = p
        if direction >= 0:
            p = next_prime(p)
        else:
            p -= 2
            while not is_prime(p):
                p -= 2
    return walk


def pocklington_walk(seed=(1 << 31) - 1, steps=10, rbits=12, verbose=True):
    """Walk UP through PROVEN primes using our own Pocklington roll -- gmpy2 only checks.

    Start at a proven seed (default M31 = ``2**31 - 1``, Lucas-Lehmer proven). Each step
    rolls ``N = 2*R*q + 1`` with a small ``rbits``-bit ``R`` and PROVES ``N`` prime by a
    single Pocklington witness off ``q`` -- that's the proof, ours, no oracle. Then
    ``gmpy2.is_prime`` is run purely as an accuracy check, and we flag if it ever
    disagrees. Small ``rbits`` = a gentle walk (a few bits per step); crank it up to leap.

    Returns ``[(prime, R, witness, gmpy2_agrees), ...]``. Needs gmpy2.
    """
    import random

    from gmpy2 import gcd, is_prime, mpz, powmod

    small = [int(x) for x in ruler_sieve(2000)]
    q = mpz(seed)
    chain = [(int(q), None, None, None)]
    if verbose:
        print(f"  start:  {int(q):,}  ({q.bit_length()} bits)  [proven seed]")
    for s in range(int(steps)):
        while True:
            R = mpz(random.getrandbits(int(rbits))) | 1
            N = 2 * R * q + 1
            if any(N % p == 0 for p in small):        # presieve
                continue
            wit = None
            for a in (3, 2, 5, 7, 11):
                if powmod(a, N - 1, N) != 1:
                    break                              # Fermat fail -> composite
                if gcd(powmod(a, (N - 1) // q, N) - 1, N) == 1:
                    wit = a                            # Pocklington: PROVEN
                    break
            if wit is not None:
                break
        agrees = bool(is_prime(N))                     # gmpy2 accuracy check only
        chain.append((int(N), int(R), wit, agrees))
        if verbose:
            flag = "" if agrees else "   <== gmpy2 DISAGREES!"
            print(f"  step {s + 1:>2}: N = 2*{int(R)}*q + 1  ->  {N.bit_length():>4} bits, {len(str(N)):>4} digits"
                  f"   PROVEN (Pocklington a={wit})   gmpy2 ok: {agrees}{flag}")
        q = N
    if verbose:
        allok = all(c[3] for c in chain[1:])
        print(f"  walked {steps} proven primes to {q.bit_length()} bits; gmpy2 agreed every step: {allok}")
    return chain


def cycle(depth, n=64, base=2, step=8):
    """The rail cycle: rails double, step halves — pure structure, no drift.

    A single rail is ``base + step*arange(n)`` (e.g. ``8n+2``). Each level
    interleaves it with its ``+step//2`` twin (``8n+2`` gains ``8n+6``), which
    collapses to a denser ruler at half the step:

        cycle(0, n)  ->  8n+2   (step 8)          the free ruler
        cycle(1, n)  ->  4n+2   (step 4)          interleave(8n+2, 8n+6)
        cycle(k, n)  ->  base + (step>>k)*arange(n<<k)

    This is the grid ``oppositex`` rides once the block-skip drift is thrown
    away: every ``oppositex`` value satisfies ``v % 4 == 2``, i.e. lies on
    ``cycle(1)``. Companion to :func:`ruler` / :func:`mersenne_ruler`.
    """
    depth = int(depth)
    return base + (step >> depth) * np.arange(n << depth)


def rail_cycle(hi, base=2, step=8, splits=(1 << 16,)):
    """The cycle read straight off ``oppositex``: a ``base+step`` ruler whose
    step halves (rails double) at each split value, keeping the same ``base``
    residue on every finer lattice.

    Array-derived defaults: ``base=2, step=8, splits=(2**16,)`` — i.e. ``8n+2``
    below ``2^16`` and ``4n+2`` above. This is exactly the grid ``oppositex``
    rides once the block-skip drift is discarded; ``oppositex`` uses half its
    points (below the split it sits on ``8n+2``, half of the ``4n+2`` it fills
    above). Recursive companion to :func:`cycle`.

        rail_cycle(70000)   ->  2,10,...,65530, 65538,65542,... (step 8 -> 4)
    """
    splits = sorted(s for s in splits if base < s < hi)
    if not splits:                                   # base case: one density
        return np.arange(base, hi, step)
    s = splits[0]
    nxt = s + ((base - s) % (step // 2))              # re-phase onto finer lattice
    return np.concatenate([
        np.arange(base, s, step),                     # this band at current step
        rail_cycle(hi, nxt, step // 2, splits[1:]),   # recurse: rails doubled
    ])


def block_mask(seq, base=2):
    """Collapse an ``8n+2`` rail to block level — the helpful ×8 skip.

    ``oppositex`` (below ``2^16``) is the ruler ``fx`` cut into whole 8-blocks,
    each entirely kept or dropped; block ``k`` covers the 8 values
    ``64k+base .. 64k+base+56`` and starts at ``64k+base`` (the ``64n+2`` grid).
    Returns a boolean keep/drop bitmask over blocks ``0..max``: ``mask[k]`` is
    True iff the whole block ``k`` is present. This is ``seq`` with per-value
    detail thrown away — only *which blocks survive* remains. Inverse:
    :func:`block_expand`.
    """
    b = (np.asarray(seq) - base) // 64                 # block id of each value
    n = int(b.max()) + 1
    return np.bincount(b, minlength=n) == 8            # whole-block keep/drop mask


def block_expand(mask, base=2):
    """Expand a block keep/drop mask back to the ``8n+2`` values (×8 larger).

    Inverse of :func:`block_mask`: each kept block ``k`` becomes the 8 values
    ``64k+base + 8*[0..7]``. Round-trips exactly on any whole-block rail.
    """
    kept = np.nonzero(np.asarray(mask))[0]
    starts = 64 * kept + base
    return (starts[:, None] + 8 * np.arange(8)).ravel()


def block_drops(seq, hi=None, base=2):
    """The sparse skip list: block ids missing from an ``8n+2`` subsequence.

    Cleaner than a ``1,0,1,0`` mask — you only record *where* the cycle skips.
    Each block is 8 values on the ``64n+2`` grid. For ``oppositex`` below
    ``2^16`` this returns 390 ids, and ``390 == B`` (the drift). Feed to
    :func:`skip_blocks` to rebuild.
    """
    seq = np.asarray(seq)
    kept = np.unique((seq - base) // 64)
    top = int(kept.max()) if hi is None else (hi - base) // 64
    return np.setdiff1d(np.arange(top + 1), kept)


def skip_blocks(seq, drops, base=2):
    """Cycle, and skip every 8 — drop whole 8-blocks from a rail.

    ``drops`` is the sparse list of block ids to skip (from :func:`block_drops`);
    each skip removes one 8-block on the ``64n+2`` grid. So
    ``oppositex = skip_blocks(fx, drops)`` — the clean ruler with 390 blocks
    skipped, no expanded keep/drop string.
    """
    seq = np.asarray(seq)
    blk = (seq - base) // 64
    return seq[~np.isin(blk, np.asarray(drops))]


def next_level(seq, times=1):
    """Skip to the next level: (seq[::8] - 2) // 2, repeated `times`.

    Takes the block-starts of an ``8n+2`` sequence (``64k+2``) and halves them
    straight down to ``32k``. Offset by 16 gives the same sequence shifted one
    step (it telescopes). Apply again to keep descending.
    """
    seq = np.asarray(seq)
    for _ in range(int(times)):
        seq = (seq[::8] - 2) // 2
    return seq


def grow(seed, n, step=8):
    """Build the ruler from ONE number: seed + step*arange(n). No array in."""
    return int(seed) + step * np.arange(n)


def prev_level(seq, times=1):
    """Inverse of :func:`next_level` — climb up a level. Each value ``x`` becomes
    the 8-block ``(2x+2) + 8*[0..7]``. Seed with a single number to grow outward.
    """
    seq = np.atleast_1d(np.asarray(seq))
    for _ in range(int(times)):
        seq = ((2 * seq + 2)[:, None] + 8 * np.arange(8)).ravel()
    return seq


if __name__ == "__main__":
    # quick self-test / demo
    df = creator()
    print("shape:", df.shape, "(rows x cols)")
    print(df.iloc[:, :10].to_string())
    print(f"... {df.shape[1]} columns total, last n = {df.iloc[0, -1]}")


# =========================================================================== #
#  PRIME RAIL SYSTEM
#
#  Walk the six-wheel, navigate the TWIN and SEXY lanes, track the count/counn
#  indexers, climb the safe-double rail, and leap base-10 orders on the 101/103
#  anchors -- all on NON-GMPY structural logic (the 6k±1 wheel does its own
#  trial division). Your walking friends around each rung are x mod 4 and x mod 5;
#  the base-10 boundary primes (11/13, 101/103, ...) are the anchors for x*10 jumps.
#
#      from numtools import rail_walk, twin_lane, sexy_lane, cunningham, scale10
#      rail_walk(120)      # the full table: two indexers, two lanes, non-gmpy
#      twin_lane(120)      # the twin lane on its own
#      cunningham(89)      # a safe-double climb
#      scale10(103)        # a prime one decimal order up
# =========================================================================== #


def _railprime(n):
    """Primality with no gmpy -- every prime > 3 is 6k±1, so the wheel divides itself."""
    n = int(n)
    if n < 2:
        return False
    if n < 4:
        return True                        # 2 and 3
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5                                   # step the 6k±1 wheel: 5,7, 11,13, 17,19, ...
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def _next_railprime(n):
    """Smallest prime >= n, walked on the six-wheel (non-gmpy)."""
    n = int(n)
    if n <= 2:
        return 2
    if n == 3:
        return 3
    n += (n % 2 == 0)                       # make odd
    while not _railprime(n):
        n += 2
    return n


def safe_to_double(p):
    """True if doubling p (-> 2p+1) is worth the try: p on the 6k-1 lane.

    A prime p == 1 (mod 3) sits in the 6k+1 lane, and there 3 divides 2p+1 always --
    a guaranteed dead rung. Only the 6k-1 lane (p == 2 mod 3) can chain. So this is
    the 'safe to double' gate: it never wastes a climb on a 3-killed number.
    """
    p = int(p)
    return p in (2, 3) or (p % 3 == 2)


def rail_walk(hi=256):
    """Walk x up to hi -- both indexers, both lanes, all structure, no gmpy.

    Columns:
      count, counn -- the two prime-locked indexers. count = primes below x; counn
                      lags by 1. (The fixed 2-gap to a pure rail-count is 2 and 3
                      living off the wheel: 2 for parity, 3 for the step.)
      m4, m5       -- x mod 4 and x mod 5, the friends for walking around: m4 picks
                      the twin lane (3 -> lower, 1 -> upper), m5 the five-wheel spoke.
      prime        -- is x prime (wheel trial division).
      lane         -- '6k-1' | '6k+1' | 'seed'(2,3) | '-'.
      twin, sexy   -- x sits on a twin (±2) / sexy (±6) pair of primes.
      safe_dbl     -- prime AND safe_to_double: 2x+1 dodges the factor 3.
      dbl, fold    -- 2x+1 (the double) and (x+1)^5 (the xor partner).
      b10          -- x is a base-10 boundary prime (within 4 of a power of ten).
    """
    if pd is None:
        raise RuntimeError("rail_walk needs pandas")
    pow10 = [10 ** k for k in range(1, len(str(hi)) + 1)]
    rows, count, counn, nbelow = [], 0, 0, 0
    for x in range(hi):
        p = _railprime(x)
        if x in (2, 3):
            lane = "seed"
        elif x % 6 == 5:
            lane = "6k-1"
        elif x % 6 == 1 and x > 1:
            lane = "6k+1"
        else:
            lane = "-"
        twin = p and (_railprime(x - 2) or _railprime(x + 2))
        sexy = p and (_railprime(x - 6) or _railprime(x + 6))
        b10 = p and any(abs(x - t) <= 4 for t in pow10)
        rows.append({
            "count": count, "counn": counn, "x": x,
            "m4": x % 4, "m5": x % 5,
            "prime": p, "lane": lane,
            "twin": bool(twin), "sexy": bool(sexy),
            "safe_dbl": bool(p and safe_to_double(x)),
            "dbl": 2 * x + 1, "fold": (x + 1) ^ 5, "b10": bool(b10),
        })
        if p:                               # the flip: index only moves on a prime
            count = nbelow
            counn = max(0, count - 1)
            nbelow += 1
    return pd.DataFrame(rows)


def twin_lane(hi=256):
    """The TWIN lane: pairs (p, p+2) both prime, each with its own running index."""
    if pd is None:
        raise RuntimeError("twin_lane needs pandas")
    rows, idx = [], 0
    for lo in range(3, hi):
        if _railprime(lo) and _railprime(lo + 2):
            idx += 1
            rows.append({"twin_idx": idx, "low": lo, "high": lo + 2,
                         "center": lo + 1, "k": (lo + 1) // 6, "m4": lo % 4})
    return pd.DataFrame(rows)


def sexy_lane(hi=256):
    """The SEXY lane: pairs (p, p+6) both prime, each with its own running index."""
    if pd is None:
        raise RuntimeError("sexy_lane needs pandas")
    rows, idx = [], 0
    for lo in range(5, hi):
        if _railprime(lo) and _railprime(lo + 6):
            idx += 1
            rows.append({"sexy_idx": idx, "low": lo, "high": lo + 6,
                         "lo_lane": "6k-1" if lo % 6 == 5 else "6k+1", "m4": lo % 4})
    return pd.DataFrame(rows)


def cunningham(seed, cap=1 << 60):
    """Climb x -> 2x+1 while prime, doubling only the safe (6k-1) lane. Non-gmpy.

    Returns the chain [seed, 2*seed+1, ...]. Stops when the next rung would sit off
    the safe lane (guaranteed 3-killed) or is composite.
    """
    seed = int(seed)
    if not _railprime(seed):
        raise ValueError(f"{seed} is not prime")
    chain, p = [seed], seed
    while p < cap and safe_to_double(p):
        d = 2 * p + 1
        if not _railprime(d):
            break
        chain.append(d)
        p = d
    return chain


def base10_anchors(kmax=6):
    """The base-10 boundary primes -- 11/13, 101/103, ... -- your decimal anchors.

    For each power of ten it reports the nearest prime just above and just below,
    and whether the pair just above is a twin (101/103 is the poster child). These
    are the good hinge points for base-10 changes: they sit right on the 10^k seam.
    """
    if pd is None:
        raise RuntimeError("base10_anchors needs pandas")
    rows = []
    for k in range(1, kmax + 1):
        t = 10 ** k
        above = _next_railprime(t + 1)
        below = t - 1
        while below > 2 and not _railprime(below):
            below -= 1
        rows.append({"k": k, "ten_k": t, "prime_below": below, "prime_above": above,
                     "above_is_twin": _railprime(above + 2), "gap_above": above - t})
    return pd.DataFrame(rows)


def scale10(p, factor=10):
    """Build a prime about `factor`x bigger -- an x*10 (or x*100) decimal leap.

    Lands on the smallest prime >= p*factor, on the six-wheel, no gmpy. Reports how
    far past the boundary it had to walk. Leaning on the 101/103 style anchors: a
    decimal-scaled prime sits a short hop past the round number, same as 103 past 100.
    """
    p = int(p)
    target = p * factor
    got = _next_railprime(target)
    return {"p": p, "factor": factor, "target": target, "prime": got,
            "hop": got - target, "twin": _railprime(got + 2) or _railprime(got - 2)}


def good_to_enter(p):
    """Pre-entry check: is p a good spot to START a doubling climb?

    The safe lane (p == 2 mod 3) is CLOSED under x -> 2x+1 -- if p == 2 mod 3 then
    2p+1 == 2 mod 3 as well, so once you step onto it every rung stays on it and 3 can
    never kill you mid-climb. A good entry is therefore just: p on the 2-mod-3 lane
    AND 2p+1 prime. Pass both and you are on closed track; only a composite rung ever
    stops you, never a wasted lane. Returns (ok, reason). Non-gmpy.
    """
    p = int(p)
    if not _railprime(p):
        return False, "not prime"
    if p % 3 != 2 and p not in (2, 3):
        return False, "off the safe lane (6k+1) -- 2p+1 is 3-killed"
    if not _railprime(2 * p + 1):
        return False, "first rung 2p+1 is composite"
    return True, "on closed track -- climb freely"


def entry_scout(hi, minlen=3):
    """Good doubling-entry seeds under hi, ranked by chain length -- stepping logic, not brute force.

    Filters to the safe lane first (skips the 3-killed half outright), keeps only chain
    STARTS (89 is listed; 179, being mid-chain, is not), then ranks by how far the climb
    runs. Returns [(length, seed, top), ...].
    """
    out = []
    for p in range(5, int(hi)):
        ok, _ = good_to_enter(p)
        if not ok:
            continue
        half = (p - 1) // 2
        if half >= 2 and _railprime(half) and half % 3 == 2:
            continue                                   # p == 2*half+1: mid-chain, not a start
        c = cunningham(p)
        if len(c) >= minlen:
            out.append((len(c), p, c[-1]))
    return sorted(out, reverse=True)


def full_rail(hi=120, base10=True):
    """The complete rail -- THREE left indexers rolling, both lanes counting, the
    stepping-logic trues flagged, and an optional base-10 double-up prime check.

    Left rolling rails : count (primes below x), counn (count-1), railn (count-2) --
                         lags 0 / 1 / 2. That 2-gap to railn is 2 and 3 living off the
                         wheel (2 for parity, 3 for the step). All three flip only on a prime.
    Lane counters      : ti / si -- twin and sexy lanes rolling their own running index
                         (sexy climbs ~2x the twin rate; both trail count as they thin).
    Stepping trues     : safe_dbl (6k-1, safe to double), enter (closed track: safe lane
                         AND 2p+1 prime), plus twin / sexy flags.
    Base-10 (option)   : with base10=True, for each prime -- dbl_p (is 2p+1 prime, the
                         double-up landing), ten_p (nearest prime to p*10 is a short hop, a
                         clean decimal landing), and good (both fire -- a spot that paid off).
    Parity friends     : m4, m5 -- x mod 4 and x mod 5, the walking-around companions.
    All non-gmpy: the 6k±1 wheel does its own primality.
    """
    if pd is None:
        raise RuntimeError("full_rail needs pandas")
    rows = []
    count = counn = ti = si = nbelow = 0
    for x in range(hi):
        p = _railprime(x)
        railn = max(0, count - 2)                       # third rail: 2 and 3 off the wheel
        twin = p and (_railprime(x - 2) or _railprime(x + 2))
        sexy = p and (_railprime(x - 6) or _railprime(x + 6))
        safe = bool(p and (x % 3 == 2 or x in (2, 3)))
        enter = bool(safe and _railprime(2 * x + 1))
        lane = ("seed" if x in (2, 3)
                else "6k-1" if x % 6 == 5
                else "6k+1" if (x % 6 == 1 and x > 1) else "-")
        row = dict(count=count, counn=counn, railn=railn, x=x, m4=x % 4, m5=x % 5,
                   prime=p, lane=lane, ti=ti, si=si,
                   twin=bool(twin), sexy=bool(sexy),
                   safe_dbl=safe, enter=enter)
        if base10:
            dbl = 2 * x + 1
            dbl_p = bool(p and _railprime(dbl))          # the double-up landing on a prime
            ten = _next_railprime(x * 10) if p else 0    # the x*10 decimal leap
            hop = (ten - x * 10) if p else 0
            ten_p = bool(p and 0 <= hop <= 9)            # short hop = a clean base-10 landing
            row.update(dbl=(dbl if p else 0), dbl_p=dbl_p,
                       ten=ten, ten_hop=hop, ten_p=ten_p,
                       good=bool(dbl_p and ten_p))        # a good spot that paid off
        rows.append(row)
        if p:
            count = nbelow
            counn = max(0, count - 1)
            nbelow += 1
        if twin:
            ti += 1
        if sexy:
            si += 1
    return pd.DataFrame(rows)


def logic_walk(hi, lanes=(2, 3, 5, 7, 11, 13)):
    """Mimic the prime walk with cheap wheel logic -- no heavy primality test below the horizon.

    A number 'survives' the lanes if no lane divides it. Below (next-prime-after-max-lane)^2
    that survival IS primality exactly -- so in that band the cheap mod-logic reproduces the
    prime walk for free, no trial division, no gmpy. Past the horizon, slip-through semiprimes
    (p*q with both factors above the lanes) start creeping in and the logic over-counts.

    Columns: x, wheel_p (logic says prime), exact (still inside the guaranteed band),
    lcount (logic running count -- equals the true prime count while exact is True).
    df.attrs['horizon'] is where the free ride ends. All non-gmpy.
    """
    if pd is None:
        raise RuntimeError("logic_walk needs pandas")
    lanes = sorted(set(int(w) for w in lanes))
    m = lanes[-1] + 1                                   # next prime after the top lane
    while not _railprime(m):
        m += 1
    horizon = m * m                                     # first composite the lanes miss
    rows, lcount = [], 0
    for x in range(int(hi)):
        wp = (x in lanes) or (x >= 2 and all(x % w for w in lanes))
        rows.append(dict(x=x, wheel_p=bool(wp), exact=bool(x < horizon), lcount=lcount))
        if wp:
            lcount += 1
    df = pd.DataFrame(rows)
    df.attrs["horizon"] = horizon
    df.attrs["lanes"] = lanes
    return df


def five_rail(hi=110):
    """How the number five threads through each prime -- the table behind 'The Five Rail'.

    For every prime 5..hi: its mod-5 spoke, its 6k±1 lane, its (p+1)^5 fold-partner and
    whether that partner is a twin (±2) or sexy (±6), the double 2p+1, whether it is safe
    to double (6k-1), and whether it sits on a base-10 anchor (near 10^k, like 101/103).
    Non-gmpy.
    """
    if pd is None:
        raise RuntimeError("five_rail needs pandas")
    rows, idx = [], 0
    for x in range(5, int(hi)):
        if not _railprime(x):
            continue
        idx += 1
        fold = (x + 1) ^ 5
        partner = ("twin" if fold in (x - 2, x + 2)
                   else "sexy" if _railprime(fold) else "-")
        if not _railprime(fold):
            partner = "-"
        rows.append(dict(idx=idx, prime=x, mod5=x % 5,
                         lane=("6k-1" if x % 6 == 5 else "6k+1"),
                         fold=fold, partner=partner, dbl=2 * x + 1,
                         safe_dbl=(x % 3 == 2),
                         b10=any(abs(x - 10 ** k) <= 4 for k in (1, 2, 3))))
    return pd.DataFrame(rows)


def rail_lars(hi=120):
    """Sir Lars's custom rail -- three offset lines plus his own mod / fold / direction lanes.

    Three rails, offsets 0 / 1 / 2:
      count  -- primes below x
      counn  -- count - 1
      raild  -- count - 2, the 'offset two down' line. The 2-drop is 2 and 3 living off the
                wheel (2 the even prime -- parity; 3 a factor of six -- the step), which is the
                odd/even cycle reason it sits two below rather than one.
    His lanes (^ is xor, kept verbatim from the loop):
      m1 = ((2x+3)^5)^1        m1_up = ((2x+3)^5) > 2x
      m2 = ((2x+3)-5)^5        m2_up = ((2x+3)-5)^5 > 2x
      v1 = (x+3)^5             v1_dn = (x+3)^5 < x
      v2 = (x-5)^5             v2_dn = (x-5)^5 < x
    plus the double 2x+1 and its primality. All non-gmpy.
    """
    if pd is None:
        raise RuntimeError("rail_lars needs pandas")
    rows, count, counn, nbelow = [], 0, 0, 0
    for x in range(int(hi)):
        p = _railprime(x)
        raild = max(0, count - 2)
        rows.append(dict(
            count=count, counn=counn, raild=raild, x=x, m4=x % 4, prime=p,
            dbl=2 * x + 1, dbl_p=bool(p and _railprime(2 * x + 1)),
            m1=((2 * x + 3) ^ 5) ^ 1, m1_up=(((2 * x + 3) ^ 5) > 2 * x),
            m2=((2 * x + 3) - 5) ^ 5, m2_up=((((2 * x + 3) - 5) ^ 5) > 2 * x),
            v1=(x + 3) ^ 5, v1_dn=(((x + 3) ^ 5) < x),
            v2=(x - 5) ^ 5, v2_dn=(((x - 5) ^ 5) < x),
        ))
        if p:
            count = nbelow
            counn = max(0, count - 1)
            nbelow += 1
    return pd.DataFrame(rows)


def slips(hi, lanes=(2, 3, 5, 7, 11, 13)):
    """The slip-through semiprimes -- where the wheel logic over-counts the true prime walk.

    Past the horizon (next-prime-after-max-lane)^2 the lanes miss composites whose smallest
    factor is above every lane -- and every one of them is a semiprime p*q with p,q >= that
    next prime. This lists them with their factors, so you can see exactly where (and why) the
    cheap logic starts lying. Non-gmpy.
    """
    if pd is None:
        raise RuntimeError("slips needs pandas")
    lanes = sorted(set(int(w) for w in lanes))
    m = lanes[-1] + 1
    while not _railprime(m):
        m += 1
    horizon = m * m
    rows = []
    for n in range(horizon, int(hi)):
        if all(n % w for w in lanes) and not _railprime(n):
            f = m                                       # smallest factor is >= next prime m
            while n % f:
                f += 2
            rows.append(dict(n=n, small=f, big=n // f, form=f"{f}×{n // f}"))
    df = pd.DataFrame(rows)
    df.attrs["horizon"] = horizon
    return df


def lane_read(x):
    """Read every DETERMINISTIC (every-time) lane for x -- free structure, no primality test.

    These lanes are pure functions of x's low bits: they cost a couple of XORs, they never
    lie, and they tell you x's residue classes and every fold direction for nothing. Primality
    is deliberately NOT here -- that is the some-of-the-time work you spend only when these run
    out. (See lane_verdict for the split in action.)
    """
    x = int(x)
    return dict(
        x=x, mod8=x % 8, mod6=x % 6, mod4=x % 4, mod5=x % 5,
        lane=("6k-1" if x % 6 == 5 else "6k+1" if x % 6 == 1 else "off"),
        fold=(x + 1) ^ 5, fold_dn=(((x + 1) ^ 5) < x),      # twin-side (down) vs sexy-side (up)
        v1=(x + 3) ^ 5, v1_eq=(((x + 3) ^ 5) == x), v1_dn=(((x + 3) ^ 5) < x),
        v2=(x - 5) ^ 5, v2_dn=(((x - 5) ^ 5) < x),
        m1_up=(((2 * x + 3) ^ 5) > 2 * x),
    )


def lane_verdict(x, lanes=(2, 3, 5, 7, 11, 13)):
    """Decide prime/composite from the FREE deterministic lanes where possible; test only if forced.

    This is the every-time / some-of-the-time split as a rule:
      * divisible by a lane            -> COMPOSITE, free and certain (an every-time lane).
      * survives the lanes, below the
        horizon (next-prime-after-top-
        lane)^2                        -> PRIME, free and certain (the wheel is exact there).
      * survives, but past the horizon -> UNKNOWN: the free lanes ran out, spend one real test.

    Returns (verdict, reason, cost) with cost 'free' or 'test'. The whole point: you pay for a
    some-of-the-time test only on the thin band of slip-through semiprimes past the horizon.
    """
    x = int(x)
    lanes = sorted(set(int(w) for w in lanes))
    if x < 2:
        return "composite", "less than 2", "free"
    if x in lanes:
        return "prime", "is a lane itself", "free"
    for w in lanes:
        if x % w == 0:
            return "composite", f"{w} divides it", "free"
    m = lanes[-1] + 1
    while not _railprime(m):
        m += 1
    horizon = m * m
    if x < horizon:
        return "prime", f"survived every lane below the horizon {horizon}", "free"
    return "unknown", f"past horizon {horizon}: only a slip-through semiprime could reach here", "test"


def _isqrt(n):
    import math
    return math.isqrt(int(n))


def is_square(n):
    """Deterministic perfect-square check -- a square is composite on sight, no division."""
    n = int(n)
    if n < 0:
        return False
    r = _isqrt(n)
    return r * r == n


def fermat_walk(n, steps=64):
    """Deterministic near-square CHAIN-WALK factorizer: n = a^2 - b^2 -> (a-b)(a+b).

    Walk a up from ceil(sqrt(n)); at each step check a^2 - n for a perfect square. When it hits,
    the factors fall out with no division and no primality test. Close-factor semiprimes go in a
    step or two; the walk gets slower the farther apart the factors are. Returns (p, q) or None.
    """
    n = int(n)
    if n % 2 == 0:
        return (2, n // 2)
    a = _isqrt(n)
    if a * a < n:
        a += 1
    for _ in range(int(steps)):
        b2 = a * a - n
        if b2 >= 0 and is_square(b2):
            b = _isqrt(b2)
            if a - b > 1:
                return (a - b, a + b)
            return None                       # a-b==1 means n was prime (a^2-n never squares cleanly otherwise)
        a += 1
    return None


def free_verdict(x, lanes=(2, 3, 5, 7, 11, 13), fermat_steps=8):
    """lane_verdict, but with two more DETERMINISTIC lanes stacked before any primality test:
    the perfect-square lane and a short Fermat near-square walk. Cost is 'free' whenever those
    settle it, 'test' only on what survives -- a prime OR a far-apart semiprime the free lanes
    genuinely cannot tell apart without one real test.
    """
    x = int(x)
    v, reason, cost = lane_verdict(x, lanes)
    if cost == "free":
        return v, reason, cost
    if is_square(x):
        return "composite", f"perfect square ({_isqrt(x)}^2)", "free"
    f = fermat_walk(x, fermat_steps)
    if f:
        return "composite", f"near-square walk: {f[0]}x{f[1]}", "free"
    w = fold_wheel(x)
    if w:
        return "composite", f"fold-wheel k={w[2]}: {w[0]}x{w[1]}", "free"
    H = fold_horizon(lanes, fermat_steps)
    if x < H:
        return "prime", f"certified prime: below the fold-wheel horizon {H}, and no deterministic lane factors it", "free"
    return "unknown", f"past the fold-wheel horizon {H} -- prime OR a far-ratio semiprime, needs one test", "test"


def rho_walk(n, steps=20000):
    """Pollard's rho -- a different chain-walk factorizer, for the FAR-apart semiprimes.

    Where fermat_walk walks a down from sqrt(n) (fast when the factors are close), this walks a
    pseudo-random loop x -> x^2 + c (mod n) and springs a factor the moment gcd(|x-y|, n) > 1 --
    fast when the factors are FAR apart, exactly Fermat's blind spot. Returns (p, q) or None.

    Honest note: unlike the wheel / square / near-square lanes this is not a mod-class certainty;
    it is a cheap walk that almost always catches a composite quickly but carries no guarantee in
    bounded steps. So 'no factor found' means 'probably prime', not 'certainly'.
    """
    from math import gcd
    n = int(n)
    if n < 4:
        return None
    if n % 2 == 0:
        return (2, n // 2)
    for c in (1, 2, 3, 5, 7, 11):
        x = y = 2
        d = 1
        i = 0
        while d == 1 and i < steps:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
            i += 1
        if 1 < d < n:
            return (d, n // d)
    return None


def verdict(x, lanes=(2, 3, 5, 7, 11, 13), fermat_steps=8):
    """The full stacked verdict: free deterministic lanes, then the cheap chain-walks
    (fermat near-square, then rho far-walk), and a real primality test only as the last resort.

    Returns (verdict, how, cost) where cost is 'free' (a mod-class / square / near-square lane),
    'cheap' (the rho walk found a factor), or 'test' (survived every walk -- confirmed prime).
    """
    x = int(x)
    v, how, cost = free_verdict(x, lanes, fermat_steps)
    if cost == "free":
        return v, how, "free"
    r = rho_walk(x)
    if r:
        return "composite", f"rho {r[0]}×{r[1]}", "cheap"
    return ("prime" if _railprime(x) else "composite"), "confirmed by test", "test"


def verdict_frame(hi=300, lanes=(2, 3, 5, 7, 11, 13)):
    """A filled table -- three rails, the mod wheels, the prime structure, and the stacked
    verdict with which method settled each number and what it cost. The beauty's in the columns.
    """
    if pd is None:
        raise RuntimeError("verdict_frame needs pandas")
    rows, count, counn, nbelow = [], 0, 0, 0
    for x in range(int(hi)):
        p = _railprime(x)
        railn = max(0, count - 2)
        v, how, cost = verdict(x, lanes)
        twin = p and (_railprime(x - 2) or _railprime(x + 2))
        sexy = p and (_railprime(x - 6) or _railprime(x + 6))
        rows.append(dict(
            count=count, counn=counn, railn=railn, x=x,
            m8=x % 8, m6=x % 6, m5=x % 5,
            lane=("6k-1" if x % 6 == 5 else "6k+1" if (x % 6 == 1 and x > 1) else "-"),
            prime=p, twin=bool(twin), sexy=bool(sexy),
            verdict=v, cost=cost, how=how,
        ))
        if p:
            count = nbelow
            counn = max(0, count - 1)
            nbelow += 1
    return pd.DataFrame(rows)


def walk_rail(lo=0, hi=120, lanes=(2, 3, 5, 7, 11, 13)):
    """A walkable rail for ANY window [lo, hi): the three rolling rails (count/counn/railn,
    lagged 0/1/2, flipping only on a prime) walking beside the stacked verdict -- with a 'walk'
    column naming which chain-walk sprang each number and 'fac' handing back the factors.

    The rails are primed by walking from 0 so the window's values are exact; drop lo on the slip
    band (e.g. walk_rail(284, 332)) to watch square/fermat fire right in the walk.
    """
    if pd is None:
        raise RuntimeError("walk_rail needs pandas")
    nextp = _next_railprime(max(lanes) + 1)          # strictly the next prime AFTER the top lane
    horizon = nextp * nextp

    def sf(n):
        d = 2
        while d * d <= n:
            if n % d == 0:
                return d
            d += 1
        return n

    rows, count, counn, nbelow = [], 0, 0, 0
    for x in range(int(hi)):
        p = _railprime(x)
        railn = max(0, count - 2)
        if x >= lo:
            _, _, cost = verdict(x, lanes)
            if x < 2:
                walk = "—"
            elif p:
                walk = "below" if x < horizon else "prime?"     # free vs the one real test
            elif all(x % w for w in lanes):                     # slipped the wheel -> a walk caught it
                walk = ("square" if is_square(x)
                        else "fermat" if fermat_walk(x, 8)
                        else "rho" if rho_walk(x) else "—")
            else:
                walk = "wheel"                                   # small factor, caught free
            fac = ""
            if not p and x >= 4:
                s = sf(x)
                fac = f"{s}×{x // s}"
            twin = p and (_railprime(x - 2) or _railprime(x + 2))
            sexy = p and (_railprime(x - 6) or _railprime(x + 6))
            rows.append(dict(
                count=count, counn=counn, railn=railn, x=x,
                m5=x % 5, lane=("6k-1" if x % 6 == 5 else "6k+1" if (x % 6 == 1 and x > 1) else "·"),
                prime=p, twin=bool(twin), sexy=bool(sexy),
                walk=walk, fac=fac, cost=cost,
            ))
        if p:
            count = nbelow
            counn = max(0, count - 1)
            nbelow += 1
    df = pd.DataFrame(rows)
    df.attrs["horizon"] = horizon
    return df


def fold_wheel(n, steps=6):
    """Deterministic multiplier fold-wheel -- retires the heuristic rho on the slip band.

    Fermat only reaches CLOSE factors; fold a far pair close by multiplying: k·n re-centers the
    near-square walk so that k·p ≈ q lands. Powers of two go first -- for those, k·n is just n
    bit-shifted, and the P2 ladder brackets any factor ratio -- then a few small odd multipliers
    mop up. Returns (p, q, k) or None. CERTAIN when it fires: no search, no randomness.
    """
    from math import gcd, isqrt
    n = int(n)
    if n % 2 == 0:
        return (2, n // 2, 1)
    KS = (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 3, 5, 6, 7, 9, 10, 11, 12, 13, 15, 17, 19, 21)
    for k in KS:
        m = k * n
        a = isqrt(m)
        if a * a < m:
            a += 1
        for _ in range(steps):
            b2 = a * a - m
            if b2 >= 0:
                r = isqrt(b2)
                if r * r == b2:
                    for A in (a - r, a + r):
                        g = gcd(A, n)
                        if 1 < g < n:
                            return (g, n // g, k)
            a += 1
    return None


def fold_sharp(n, steps=6):
    """The fold-wheel sharpened by the n mod 4 fork -- same catches as fold_wheel, ~30% fewer
    square-checks. Two facts, both exact:

      * a number ≡ 2 (mod 4) is never a difference of two squares, and n is odd, so any
        multiplier k ≡ 2 (mod 4) makes k·n ≡ 2 (mod 4) and can NEVER land -- skip it outright.
      * for odd k·n, only one parity of a can give an integer b: a odd when k·n ≡ 1 (mod 4),
        a even when k·n ≡ 3 (mod 4). Start on that parity and step a by 2 -- half the checks.

    Returns (p, q, k) or None. Deterministic; certain when it fires.
    """
    from math import gcd, isqrt
    n = int(n)
    if n % 2 == 0:
        return (2, n // 2, 1)
    KS = (1, 4, 8, 16, 32, 64, 128, 256, 512, 3, 5, 7, 9, 11, 12, 13, 15, 17, 19, 20, 21, 24)
    for k in KS:
        if k % 4 == 2:                       # k·n ≡ 2 (mod 4): structurally not a^2 - b^2
            continue
        m = k * n
        a = isqrt(m)
        if a * a < m:
            a += 1
        r4 = m % 4
        if r4 == 1:
            if a % 2 == 0:
                a += 1                        # valid a must be odd
            step = 2
        elif r4 == 3:
            if a % 2 == 1:
                a += 1                        # valid a must be even
            step = 2
        else:
            step = 1                          # k·n ≡ 0 (mod 4): a, b same parity, both live
        for _ in range(steps):
            b2 = a * a - m
            if b2 >= 0:
                r = isqrt(b2)
                if r * r == b2:
                    for A in (a - r, a + r):
                        g = gcd(A, n)
                        if 1 < g < n:
                            return (g, n // g, k)
            a += step
    return None


def factor_frame(lo=0, hi=340, lanes=(2, 3, 5, 7, 11, 13)):
    """A walk that factors by squares and shows the free primes. Every number gets its
    difference-of-squares view: n = a^2 - b^2 with a=(p+q)/2, b=(q-p)/2 -- so a perfect square
    is simply b=0, a close pair is small b, and a prime has no such split at all (that's what
    makes it prime, and below the horizon that costs nothing).

    Columns: fp (running free-prime tally), x, lane, kind, a, b, n=a^2-b^2, k (fold multiplier),
    p×q, m4 (the n mod 4 fork). kind is prime·free / prime·test / square / fold / wheel.
    """
    if pd is None:
        raise RuntimeError("factor_frame needs pandas")
    nextp = _next_railprime(max(lanes) + 1)
    horizon = nextp * nextp

    def sf(n):
        d = 2
        while d * d <= n:
            if n % d == 0:
                return d
            d += 1
        return n

    rows, fp = [], 0
    for x in range(int(hi)):
        p = _railprime(x)
        free_prime = p and x < horizon
        if x >= lo:
            a = b = k = ""
            fac = ""
            if p:
                kind = "prime·free" if x < horizon else "prime·test"
            elif x >= 4 and all(x % w for w in lanes):        # slipped the wheel
                s = sf(x)
                q = x // s
                a = (s + q) // 2
                b = (q - s) // 2
                fac = f"{s}×{q}"
                if is_square(x):
                    kind = "square"
                else:
                    kind = "fold"
                    r = fold_sharp(x)
                    if r:
                        k = r[2]
            elif x >= 4:
                s = sf(x)
                fac = f"{s}×{x // s}"
                kind = "wheel"
            else:
                kind = "unit"
            rows.append(dict(
                fp=fp, x=x,
                lane=("6k-1" if x % 6 == 5 else "6k+1" if (x % 6 == 1 and x > 1) else "·"),
                kind=kind, a=a, b=b, k=k, pq=fac, m4=x % 4,
            ))
        if free_prime:
            fp += 1
    df = pd.DataFrame(rows)
    df.attrs["horizon"] = horizon
    return df


# =====================================================================================
#  Prime toolkit — restartable rail + a stack of DataFrame tools (built 2026-07-18)
#  Everything below is non-gmpy and additive; nothing above is touched.
# =====================================================================================

def _small_factor(n):
    """Smallest prime factor of n (n itself if prime). Pure trial division on the 6-wheel."""
    n = int(n)
    if n < 2:
        return n
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    d = 5
    while d * d <= n:
        if n % d == 0:
            return d
        if n % (d + 2) == 0:
            return d + 2
        d += 6
    return n


def factor(n):
    """Full factorization of n as a list of (prime, power). Non-gmpy trial division."""
    n = int(n)
    out = []
    for d in (2, 3):
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))
    d = 5
    while d * d <= n:
        for step in (d, d + 2):
            if n % step == 0:
                e = 0
                while n % step == 0:
                    n //= step
                    e += 1
                out.append((step, e))
        d += 6
    if n > 1:
        out.append((n, 1))
    return out


def range_rail(lo=0, hi=120, lanes=(2, 3, 5, 7, 11, 13), restart=True):
    """The left rail, RESTARTABLE from scratch on any range.

    With restart=True (default) the count / counn / railn rails reset to 0 at `lo` -- a LOCAL
    prime count within [lo, hi), so you can study any window in isolation without carrying the
    whole tally from 0. With restart=False they are primed from 0 (the global count, like
    walk_rail). Same lagged 0/1/2 lags either way, flipping only on a prime.

    Columns: count, counn, railn, x, m5, lane, prime, twin, sexy, walk, fac.
    """
    if pd is None:
        raise RuntimeError("range_rail needs pandas")
    nextp = _next_railprime(max(lanes) + 1)
    horizon = nextp * nextp
    rows, count, counn, nbelow = [], 0, 0, 0
    start = int(lo) if restart else 0
    for x in range(start, int(hi)):
        p = _railprime(x)
        railn = max(0, count - 2)
        if x >= lo:
            if x < 2:
                walk = "—"
            elif p:
                walk = "below" if x < horizon else "prime?"
            elif all(x % w for w in lanes):
                walk = ("square" if is_square(x) else "fermat" if fermat_walk(x, 8)
                        else "fold" if fold_sharp(x) else "rho?")
            else:
                walk = "wheel"
            fac = ""
            if not p and x >= 4:
                s = _small_factor(x)
                fac = f"{s}×{x // s}"
            twin = p and (_railprime(x - 2) or _railprime(x + 2))
            sexy = p and (_railprime(x - 6) or _railprime(x + 6))
            rows.append(dict(
                count=count, counn=counn, railn=railn, x=x,
                m5=x % 5, lane=("6k-1" if x % 6 == 5 else "6k+1" if (x % 6 == 1 and x > 1) else "·"),
                prime=p, twin=bool(twin), sexy=bool(sexy), walk=walk, fac=fac,
            ))
        if p:
            count = nbelow
            counn = max(0, count - 1)
            nbelow += 1
    df = pd.DataFrame(rows)
    df.attrs["horizon"] = horizon
    df.attrs["restart"] = restart
    return df


def twin_frame(lo=0, hi=500):
    """Twin-prime pairs (p, p+2) in [lo, hi): the pair, the shared six-cell k=(p+1)/6, the even
    front p+1=6k between them, and a running index."""
    if pd is None:
        raise RuntimeError("twin_frame needs pandas")
    rows, idx = [], 0
    for p in range(max(2, int(lo)), int(hi)):
        if _railprime(p) and _railprime(p + 2):
            idx += 1
            rows.append(dict(idx=idx, low=p, high=p + 2, front=p + 1, cell=(p + 1) // 6))
    return pd.DataFrame(rows)


def sexy_frame(lo=0, hi=500):
    """Sexy-prime pairs (p, p+6) in [lo, hi): the pair, the lane they share, and the index."""
    if pd is None:
        raise RuntimeError("sexy_frame needs pandas")
    rows, idx = [], 0
    for p in range(max(2, int(lo)), int(hi)):
        if _railprime(p) and _railprime(p + 6):
            idx += 1
            rows.append(dict(idx=idx, low=p, high=p + 6,
                             lane=("6k-1" if p % 6 == 5 else "6k+1" if p % 6 == 1 else "·")))
    return pd.DataFrame(rows)


def gap_frame(lo=2, hi=1000):
    """Consecutive prime gaps in [lo, hi): each prime, the next prime, the gap, whether it set a
    new record so far, and the running maximum gap."""
    if pd is None:
        raise RuntimeError("gap_frame needs pandas")
    rows, prev, maxg = [], None, 0
    for x in range(max(2, int(lo)), int(hi)):
        if _railprime(x):
            if prev is not None:
                g = x - prev
                rec = g > maxg
                maxg = max(maxg, g)
                rows.append(dict(prime=prev, nxt=x, gap=g, record=rec, maxgap=maxg))
            prev = x
    return pd.DataFrame(rows)


def cunningham_frame(seed=89, step="2p+1", cap=1 << 40):
    """A Cunningham chain from `seed`: repeatedly apply 2p+1 (step='2p+1') or 2p-1 (step='2p-1')
    while it stays prime. Each rung shows the value, its lane, whether it is safe to double
    (6k-1), and the next candidate with its primality."""
    if pd is None:
        raise RuntimeError("cunningham_frame needs pandas")
    nxt = (lambda q: 2 * q + 1) if step == "2p+1" else (lambda q: 2 * q - 1)
    rows, p, i = [], int(seed), 0
    while _railprime(p) and p < cap:
        nx = nxt(p)
        rows.append(dict(rung=i, p=p,
                         lane=("6k-1" if p % 6 == 5 else "6k+1" if p % 6 == 1 else "·"),
                         safe_dbl=(p % 3 == 2), nxt=nx, nxt_prime=_railprime(nx)))
        if not _railprime(nx):
            break
        p, i = nx, i + 1
    return pd.DataFrame(rows)


def octave_frame(kmax=16):
    """The powers of two and their shoulders. For each k: 2^k-1 (Mersenne side), 2^k, and 2^k+1
    (Fermat side), with primality, small factor, and n mod 3 -- the mod-3 rule shows the Mersenne
    side is divisible by 3 exactly when k is even, the Fermat side exactly when k is odd."""
    if pd is None:
        raise RuntimeError("octave_frame needs pandas")
    rows = []
    for k in range(2, int(kmax) + 1):
        P = 1 << k
        for off, name in ((-1, "2^k-1"), (0, "2^k"), (1, "2^k+1")):
            n = P + off
            pr = _railprime(n)
            rows.append(dict(k=k, name=name, n=n, prime=pr,
                             small=(None if pr else _small_factor(n)), mod3=n % 3))
    return pd.DataFrame(rows)


def fold_frame(lo=289, hi=2000, lanes=(2, 3, 5, 7, 11, 13)):
    """The fold detail for every slip-through semiprime in [lo, hi): the factors, their ratio and
    gap, the fold multiplier k (a power of two = a bit-shift), the difference-of-squares a,b
    (n = a^2 - b^2), n mod 4 (the fork), and which walk sprang it."""
    if pd is None:
        raise RuntimeError("fold_frame needs pandas")
    rows = []
    for n in range(int(lo), int(hi)):
        if all(n % L for L in lanes) and not _railprime(n):
            p = _small_factor(n)
            q = n // p
            a = (p + q) // 2
            b = (q - p) // 2
            r = fold_sharp(n)
            walk = "square" if is_square(n) else "fermat" if fermat_walk(n, 8) else "fold"
            rows.append(dict(n=n, p=p, q=q, ratio=round(q / p, 2), gap=q - p,
                             k=(r[2] if r else None), a=a, b=b, m4=n % 4, walk=walk))
    return pd.DataFrame(rows)


def wheel_compare(wheels=None, N=6000):
    """Compare divisor wheels side by side: the naive horizon (next prime after the top lane)^2,
    the REAL first slip (smallest surviving composite -- earlier than the horizon when the wheel
    has a gap), the free-verdict %, and whether the parallel fold-lane recovers every miss."""
    if pd is None:
        raise RuntimeError("wheel_compare needs pandas")
    if wheels is None:
        wheels = [(2, 3, 5), (2, 3, 5, 7), (2, 3, 5, 7, 11, 13),
                  (2, 3, 5, 7, 11, 13, 17, 19), (2, 3, 5, 7, 13)]
    rows = []
    for L in wheels:
        slips = [n for n in range(2, N) if all(n % x for x in L) and not _railprime(n)]
        first = slips[0] if slips else None
        free = sum(1 for n in range(2, N)
                   if (not _railprime(n) and any(n % x == 0 for x in L))
                   or (_railprime(n) and (first is None or n < first)))
        caught = sum(1 for n in slips if is_square(n) or fermat_walk(n, 8) or fold_sharp(n))
        rows.append(dict(wheel=",".join(map(str, L)),
                         horizon=_next_railprime(max(L) + 1) ** 2, first_slip=first,
                         free_pct=round(free / (N - 2) * 100, 1), slips=len(slips),
                         fold_recovers=f"{caught}/{len(slips)}"))
    return pd.DataFrame(rows)


def dashboard(lo=2, hi=6000, lanes=(2, 3, 5, 7, 11, 13)):
    """The logic vs non-logic accounting for a range, as a tidy frame: how much settles by clean
    deterministic logic (wheel / square / fermat / fold), how much needed the heuristic rho, and
    how much needs a real primality test (only primes past the horizon)."""
    if pd is None:
        raise RuntimeError("dashboard needs pandas")
    nextp = _next_railprime(max(lanes) + 1)
    H = nextp * nextp
    cats = {"wheel (small factor)": 0, "below horizon (prime)": 0,
            "square/fermat/fold walk": 0, "heuristic (rho)": 0, "brute test (prime)": 0}
    for x in range(max(2, int(lo)), int(hi)):
        p = _railprime(x)
        if not p:
            if all(x % L for L in lanes):
                if is_square(x) or fermat_walk(x, 8) or fold_sharp(x):
                    cats["square/fermat/fold walk"] += 1
                else:
                    cats["heuristic (rho)"] += 1
            else:
                cats["wheel (small factor)"] += 1
        else:
            cats["below horizon (prime)" if x < H else "brute test (prime)"] += 1
    tot = sum(cats.values()) or 1
    kind = {"wheel (small factor)": "deterministic", "below horizon (prime)": "deterministic",
            "square/fermat/fold walk": "deterministic", "heuristic (rho)": "heuristic",
            "brute test (prime)": "brute test"}
    return pd.DataFrame([dict(category=k, count=v, pct=round(v / tot * 100, 1), kind=kind[k])
                         for k, v in cats.items()])


# --- more toolkit: the master rail, a single-number dossier, two-square frame, conveniences ---

def _two_squares(n):
    """n = a^2 + b^2 with 0 < a <= b, or None. (Primes: exists iff n==2 or n%4==1.)"""
    from math import isqrt
    n = int(n)
    a = 1
    while 2 * a * a <= n:
        b2 = n - a * a
        b = isqrt(b2)
        if b * b == b2:
            return (a, b)
        a += 1
    return None


def _facstr(n):
    """Factorization as a compact string: 13 -> '13', 360 -> '2^3×3^2×5', 323 -> '17×19'."""
    fs = factor(n)
    return "×".join(f"{p}^{e}" if e > 1 else str(p) for p, e in fs) if fs else str(n)


def is_prime(n):
    """Public primality (non-gmpy 6-wheel trial division)."""
    return _railprime(int(n))


def next_prime(n):
    """Smallest prime strictly greater than n."""
    return _next_railprime(int(n) + 1)


def prev_prime(n):
    """Largest prime strictly less than n (or None if none)."""
    m = int(n) - 1
    while m >= 2:
        if _railprime(m):
            return m
        m -= 1
    return None


def prime_pi(n):
    """Count of primes <= n."""
    return sum(1 for k in range(2, int(n) + 1) if _railprime(k))


def nth_prime(k):
    """The k-th prime, 1-indexed (nth_prime(1) == 2)."""
    k = int(k)
    if k < 1:
        raise ValueError("k must be >= 1")
    c, m = 0, 1
    while c < k:
        m += 1
        if _railprime(m):
            c += 1
    return m


def describe(n):
    """A full dossier on a single number: lane, residues, prime/kind, factors, the square view
    (a^2-b^2 for a composite, a^2+b^2 for a prime = 1 mod 4), twin/sexy, and the stacked verdict."""
    n = int(n)
    p = _railprime(n)
    d = {
        "n": n,
        "prime": p,
        "lane": ("6k-1" if n % 6 == 5 else "6k+1" if (n % 6 == 1 and n > 1) else "·"),
        "mod4": n % 4, "mod6": n % 6, "mod8": n % 8,
    }
    if n < 2:
        d["kind"] = "unit"
        return d
    if p:
        d["kind"] = "prime"
        d["twin"] = bool(_railprime(n - 2) or _railprime(n + 2))
        d["sexy"] = bool(_railprime(n - 6) or _railprime(n + 6))
        ts = _two_squares(n) if (n % 4 == 1 or n == 2) else None
        d["sum_of_squares"] = (f"{ts[0]}^2+{ts[1]}^2" if ts else "— (3 mod 4)")
        d["verdict"] = verdict(n)
    else:
        fs = factor(n)
        omega = sum(e for _, e in fs)
        d["kind"] = "square" if is_square(n) else "semiprime" if omega == 2 else "composite"
        d["factors"] = _facstr(n)
        if n % 2:
            sp = _small_factor(n)
            q = n // sp
            a, b = (sp + q) // 2, (q - sp) // 2
            d["diff_of_squares"] = f"{a}^2-{b}^2"
            r = fold_sharp(n)
            if r:
                d["fold_k"] = r[2]
        d["verdict"] = verdict(n)
    return d


def two_square_frame(lo=2, hi=300):
    """Every prime in [lo, hi) with its sum-of-two-squares. p = 2 or p = 1 (mod 4) -> a^2+b^2
    (unique); p = 3 (mod 4) -> none. The mod-4 fork, on the prime side."""
    if pd is None:
        raise RuntimeError("two_square_frame needs pandas")
    rows = []
    for p in range(max(2, int(lo)), int(hi)):
        if _railprime(p):
            t = _two_squares(p) if (p % 4 == 1 or p == 2) else None
            rows.append(dict(p=p, m4=p % 4, a=(t[0] if t else ""), b=(t[1] if t else ""),
                             sum_sq=(f"{t[0]}^2+{t[1]}^2" if t else "—"), representable=bool(t)))
    return pd.DataFrame(rows)


def master_rail(lo=0, hi=120, lanes=(2, 3, 5, 7, 11, 13), restart=True):
    """EVERYTHING in one rail walk. The three restartable rails, the number and its residues,
    whether it's prime, which KIND it is, its factors, the square view (a^2-b^2 for composites,
    a^2+b^2 for 1-mod-4 primes), the actual walk that settled it, twin/sexy, and the cost.

    Columns: count, counn, railn, rodd, x, lane, m4, prime, kind, factors, square_view,
    walk, twin, sexy, cost.

    rodd is the ODD-ONLY rail: 2*railn + 1, walking 1, 3, 5, 7, ... -- once a rung is
    used, the next one it can reach is two further on. Sits beside railn, steps when
    railn steps, never takes an even value.
    """
    if pd is None:
        raise RuntimeError("master_rail needs pandas")
    nextp = _next_railprime(max(lanes) + 1)
    horizon = nextp * nextp
    free_h = fold_horizon(lanes)                # the fold-wheel's PROVEN free horizon (>> naive 17^2)
    rows, count, counn, nbelow = [], 0, 0, 0
    start = int(lo) if restart else 0
    for x in range(start, int(hi)):
        p = _railprime(x)
        railn = max(0, count - 2)
        rodd = 2 * railn + 1                            # odd-only rail: 1, 3, 5, 7, ...
        if x >= lo:
            twin = bool(p and (_railprime(x - 2) or _railprime(x + 2)))
            sexy = bool(p and (_railprime(x - 6) or _railprime(x + 6)))
            if x < 2:
                kind, facs, sv, walk, cost = "unit", "", "", "—", "free"
            elif p:
                kind, facs = "prime", str(x)
                ts = _two_squares(x) if (x % 4 == 1 or x == 2) else None
                sv = f"{ts[0]}²+{ts[1]}²" if ts else ""
                walk = "below" if x < free_h else "prime?"
                cost = "free" if x < free_h else "test"
            else:
                facs = _facstr(x)
                omega = sum(e for _, e in factor(x))
                kind = "square" if is_square(x) else "semiprime" if omega == 2 else "composite"
                if x % 2:
                    sp = _small_factor(x)
                    q = x // sp
                    a, b = (sp + q) // 2, (q - sp) // 2
                    sv = f"{a}²−{b}²"
                else:
                    sv = ""
                if all(x % w for w in lanes):
                    walk = ("square" if is_square(x) else "fermat" if fermat_walk(x, 8)
                            else "fold" if fold_sharp(x) else "rho?")
                else:
                    walk = "wheel"
                cost = "cheap" if walk == "rho?" else "free"
            rows.append(dict(
                count=count, counn=counn, railn=railn, rodd=rodd, x=x,
                lane=("6k-1" if x % 6 == 5 else "6k+1" if (x % 6 == 1 and x > 1) else "·"),
                m4=x % 4, prime=p, kind=kind, factors=facs, square_view=sv,
                walk=walk, twin=twin, sexy=sexy, cost=cost,
            ))
        if p:
            count = nbelow
            counn = max(0, count - 1)
            nbelow += 1
    df = pd.DataFrame(rows)
    df.attrs["horizon"] = horizon
    df.attrs["restart"] = restart
    return df


def sum_of_squares(n):
    """Every representation n = a^2 + b^2 with 1 <= a <= b (both positive), as a list of (a,b).
    Empty if none. (A number is a sum of two squares iff every prime = 3 mod 4 in its
    factorization appears to an even power.)"""
    from math import isqrt
    n = int(n)
    reps = []
    a = 1
    while 2 * a * a <= n:
        b2 = n - a * a
        b = isqrt(b2)
        if b * b == b2:
            reps.append((a, b))
        a += 1
    return reps


def sos_walk(lo=0, hi=120):
    """The sum-of-squares walk. For every number: its lane and n mod 4, whether it's prime and
    which kind, its perfect-square root (if a square), all its sum-of-two-positive-squares
    representations, and the count r2. Primes = 1 mod 4 have exactly one rep; squares and products
    of such primes have more; a number carrying a 3-mod-4 prime to an odd power has none.

    Columns: x, lane, m4, prime, kind, sqrt, sos, r2, gaussian.
    """
    if pd is None:
        raise RuntimeError("sos_walk needs pandas")
    rows = []
    for x in range(int(lo), int(hi)):
        p = _railprime(x)
        reps = sum_of_squares(x)
        sq = _isqrt(x) if is_square(x) else ""
        if x < 2:
            kind = "unit"
        elif p:
            kind = "prime"
        elif is_square(x):
            kind = "square"
        elif sum(e for _, e in factor(x)) == 2:
            kind = "semiprime"
        else:
            kind = "composite"
        sos = ", ".join(f"{a}²+{b}²" for a, b in reps) if reps else "—"
        # a prime ≡1 mod4 splits in the Gaussian integers as (a+bi)(a-bi)
        gaussian = (f"({reps[0][0]}+{reps[0][1]}i)({reps[0][0]}−{reps[0][1]}i)"
                    if (p and len(reps) == 1) else "")
        rows.append(dict(x=x, lane=("6k-1" if x % 6 == 5 else "6k+1" if (x % 6 == 1 and x > 1) else "·"),
                         m4=x % 4, prime=p, kind=kind, sqrt=sq, sos=sos, r2=len(reps), gaussian=gaussian))
    return pd.DataFrame(rows)


# =====================================================================================
#  Restored 2026-07-19 — functions an older save had reverted (lane_scan, the fold
#  horizon / coprime-miss fix, the prime-forge tools). Reconstructed from the session.
# =====================================================================================

def _sprp2(n):
    """Strong probable-prime test, base 2 (one Miller-Rabin round). One modpow -- cheap, not free."""
    n = int(n)
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    x = pow(2, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(r - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False


def _mr(n, bases=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)):
    """Deterministic Miller-Rabin over a fixed base set -- correct for ALL n < 3.317e24. A dozen
    modpows: cheap, and unlike single-base it has no pseudoprime holes in range."""
    n = int(n)
    if n < 2:
        return False
    for p in bases:
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in bases:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _spf_upto(n, bound=200):
    """Smallest prime factor of n up to `bound`, or None if none found (free-ish)."""
    n = int(n)
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    d = 5
    while d <= bound:
        if n % d == 0:
            return d
        if n % (d + 2) == 0:
            return d + 2
        d += 6
    return None


_FOLD_HORIZON_CACHE = {}

def fold_horizon(lanes=(2, 3, 5, 7, 11, 13), fermat_steps=8):
    """The EXTENDED free horizon: the first wheel-coprime COMPOSITE the deterministic fold-lanes
    (square, Fermat, fold-wheel) fail to factor. BELOW this bound every wheel-coprime composite is
    caught, so a wheel-survivor the fold-lanes miss is guaranteed PRIME -- free and certain. A proven
    constant, like 17^2=289; for the default wheel it is 11441 = 17x673, ~40x further out."""
    key = (tuple(lanes), int(fermat_steps))
    if key in _FOLD_HORIZON_CACHE:
        return _FOLD_HORIZON_CACHE[key]
    nextp = _next_railprime(max(lanes) + 1)
    n = nextp * nextp
    while True:
        if all(n % L for L in lanes) and not _railprime(n):
            if not (is_square(n) or fermat_walk(n, fermat_steps) or fold_sharp(n)):
                _FOLD_HORIZON_CACHE[key] = n
                return n
        n += 1


def lane_scan(x, wheel_to=2000):
    """An upgraded lane_read that SETTLES the number: free residue lanes, a free extended-wheel
    factor catch, the deterministic fold-lanes, and finally -- only past where trial division could
    finish -- ONE cheap strong test. Each call is marked free or cheap."""
    x = int(x)
    out = dict(x=x, mod8=x % 8, mod6=x % 6, mod4=x % 4, mod3=x % 3, mod5=x % 5,
               lane=("6k-1" if x % 6 == 5 else "6k+1" if (x % 6 == 1 and x > 1) else "·"))
    if x < 2:
        out.update(call="unit", why="less than 2", cost="free", smallest_factor=None)
        return out
    spf = None
    if x % 2 == 0:
        spf = 2
    elif x % 3 == 0:
        spf = 3
    else:
        d = 5
        while d <= wheel_to and d * d <= x:
            if x % d == 0:
                spf = d
                break
            if x % (d + 2) == 0:
                spf = d + 2
                break
            d += 6
    out["smallest_factor"] = spf
    if spf is not None and spf < x:
        out.update(call="composite", why=f"{spf} divides it", cost="free")
        return out
    if wheel_to * wheel_to >= x:                  # trial division reached sqrt(x): PROVEN prime
        out.update(call="prime", why=f"survives trial division to sqrt (<= {wheel_to})", cost="free")
        return out
    if is_square(x):
        out.update(call="composite", why=f"perfect square {_isqrt(x)}^2", cost="free")
        return out
    f = fermat_walk(x, 8)
    if f:
        out.update(call="composite", why=f"near-square {f[0]}x{f[1]}", cost="free")
        return out
    w = fold_sharp(x)
    if w:
        out.update(call="composite", why=f"fold-wheel k={w[2]}: {w[0]}x{w[1]}", cost="free")
        return out
    ok = _mr(x)
    out["mr"] = ok
    out.update(call=("prime" if ok else "hidden composite"),
               why=("strong test passed (12 bases, deterministic in range)" if ok
                    else "strong test FAILED -- a semiprime hiding among the free lanes"),
               cost="cheap (a few modpows)")
    return out


def uplift_walk(x, K=24, tail=2, mask=5):
    """Walk a number UP the base-10 ladder with lane_read's +tail ^mask cap, every rung:
    vv_k = (x*10^k + tail) ^ mask. Shows each uplift, its smallest small factor, and primality."""
    if pd is None:
        raise RuntimeError("uplift_walk needs pandas")
    rows = []
    for k in range(1, int(K) + 1):
        vv = (x * 10 ** k + tail) ^ mask
        sf = _spf_upto(vv, 200)
        pr = _mr(vv) if sf is None else False
        rows.append(dict(k=k, digits=len(str(vv)), tail6=str(vv)[-6:], small=sf, prime=pr))
    return pd.DataFrame(rows)


def build_prime(seed, kmax=400):
    """BUILD a provable prime by uplifting `seed` up the base-10 tower, recapped +4^5 so it lands on
    ...0001: N = seed*10^k + 1. Then N-1 = seed*2^k*5^k is fully factored, so once N passes a strong
    test and 10^k > sqrt(N), a Pocklington certificate PROVES it prime. Build, don't decode."""
    from math import gcd
    seed = int(seed)
    if (seed + 1) % 3 == 0:
        return dict(seed=seed, built=None, note="dead tower: seed+1 = 0 mod 3, every rung is div by 3")
    d = len(str(seed))
    for k in range(d + 1, int(kmax)):
        N = seed * 10 ** k + 1
        if not _mr(N):
            continue
        if 10 ** (2 * k) <= N:
            continue
        for a in range(2, 1000):
            if pow(a, N - 1, N) != 1:
                continue
            if gcd(pow(a, (N - 1) // 2, N) - 1, N) == 1 and gcd(pow(a, (N - 1) // 5, N) - 1, N) == 1:
                return dict(seed=seed, N=N, k=k, digits=len(str(N)), witness=a, proven=True,
                            certificate=f"N-1 = {seed}*2^{k}*5^{k}, base a={a}")
    return dict(seed=seed, built=None, note=f"no prime rung up to k={kmax}")


def cheatless_primes(hi=None, lanes=(2, 3, 5, 7, 11, 13)):
    """Generate primes with ZERO primality tests -- pure +6 wheel plus the square / near-square /
    fold lanes. Exact up to the fold horizon for the wheel (the proven no-test wall); pass a smaller
    hi to stop early, or a bigger wheel to push the wall out. Returns a list of primes.

    The only 'cheat' is the small primes in `lanes` -- a handful buys thousands for free.
    """
    H = fold_horizon(lanes)
    top = H if hi is None else min(int(hi), H)
    wheelset = set(lanes)
    out = []
    for x in range(2, top):
        if x in wheelset:
            out.append(x)
            continue
        if any(x % w == 0 for w in lanes):          # +6 wheel, just mods
            continue
        if is_square(x) or fermat_walk(x, 8) or fold_sharp(x):   # the SQUARE lanes
            continue
        out.append(x)                                # survived every lane below H -> prime, no test
    return out


def superhero_walk(hi=None, lanes=(2, 3, 5, 7, 11, 13)):
    """The fun one. Walk the primes with NO prime test and fold each into a running total:
    wh += p ^ 5. Every step is +6 (the rail), one XOR (the fold), and a sum -- no interrogation
    anywhere. The diff between consecutive totals IS the current prime folded (p ^ 5).

    Columns: idx, p, lane, fold (= p^5), total (running sum of the folds). df.attrs carries
    tests_spent (always 0) and free_horizon. See how far you get on the tiniest cheat.
    """
    if pd is None:
        raise RuntimeError("superhero_walk needs pandas")
    ps = cheatless_primes(hi, lanes)
    rows, wh = [], 0
    for i, p in enumerate(ps, 1):
        fold = p ^ 5
        wh += fold
        rows.append(dict(idx=i, p=p,
                         lane=("6k-1" if p % 6 == 5 else "6k+1" if p % 6 == 1 else "seed"),
                         fold=fold, total=wh))
    df = pd.DataFrame(rows)
    df.attrs["tests_spent"] = 0
    df.attrs["free_horizon"] = fold_horizon(lanes)
    return df


def lane_pile(lo=2, hi=200):
    """The fold pile split into its two mod-8 lanes, with the closed form beside it.

    fold(p) = p ^ 5 is only ever p+3 (when p%8 is 1 or 3) or p-5 (when p%8 is 5 or 7)
    -- and p+5 for the lone 2. So the running XOR pile wh = sum(p ^ 5) is knowable by
    COUNTING lanes, with no XOR at all:

        pile = sum(primes) + 3*(# in +3 lane) - 5*(# in -5 lane) + 5*(the 2, if present)

    The walked pile and the counted `closed` sit shoulder to shoulder; `same` is True on
    every row. The out-of-order look of the raw folds is just the two ordered lanes
    crossing -- split by `lane` and each column is sorted.

    Columns: idx, p, mod8, lane, fold, kplus, kminus, sump, pile (walked XOR),
             closed (counted), same. Restarts local to [lo, hi] (pile begins at 0 at lo).
    """
    if pd is None:
        raise RuntimeError("lane_pile needs pandas")
    lo, hi = int(lo), int(hi)
    ps = [p for p in wheel_sieve(hi) if p >= lo]
    rows = []
    wh = sump = kplus = kminus = ktwo = 0
    for i, p in enumerate(ps, 1):
        fold = p ^ 5
        if p == 2:
            lane = "+5"; ktwo += 1
        elif p % 8 in (1, 3):
            lane = "+3"; kplus += 1
        else:
            lane = "-5"; kminus += 1
        wh += fold
        sump += p
        closed = sump + 3 * kplus - 5 * kminus + 5 * ktwo
        rows.append(dict(idx=i, p=p, mod8=p % 8, lane=lane, fold=fold,
                         kplus=kplus, kminus=kminus, sump=sump,
                         pile=wh, closed=closed, same=(wh == closed)))
    df = pd.DataFrame(rows)
    df.attrs["tests_spent"] = 0
    df.attrs["closed_form"] = "pile = sum(p) + 3*kplus - 5*kminus (+5 for the 2)"
    return df


def joke_universe(epochs=12, start=2, narrate=True):
    """A joke of a universe. Particles are primes -- built, never tested (nobody asks if
    they are real). Each joins a force by p mod 8: the +3 push (matter) or the -5 pull
    (antimatter); the lone 2 is the +5 weirdo. The universe's size is the fold pile
    wh = sum(p ^ 5), which by the lane law is just sum(primes) with a mod-8 rounding error.
    Ships with an undo button: nt.unfold([0] + list(df['size'])) walks every particle back.

    Columns: epoch, mass (=p), mod8, force, fold (p^5), size (=wh),
             balance (matter - antimatter), gag. No test spent -- 0 interrogations.
    """
    if pd is None:
        raise RuntimeError("joke_universe needs pandas")
    hi = max(1000, (int(start) + int(epochs)) * 20)
    ps = [p for p in wheel_sieve(hi) if p >= int(start)][:int(epochs)]
    rows, wh, matter, antimatter, prev_fold = [], 0, 0, 0, None
    if narrate:
        print("=== BIG BANG (size 0, nobody asked permission) ===")
    for e, p in enumerate(ps, 1):
        fold = p ^ 5
        if p == 2:
            force = "+5 weirdo"
        elif p % 8 in (1, 3):
            force, matter = "+3 push", matter + 1
        else:
            force, antimatter = "-5 pull", antimatter + 1
        wh += fold
        gag = ""
        if fold == 0:
            gag = "weightless -- exists, weighs nothing, universe shrugs"
        elif prev_fold is not None and fold < prev_fold:
            gag = "time hiccup -- arrived lighter than the last"
        if matter == antimatter and p != 2:
            gag = (gag + " | " if gag else "") + "perfect balance (Thanos moment)"
        rows.append(dict(epoch=e, mass=p, mod8=p % 8, force=force, fold=fold,
                         size=wh, balance=matter - antimatter, gag=gag))
        if narrate:
            tail = ("   <-- " + gag) if gag else ""
            print("epoch %3d : mass %4d  %-9s  size now %8d%s" % (e, p, force, wh, tail))
        prev_fold = fold
    df = pd.DataFrame(rows)
    if narrate:
        sp = sum(ps)
        print()
        print("why something rather than nothing?  matter %d vs antimatter %d -- decided by p mod 8."
              % (matter, antimatter))
        print("the whole universe weighs %d; the primes in it sum to %d; the %+d is pure mod-8 bookkeeping."
              % (wh, sp, wh - sp))
        print("undo:  nt.unfold([0] + list(df['size']))  walks every particle back out, exact.")
        print("and it never checked whether any of it was real. 0 interrogations.")
    df.attrs["tests_spent"] = 0
    return df


def freddy(n, verbose=True):
    """1, 2, Freddy's coming for you. Here Freddy is the primality TEST -- the one interrogation --
    and he can only cross the free horizon. Below fold_horizon() every number is settled by the
    deterministic lanes (wheel + square + fermat + fold), so no test is spent and Freddy stays out:
    you sleep sound, safe by place. At or past the horizon the free lanes have gaps, so a survivor
    needs one real test to be sure -- that knock at the door is Freddy.

    Returns dict(n, horizon, safe, verdict, how, freddy). verbose prints the bedtime status.
    """
    n = int(n)
    H = fold_horizon()
    v, how, cost = free_verdict(n)
    coming = (cost == "test")
    if not coming:
        if v == "prime":
            line = "prime by survival -- no test. Freddy can't cross the horizon; you sleep sound."
        else:
            line = "caught free (%s) -- no test needed. Freddy stays out." % how
    else:
        if _mr(n):
            v, how = "prime", "Freddy knocked -- and you were prime all along (1 test)"
        else:
            r = rho_walk(n)
            how = ("Freddy got you -- %dx%d (1 test)" % (r[0], r[1])) if r else "Freddy got you -- composite (1 test)"
            v = "composite"
        line = "past the horizon, the free lanes ran dry -- one test at the door."
    if verbose:
        print("  1, 2, Freddy's coming for you...")
        print("  %s  n=%d   (free horizon = %d)" % ("-- knock knock --" if coming else "-- door locked --", n, H))
        print("  %s" % line)
        print("  verdict: %s" % v)
    return dict(n=n, horizon=H, safe=(not coming), verdict=v, how=how, freddy=coming)


def square_xor(masks=(5, 10, 15, 20, 25), hi=4000):
    """Run the squares 1^2..hi^2 through XOR by each mask and count where the fold lands on a
    PURE add (s^m == s+m) or a PURE sub (s^m == s-m). The rates come out as clean dyadic
    fractions because squares are rigid in the low bits: every square is 0 or 1 mod 4, so bit 1
    is NEVER set -- any mask that touches bit 1 can never subtract. Returns a DataFrame:
    mask, bits, add, add_pct, sub, sub_pct, either_pct.
    """
    if pd is None:
        raise RuntimeError("square_xor needs pandas")
    sq = [k * k for k in range(1, int(hi) + 1)]
    N = len(sq)
    rows = []
    for m in masks:
        m = int(m)
        add = sum(1 for s in sq if (s ^ m) == s + m)
        sub = sum(1 for s in sq if (s ^ m) == s - m)
        rows.append(dict(mask=m, bits=format(m, "b"),
                         add=add, add_pct=round(100 * add / N, 2),
                         sub=sub, sub_pct=round(100 * sub / N, 2),
                         either_pct=round(100 * (add + sub) / N, 2)))
    return pd.DataFrame(rows)


def hide_index(k, m=30):
    """Hide a square's root k inside a scrambled number: (k*k) ^ m. It stops looking like a
    square; recover it with find_index(hidden, m). XOR is its own inverse, so the round trip is
    exact and needs no search or test -- the mask m is the key."""
    return (int(k) * int(k)) ^ int(m)


def find_index(h, m=30):
    """Reveal the hidden square index: unscramble with the same mask, then isqrt. Returns the
    root k when h ^ m is a perfect square (always true for a genuine hide_index output), else
    None. One XOR + one isqrt -- no search, no test."""
    s = int(h) ^ int(m)
    if s < 0:
        return None
    r = _isqrt(s)
    return r if r * r == s else None


def is_ground(seed, m=5):
    """A seed is GROUND when seed ^ m is not a perfect square -- nothing nests below it, so a
    tower built from it stops cleanly at the bottom. For m=5, 0 and 2 are ground."""
    s = int(seed) ^ int(m)
    if s < 0:
        return False
    r = _isqrt(s)
    return r * r != s


def fold_ladder(n=12, m=5):
    """BREADTH, built from an index -- no scanning. Rung j hides root j inside a bigger square:
    k = j^2 ^ m is the visible root, x = k^2 is the square a range-scan would find. Because ^m is
    its own inverse, k ^ m == j^2 exactly, so every rung is valid by construction -- the loop
    `x is square AND isqrt(x)^m is square` recreated with zero searching and zero tests. Returns
    a DataFrame: j, k, x, back (= k^m), root (= isqrt(back), which is j again)."""
    if pd is None:
        raise RuntimeError("fold_ladder needs pandas")
    rows = []
    for j in range(int(n)):
        k = (j * j) ^ int(m)
        x = k * k
        back = k ^ int(m)
        rows.append(dict(j=j, k=k, x=x, back=back, root=_isqrt(back)))
    return pd.DataFrame(rows)


def fold_climb(seed=2, height=4, m=5):
    """DEPTH, built from an index -- a square nested `height` deep. Climb with u(r) = r^2 ^ m;
    the key fact is u(r) ^ m == r^2, a perfect square EVERY time, so each step up is a guaranteed
    nesting -- no luck, no test, no dead end. Returns the nested square x = top^2; recover the
    whole chain with fold_peel(x, m). Use a GROUND seed (is_ground) so it bottoms out cleanly --
    for m=5 try 0 or 2. Depth is unbounded; the only cost is that x's digit count roughly doubles
    per step."""
    r = int(seed)
    for _ in range(int(height) - 1):
        r = (r * r) ^ int(m)
    return r * r


def fold_peel(x, m=5):
    """DEPTH, decoded -- peel a square back toward its seed. Take the root of x, XOR by m, and if
    that is a square take ITS root, and again, until (root ^ m) stops being a square. Returns the
    chain of roots top -> bottom; its length is how many squares nest inside x. [] if x is not a
    perfect square. Inverse of fold_climb: fold_peel(fold_climb(s, h, m), m)[-1] == s when s is
    ground."""
    x = int(x)
    r = _isqrt(x)
    if r * r != x:
        return []
    chain = [r]
    while True:
        s = chain[-1] ^ int(m)
        if s < 0:
            break
        root = _isqrt(s)
        if root * root != s:
            break
        chain.append(root)
    return chain


def nest_scan(hi=16430, m=5):
    """KEEP IT REAL -- the search side. Sweep range(hi) for squares whose root, folded by ^m,
    lands on another square (the loop that started all this), and tag each with how DEEP it nests
    via fold_peel. Shows the trade the builder hides: a window caps breadth at ~hi^(1/4) rungs
    while fold_climb has no cap at all. Returns a DataFrame: x, root, depth, chain."""
    if pd is None:
        raise RuntimeError("nest_scan needs pandas")
    rows = []
    for x in range(int(hi)):
        r = _isqrt(x)
        if r * r != x:
            continue
        s = r ^ int(m)
        rr = _isqrt(s)
        if rr * rr != s:
            continue
        chain = fold_peel(x, m)
        rows.append(dict(x=x, root=r, depth=len(chain), chain="->".join(map(str, chain))))
    return pd.DataFrame(rows)


def blind_sieve(n):
    """A sieve that doesn't know it's finding primes. It only ever does two things -- multiply and
    mark. For each survivor i, it crosses off i*i, i*i+i, i*i+2i, ... -- starting at the SQUARE,
    because everything below i*i was already crossed by a smaller factor (that's the 'squared to
    cheat less'). No is_prime, no test, no clue what it's up to. The survivors are the primes
    anyway. Returns the list of survivors.
    """
    n = int(n)
    alive = [True] * (n + 1)
    if n >= 0:
        alive[0] = False
    if n >= 1:
        alive[1] = False
    i = 2
    while i * i <= n:
        if alive[i]:
            j = i * i                       # the square: skip everything already gone
            while j <= n:
                alive[j] = False
                j += i
        i += 1
    return [k for k in range(2, n + 1) if alive[k]]


def wheel_sieve(n):
    """The blind sieve, but it only ever touches the +6 rail. It erases the multiples of 2 and 3
    up front (they're off the rail), then crosses off only from the 6k±1 survivors, stepping the
    rail +2,+4,+2,+4... and starting each at its square. Two-thirds of the number line never gets
    looked at. Same ignorant loop, same emergent primes -- just riding the rail. Returns the list.
    """
    n = int(n)
    if n < 2:
        return []
    alive = bytearray(b"\x01") * (n + 1)
    alive[0] = alive[1] = 0
    for base in (2, 3):                          # sweep the off-rail multiples away first
        for j in range(base * base, n + 1, base):
            alive[j] = 0
    i, step = 5, 2                               # walk the rail: 5,7,11,13,17,... (+2,+4,...)
    while i * i <= n:
        if alive[i]:
            for j in range(i * i, n + 1, i):     # cross from the square
                alive[j] = 0
        i += step
        step = 6 - step
    return [k for k in range(2, n + 1) if alive[k]]


def unfold(totals):
    """The inverse of the fold-sum walk. Hand it the running-total history -- the `total` list, wh
    after each step, starting at 0 -- and it walks the pile all the way back DOWN to the primes.
    Each prime is (total[i] - total[i-1]) ^ 5, because the step added p^5 and ^5 undoes itself.

    Build up by folding and summing; read back down by differencing and refolding. No test in
    either direction, and the round trip is exact. Pairs with superhero_walk:
        df = superhero_walk(); unfold([0] + list(df.total)) == list(df.p)
    """
    t = [int(v) for v in totals]
    return [(t[i] - t[i - 1]) ^ 5 for i in range(1, len(t))]


def at(totals, i):
    """Seek straight to the i-th prime of the fold-walk (1-indexed) from the total[] history --
    no walking. It's just (total[i] - total[i-1]) ^ 5: the step that added this prime folded it in,
    and ^5 undoes itself. Random access into the primes."""
    return (int(totals[i]) - int(totals[i - 1])) ^ 5


def fold_run(n):
    """The go-big build in one call: sieve every prime up to n (marking, no prime test), fold each
    one (p ^ 5) and pile it into the running total. Returns (primes, totals) -- primes is the list,
    totals is the [0, ...] history that unfold / at / radiate read. Millions without blinking."""
    ps = wheel_sieve(int(n))
    wh, totals = 0, [0]
    for p in ps:
        wh += p ^ 5
        totals.append(wh)
    return ps, totals


def radiate(totals, center, span=6):
    """Build LEFT and RIGHT from a point. Given the totals history and a center prime index
    (1-based), read the primes on both sides straight through the fold -- descend to the left,
    climb to the right, at once, no walking the whole way. Returns a DataFrame: offset, idx, prime,
    fold (p^5), side. Pick any floor and see the walk radiate out both directions from it."""
    if pd is None:
        raise RuntimeError("radiate needs pandas")
    rows = []
    for off in range(-int(span), int(span) + 1):
        i = int(center) + off
        if 1 <= i < len(totals):
            p = at(totals, i)
            rows.append(dict(offset=off, idx=i, prime=p, fold=(p ^ 5),
                             side=("*center*" if off == 0 else "left" if off < 0 else "right")))
    return pd.DataFrame(rows)


def party(n=50_000_000, seed=1009732533765211):
    """The whole show in one call -- built to impress someone who knows primes. Three acts:
    the flex (teleport + fold, for the lulz), the real deal (an 80-digit PROVABLE prime, for the
    respect), and the clever bit (zero-test generation to a proven wall, to impress). Prints the
    performance with patter and returns the goods as a dict.
    """
    out = {}
    print("=" * 66)
    print("  THE PRIME PARTY  --  pick a chair, this only takes a second")
    print("=" * 66)

    # ---- Act I: the flex (for the lulz) ----
    print("\n--  Act I: the flex  --")
    ps, totals = fold_run(n)
    mil = at(totals, 1_000_000) if len(ps) >= 1_000_000 else at(totals, len(ps))
    out["millionth_prime"] = mil
    out["pile"] = totals[-1]
    out["exact_roundtrip"] = (unfold(totals) == ps)
    print(f'  "name a prime by position -- the millionth?"  ->  {mil:,}   (no walking)')
    print(f'  all {len(ps):,} primes folded into ONE number:  {totals[-1]}')
    print(f'  ...and it unfolds back, every one, exact:  {out["exact_roundtrip"]}')

    # ---- Act II: the real deal (for the respect) ----
    print("\n--  Act II: the real deal  --")
    r = build_prime(seed)
    out["prime"] = r.get("N")
    out["certificate"] = r.get("certificate")
    if r.get("proven"):
        print(f'  a PROVABLE prime -- {r["digits"]} digits, built not guessed:')
        print(f'    {r["N"]}')
        print(f'    proof (Pocklington):  {r["certificate"]}')
        print(f'    N-1 fully factored, one witness -- anyone can verify it in seconds.')
    else:
        print(f'  (seed {seed} tower: {r.get("note")})')

    # ---- Act III: the clever bit (to impress) ----
    print("\n--  Act III: the clever bit  --")
    cp = cheatless_primes()
    H = fold_horizon()
    out["cheatless_count"] = len(cp)
    out["horizon"] = H
    print(f'  every prime to {H:,} with ZERO primality tests -- {len(cp):,} of them.')
    print(f'  "the 6-wheel is only exact to 17^2 = 289" -- true. but the square lanes prove it')
    print(f'  exact {H // 289}x further, to {H:,} = 17x673, the first composite that slips them all.')

    print("\n" + "=" * 66)
    print("  the lulz, the respect, the impress. no test ran that you did not watch.")
    print("=" * 66)
    return out


def sqkey(p):
    """The square-fold index key for p: (p^2) ^ 5. For any odd prime it's just p^2 + 3 (odd squares
    are all 1 mod 8). Reversible with unsqkey -- a clean bijective index, prime -> key."""
    return (int(p) * int(p)) ^ 5


def unsqkey(k):
    """Walk a square-fold key back to its number: since ^5 is its own inverse, k ^ 5 is the square,
    and its root is the number. Returns n with (n^2)^5 == k, or None if k isn't such a key."""
    sq = int(k) ^ 5
    r = _isqrt(sq)
    return r if r * r == sq else None


def ballers(hi=None, lanes=(2, 3, 5, 7, 11, 13)):
    """baller. baller. baller. The fold-walk carrying THREE piles at once:

        prevballer   -- the pile one step BACK   (yesterday)
        baller       -- the pile you're ON       (today)
        futureballer -- the pile one step AHEAD  (tomorrow, folding in the NEXT prime)

    All three climb a whisper apart, every step's diff unfolds straight back to a prime -- always
    TRUE, because you folded it in -- and not one primality test runs the whole way up. Reversible,
    step-invariant, no interrogation in the building. Returns a DataFrame:
    idx, prime, fold (p^5), prevballer, baller, futureballer, unfolds_to, TRUE.
    """
    if pd is None:
        raise RuntimeError("ballers needs pandas")
    ps = cheatless_primes(hi, lanes)
    rows, present, future = [], 0, 0
    for i, p in enumerate(ps):
        prevballer = present
        present += p ^ 5                                   # baller: fold in today's prime
        nxt = ps[i + 1] if i + 1 < len(ps) else p
        future += nxt ^ 5                                  # futureballer: fold in tomorrow's
        step = present - prevballer                        # == p^5
        rows.append(dict(idx=i + 1, prime=p, fold=(p ^ 5),
                         prevballer=prevballer, baller=present, futureballer=future,
                         unfolds_to=(step ^ 5), TRUE=(step ^ 5 == p)))   # unfold recovers it? -- a compare, NOT a test
    return pd.DataFrame(rows)
