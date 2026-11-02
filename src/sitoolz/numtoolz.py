# numtoolz - the disjoint half of the numtools family. Import as ntz.
#
#     import numtools as nt      # 162 defs, the main library
#     import numtoolz as ntz     #  50 defs, this one
#
# These two share ZERO function names. Verified by AST across every numtools
# variant on this machine: pylib's numtools has 50 top-level defs and not one
# of them appears in sitoolz/numtools.py. So nt and ntz complement rather than
# shadow each other - nothing to choose between, you can hold both.
#
# Provenance: copied verbatim from
#     Downloads/pylib-20260710T040032Z-2-001/pylib/numtools.py
# Depends only on numpy.
#
# The other variants on this machine are all strict SUBSETS of numtools.py and
# carry nothing unique - primebase/src (80 defs), Documents/larsdoku-3.1.9 (68),
# larsdoku-3.1.7 (58), larsdokuzonehtml (40). sitoolz/numtools_forlars.py is
# byte-identical to sitoolz/numtools.py.

"""numtools — small integer-set helpers.

Lars Rocha's number-theory playground, the half that depends on numpy alone.
Lives in ~/pylib and is auto-imported by IPython startup, so in any session you
can just call powers_and_sixes(...) directly. Also importable anywhere ~/pylib
is on sys.path:  `from numtools import powers_and_sixes`.
"""

import numpy as np

__all__ = ["powers_of_two", "multiples_of", "powers_and_sixes", "structured_index",
           "first_n_primes", "xor_prefix", "buildsystem", "recover_index",
           "xor_scan", "structured_center2", "unstructure", "restructure",
           "prime_col", "xor_outer", "primes_to_deltas", "deltas_to_primes",
           "cycle8", "cycle8_square", "paint",
           "MERSENNE_EXPONENTS", "mersenne_primes", "mersenne_via_cycle8",
           "perfect_numbers", "mersenne_sequence",
           "digital_root", "perfect_facts", "mersenne_facts", "mersenne_array",
           "even_sieve", "stripe", "xor_twos", "nearest_even", "zero_start",
           "nearest_pair", "sixes", "divchain", "divmatrix",
           "coprime_wheel", "dive_primes", "self_sieve", "RT", "rt_sieve",
           "rt_survivors", "extend_rt", "pure_prime_cap", "linear_sieve",
           "next_lane_prime", "prime_chain"]


def mersenne_array(count=20):
    """Return a numpy array of the first `count` Mersenne primes (2^p - 1).

    dtype=object so it holds the arbitrarily-large ones. There are 52 known;
    default 20 keeps it small (up to ~1300 digits). Pass count=52 for ALL --
    but the last is ~41 million digits, so only do that on purpose.

    >>> mersenne_array(5)
    array([3, 7, 31, 127, 8191], dtype=object)
    """
    return np.array([(1 << p) - 1 for p in MERSENNE_EXPONENTS[:count]],
                    dtype=object)


def digital_root(n):
    """Repeated digit-sum collapsed to 1..9 (0 for 0). Uses the mod-9 formula,
    so it's instant even for million-digit numbers.

    >>> digital_root(8128)
    1
    """
    n = abs(int(n))
    return 0 if n == 0 else 1 + (n - 1) % 9


def perfect_facts(count=6):
    """The party tricks of each perfect number 2^(p-1)*(2^p-1) -- yields a dict:

      * triangle_of : it equals 1+2+...+M  (the triangle of its Mersenne prime M)
      * odd_cubes_upto : (p>2) it equals 1^3 + 3^3 + ... + this^3
      * digital_root : always 1, except for 6

    >>> next(perfect_facts())['perfect']
    6
    """
    for p in MERSENNE_EXPONENTS[:count]:
        M = (1 << p) - 1
        perfect = M << (p - 1)
        d = {"p": p, "mersenne": M, "perfect": perfect,
             "triangle_of": M, "digital_root": digital_root(perfect)}
        if p > 2:
            d["odd_cubes_upto"] = (1 << ((p + 1) // 2)) - 1
        yield d


def mersenne_facts(count=12):
    """The party tricks of each Mersenne prime 2^p - 1 -- yields a dict:

      * ones : its binary is `p` ones, ALWAYS a palindrome
      * last_digit : only ever 1 or 7 (for p > 2), set by p mod 4
      * digital_root : only ever 1 or 4 (for p > 3), set by p mod 6

    So you know a Mersenne prime's last digit and digital root from its
    exponent alone -- no need to build the (possibly million-digit) number.

    >>> [d["last_digit"] for d in mersenne_facts(8) if d["p"] > 2]
    [7, 1, 7, 1, 1, 7, 7]
    """
    for p in MERSENNE_EXPONENTS[:count]:
        M = (1 << p) - 1
        yield {"p": p, "mersenne": M, "ones": p, "binary_palindrome": True,
               "last_digit": M % 10, "digital_root": digital_root(M)}


def mersenne_sequence(stop=127):
    """Reproduce the Mersenne sequence up to exponent `stop` (a size guard, so
    you never accidentally build a million-digit number and crash).

    For each Mersenne exponent p <= stop, yields a tuple:
        (shift, exponent, mersenne, perfect)
    where shift = p-1, mersenne = 2^p - 1, perfect = 2^(p-1)*(2^p - 1).

    stop is the MAX exponent. Safe sizes:
        stop<=127  -> up to ~39 digits     stop<=607  -> up to ~183 digits
        stop<=4423 -> up to ~1300 digits   stop>20000 -> thousands of digits+

    >>> [p for _, p, _, _ in mersenne_sequence(31)]
    [2, 3, 5, 7, 13, 17, 19, 31]
    """
    for p in MERSENNE_EXPONENTS:
        if p > stop:
            break
        m = (1 << p) - 1
        yield p - 1, p, m, m << (p - 1)


# The 52 known Mersenne-prime exponents p (i.e. 2**p - 1 is prime), in order.
MERSENNE_EXPONENTS = [
    2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281,
    3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243,
    110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377,
    6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657,
    37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933,
    136279841,
]


def mersenne_primes(count=None):
    """Yield (exponent, value) for the known Mersenne primes 2**p - 1.

    There are 52 known. Lazy on purpose -- the big ones have millions of
    digits, so you only pay for the ones you actually consume.

    >>> [v for _, v in mersenne_primes(5)]
    [3, 7, 31, 127, 8191]
    """
    exps = MERSENNE_EXPONENTS[:count] if count else MERSENNE_EXPONENTS
    for p in exps:
        yield p, (1 << p) - 1


def mersenne_via_cycle8(count=None):
    """Generate the known Mersenne primes the cycle8 way.

    cycle8 is 8*N (multiples of 8); at index N = 2^(p-3) - 1 the value is
    2^p - 8, and XOR 7 lifts it to 2^p - 1.  Same numbers as `mersenne_primes`,
    reached through the cycle8 corner identity.

    >>> [v for _, v in mersenne_via_cycle8(5)]
    [3, 7, 31, 127, 8191]
    """
    exps = MERSENNE_EXPONENTS[:count] if count else MERSENNE_EXPONENTS
    for p in exps:
        if p < 3:
            yield p, (1 << p) - 1
        else:
            yield p, (8 * ((1 << (p - 3)) - 1)) ^ 7


def perfect_numbers(count=None):
    """Yield (exponent, perfect_number) for the known even perfect numbers.

    Euclid-Euler: each comes from a Mersenne prime as 2^(p-1) * (2^p - 1).
    One per Mersenne prime, so 52 known.

    >>> [n for _, n in perfect_numbers(4)]
    [6, 28, 496, 8128]
    """
    exps = MERSENNE_EXPONENTS[:count] if count else MERSENNE_EXPONENTS
    for p in exps:
        yield p, (1 << (p - 1)) * ((1 << p) - 1)


def _make_stripe(kind, rows, mask=0xFFF):
    """Resolve a stripe name (or pass-through array) to a bounded int64 sequence."""
    if not isinstance(kind, str):
        return (np.asarray(kind, np.int64) & mask)
    i = np.arange(rows)
    if kind == "primes":
        seq = first_n_primes(rows)
    elif kind == "evens":
        seq = i * 2
    elif kind == "odds":
        seq = i * 2 + 1
    elif kind == "squares":
        seq = (i + 1) ** 2
    elif kind == "triangular":
        seq = (i + 1) * (i + 2) // 2
    elif kind == "fibonacci":
        f = [1, 2]
        while len(f) < rows:
            f.append(f[-1] + f[-2])
        return np.array([v & mask for v in f[:rows]], dtype=np.int64)
    else:
        raise ValueError(f"unknown stripe {kind!r}")
    return (np.asarray(seq) & mask).astype(np.int64)


def _prime_mask(arr):
    """Boolean array: which entries of `arr` are prime (sieve to max)."""
    arr = np.asarray(arr)
    hi = int(arr.max()) if arr.size else 1
    sv = np.ones(hi + 1, bool)
    sv[:2] = False
    for k in range(2, int(hi ** 0.5) + 1):
        if sv[k]:
            sv[k * k::k] = False
    return sv[arr]


def paint(stripe="primes", comb=None, tint=7, cmap="inferno", rows=384,
          mask=0xFFF, figsize=(16, 8), dpi=120, save=None, show=False, title=None,
          show_primes=False, prime_xor=3, prime_min=200, star_cmap="cool"):
    """Render a stripe x comb XOR canvas:  design[i,k] = (stripe[i] ^ comb[k]) ^ tint.

    stripe : 'primes'|'evens'|'odds'|'squares'|'triangular'|'fibonacci' or an array
    comb   : the column sequence (e.g. your ssx); defaults to structured_index(2048)
    tint   : XOR constant applied to the whole canvas
    cmap   : any matplotlib colormap name
    show_primes : if True, dim the canvas and light every cell whose
                  (value ^ prime_xor) is a prime > prime_min, as stars colored
                  by size with star_cmap.
    Returns the design array; pass save='foo.png' to write it.
    """
    import matplotlib.pyplot as plt
    s = _make_stripe(stripe, rows, mask)
    c = np.asarray(structured_index(2048) if comb is None else comb, np.int64)
    design = (s[:, None] ^ c[None, :]) ^ tint
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(design, cmap=cmap, interpolation="nearest", aspect="auto",
              alpha=0.55 if show_primes else 1.0)
    if show_primes:
        v = design ^ prime_xor
        hi = int(v.max())                       # inline sieve -> no external helper
        sv = np.ones(hi + 1, bool)
        sv[:2] = False
        for k in range(2, int(hi ** 0.5) + 1):
            if sv[k]:
                sv[k * k::k] = False
        star = sv[v] & (v > prime_min)
        ys, xs = np.where(star)
        ax.scatter(xs, ys, s=4, c=v[star], cmap=star_cmap, marker="*", linewidths=0)
    ax.set_xticks([]); ax.set_yticks([])
    if title:
        ax.set_title(title, color="0.9")
    fig.patch.set_facecolor("#08080f")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=dpi, facecolor="#08080f")
    if show:
        plt.show()
    return design


def cycle8(side, step=8):
    """Flat structured array of `side*side` multiples of `step` (default 8):
    0, step, 2*step, ...   Length is a perfect square, so it reshapes clean.

    >>> cycle8(32).shape           # 1024 = 32^2 multiples of 8
    (1024,)
    >>> int(cycle8(32)[-1])        # 1023 * 8
    8184
    """
    return np.arange(side * side, dtype=np.int64) * step


def cycle8_square(side, step=8):
    """`cycle8` reshaped to a (side, side) square of multiples of `step`.

    cycle8_square(32) -> 32x32 of 0,8,16,...,8184  (range(0, 8192, 8)).
    """
    return cycle8(side, step).reshape(side, side)


def primes_to_deltas(primes):
    """Forward half of the prime codec: the consecutive-XOR delta stream.

    delta[k] = primes[k] ^ primes[k+1].  Inverse of `deltas_to_primes`.

    >>> primes_to_deltas([2, 3, 5, 7, 11])
    array([ 1,  6,  2, 12])
    """
    p = np.asarray(primes, np.int64)
    return p[:-1] ^ p[1:]


def deltas_to_primes(a, seed=2):
    """Back half of the prime codec: replay the primes from a seed + deltas.

    Starts at `seed`, running-XORs through the delta stream `a`, and prepends
    the seed -- reconstructing the full prime sequence.

    >>> deltas_to_primes([1, 6, 2, 12]).tolist()
    [2, 3, 5, 7, 11]
    """
    return np.concatenate([[seed], xor_scan(a, seed)])


def prime_col(n):
    """The X first-column generator:  X[:n, 0] == first_n_primes(n) ^ 1.

    >>> prime_col(6)
    array([ 3,  2,  4,  6, 10, 12])
    """
    return first_n_primes(n) ^ 1


def xor_outer(a, b, sort='rows'):
    """Outer XOR table:  M[i, j] = a[i] ^ b[j].

    sort = 'rows' (sort within each row, default), 'cols', 'flat'
    (one fully-sorted 1D array), or None (raw).
    """
    M = np.asarray(a, np.int64)[:, None] ^ np.asarray(b, np.int64)[None, :]
    if sort == 'rows':
        return np.sort(M, axis=1)
    if sort == 'cols':
        return np.sort(M, axis=0)
    if sort == 'flat':
        return np.sort(M, axis=None)
    return M


def structured_center2(hi):
    """The structured lattice `sz`: structured_index(hi) with the center bumped
    1 -> 2 and the Mersenne cap (hi - 1) appended.  Length = len(idx) + 1.

    >>> structured_center2(128)[:6].tolist()
    [0, 2, 2, 4, 6, 8]
    """
    sz = structured_index(hi).copy()
    sz[1] = 2                                  # the '2 in the center'
    return np.append(sz, hi - 1)


def unstructure(primes, sz):
    """Build the unstructured list `rr` from the (unstructured) prime list and
    the structured lattice:  rr = primes ^ sz.  Pure XOR, fully invertible.
    """
    return np.asarray(primes, np.int64) ^ np.asarray(sz, np.int64)


def restructure(rr, sz):
    """Inverse of `unstructure`:  primes = rr ^ sz."""
    return np.asarray(rr, np.int64) ^ np.asarray(sz, np.int64)


def xor_scan(arr, seed=0):
    """Running (inclusive) XOR from a seed:

        out[k] = seed ^ arr[0] ^ arr[1] ^ ... ^ arr[k]

    Start at `seed`, XOR in each element, emit the running result. The
    exclusive cousin (output BEFORE consuming each element, leading 0) is
    `xor_prefix`.

    >>> xor_scan([2, 4, 6, 8, 10], seed=3)
    array([ 1,  5,  3, 11,  1])
    """
    return seed ^ np.bitwise_xor.accumulate(np.asarray(arr, dtype=np.int64))


def _is_prime(n):
    """True if int n is prime. Prefers gmpy2.is_prime (fast, big-int safe),
    falls back to a plain trial-division check."""
    n = int(n)
    if n < 2:
        return False
    try:
        import gmpy2
        return bool(gmpy2.is_prime(n))
    except ImportError:
        if n < 4:
            return True
        if n % 2 == 0:
            return False
        import math
        for k in range(3, math.isqrt(n) + 1, 2):
            if n % k == 0:
                return False
        return True


def divchain(n, width=None, dtype=object, tag_prime=False):
    """Halve-and-mark chain of `n`: interleave [value, flag, value, flag, ...].

    At each step you take the current value v and divide by 2 -- but if v is ODD
    you first do v ^ 1 (clear the low bit). Since (v ^ 1)//2 == v//2 for odd v,
    the next value is always v >> 1, and the flag BETWEEN the values records
    whether the ^1 was needed:  0 = even (no ^1),  1 = odd (^1 applied).  Those
    flags are exactly the binary digits of n, least-significant first.

    When v hits 0, every remaining slot is filled with the DIVISION COUNT (how
    many halvings it took to get there). `dtype=object` by default so it holds
    arbitrarily large n (2**51+1 and friends); pass dtype=np.int64 for speed
    when n stays under 2**63. width defaults to 2*(bit_length+1), leaving a
    small zero-filled tail; pass an explicit width to align many rows (see
    `divmatrix`).

    tag_prime=True prepends a leading primality tag to the row: 1 if n is prime,
    else 0 -- so prime and composite chains are easy to tell apart at a glance.
    With the tag on, the value/flag pairs start at index 1 instead of 0.

    >>> divchain(20)
    array([20, 0, 10, 0, 5, 1, 2, 0, 1, 1, 5, 5], dtype=object)
    >>> divchain(20)[1::2][:5].tolist()      # the flags == bits of 20, LSB-first
    [0, 0, 1, 0, 1]
    >>> divchain(19, tag_prime=True)[0]      # 19 is prime -> leading 1
    1
    >>> divchain(20, tag_prime=True)[0]      # 20 is composite -> leading 0
    0
    """
    n = int(n)
    steps = n.bit_length()                     # halvings until 0 (0 when n == 0)
    if width is None:
        width = 2 * (steps + 1)                # room for the zero-fill tail
    lead = 1 if tag_prime else 0
    out = np.empty(width + lead, dtype=dtype)
    if tag_prime:
        out[0] = int(_is_prime(n))             # 1 = prime, 0 = not
    v, count, i = n, 0, lead
    while v > 0 and i + 1 < width + lead:
        out[i] = v                             # the value
        out[i + 1] = v & 1                     # 1 if a ^1 was needed (odd) else 0
        v >>= 1                                # (v ^ (v & 1)) // 2  ==  v // 2
        count += 1
        i += 2
    out[i:] = count                            # hit zero: remaining = # divisions
    return out


def divmatrix(numbers, width=None, dtype=object, tag_prime=False):
    """Stack `divchain` over many numbers into one 2D array (one row each).

    Every row shares the same width so they line up -- handy for feeding a whole
    sweep (like your r + range(...) loop) straight into a numpy grid. If width
    is None it's sized to the largest number's chain. tag_prime=True adds the
    leading prime tag to every row (column 0 is then 1/0 = prime/not).

    >>> divmatrix([4, 20]).shape[0]
    2
    >>> divmatrix([19, 20], tag_prime=True)[:, 0].tolist()   # 19 prime, 20 not
    [1, 0]
    """
    nums = [int(x) for x in numbers]
    if width is None:
        steps = max((x.bit_length() for x in nums), default=0)
        width = 2 * (steps + 1)
    return np.array([divchain(x, width, dtype, tag_prime) for x in nums],
                    dtype=dtype)


def coprime_wheel(lo, hi, m=6):
    """The wheel candidates in [lo, hi): every n with gcd(n, m) == 1, PLUS the
    prime factors of m themselves if they land in range.

    This is the "dive" -- the only slots a prime can occupy once you've thrown
    away the multiples of m's prime factors. m=6 keeps 6k+/-1 (drops all evens
    and multiples of 3, 2/3 of the line gone). m=30 keeps the 8 residues
    {1,7,11,13,17,19,23,29} -- your cycle8 -- and drops 2,3,5's multiples
    (73% gone). dtype=object so it's safe out near 2**51.

    It does NOT test primality -- composite candidates (25, 35, 49, ...) still
    ride along. Feed the result to `dive_primes` (or gmpy2.is_prime) to finish.

    >>> coprime_wheel(2, 30).tolist()
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 25, 29]
    >>> coprime_wheel(2, 30, m=30).tolist()
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    """
    import math
    lo, hi = int(lo), int(hi)
    res = [r for r in range(m) if math.gcd(r, m) == 1]     # coprime residues
    factors = [p for p in range(2, m + 1) if m % p == 0    # 2,3 (,5) themselves
               and all(p % q for q in range(2, p))]
    keep = sorted(set(factors) | set())
    out = [p for p in keep if lo <= p < hi]
    start = (lo // m) * m
    for base in range(start, hi + m, m):
        for r in res:
            n = base + r
            if lo <= n < hi and n > 1:
                out.append(n)
    return np.array(sorted(set(out)), dtype=object)


def dive_primes(lo, hi, m=6):
    """Primes in [lo, hi), reached by diving with the mod-`m` wheel and then
    confirming each survivor with a primality test (gmpy2.is_prime, falling back
    to trial division). Big-int safe -- builds no size-`hi` array -- so it works
    on windows way out near 2**51 the way your r + range(...) sweep does.

    m is the wheel modulus: 6 (6k+/-1) or 30 (the cycle8 residues) dive hardest
    for the least overhead.

    >>> dive_primes(2, 30).tolist()
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    >>> dive_primes(100, 130).tolist()
    [101, 103, 107, 109, 113, 127]
    """
    cand = coprime_wheel(lo, hi, m)
    return np.array([int(n) for n in cand if _is_prime(int(n))], dtype=object)


def self_sieve(nums, cap=8192, inclusive=False):
    """Give every number its OWN cycle: for each n in `nums`, the multiples
    n, 2n, 3n, ...  up to `cap`. This is exactly the slice np.arange(cap)[n::n]
    -- n striding by itself, keeping n as the first term.

    Returns a dict {n: array-of-its-multiples}. `cap` is the ceiling you set
    (default 8192); by default it is EXCLUSIVE, matching the arange(cap)[n::n]
    slice (so the 4-cycle ends at 8188). Pass inclusive=True to let a multiple
    that lands exactly on `cap` be kept (4-cycle then ends at 8192).

    >>> s = self_sieve([4, 6, 10], cap=32)
    >>> s[4].tolist()
    [4, 8, 12, 16, 20, 24, 28]
    >>> s[6].tolist()
    [6, 12, 18, 24, 30]
    >>> int(self_sieve([4], cap=8192)[4][-1])     # exclusive: stops before 8192
    8188
    >>> int(self_sieve([4], cap=8192, inclusive=True)[4][-1])
    8192
    """
    top = int(cap) + 1 if inclusive else int(cap)
    return {int(n): np.arange(int(n), top, int(n), dtype=np.int64)
            for n in nums}


# Your saved starting set `rt` -- the even_sieve(256) survivors, frozen here so
# rt_sieve can lean on it without you passing anything.
RT = np.array([4, 6, 10, 14, 22, 26, 34, 38, 46, 58, 62, 74, 82, 86, 94, 106,
               118, 122, 134, 142, 146, 158, 166, 178, 194, 202, 206, 214, 218,
               226, 254], dtype=np.int64)


def rt_sieve(cap=8192, inclusive=False, nums=None):
    """Self-sieve the saved `RT` set: each of its numbers cycles by ITSELF
    (n, 2n, 3n, ...) up to `cap`, and the original number is kept as the first
    term. This is `self_sieve` with your array baked in, so you just call it.

    Pass `nums` to sieve a different set instead of RT; `cap` (default 8192) and
    `inclusive` behave exactly as in `self_sieve`.

    >>> s = rt_sieve(cap=32)
    >>> sorted(s)[:5]                     # the original numbers are the keys
    [4, 6, 10, 14, 22]
    >>> s[4].tolist()
    [4, 8, 12, 16, 20, 24, 28]
    >>> s[6].tolist()
    [6, 12, 18, 24, 30]
    """
    return self_sieve(RT if nums is None else nums, cap, inclusive)


def rt_survivors(cap=8192, inclusive=False, nums=None):
    """Overlay every RT cycle onto one number line and return the EVENS that no
    cycle ever hits -- the evens that are not a multiple of any number in RT.

    Each n in RT strikes n, 2n, 3n, ... (its `rt_sieve` cycle). Union all those
    strikes; whatever even number is left standing (never hit) is a survivor.
    Since every RT number strikes itself, none of them survive -- the survivors
    are a disjoint set of evens that the whole array steps over.

    Pass `nums` to sieve by a different set; `cap`/`inclusive` as in rt_sieve.
    Only 2 survives small caps (the 4-cycle kills every other power of two, the
    small primes cover the rest); real survivors start at 2*131 == 262, the
    first even whose half dodges all primes <= 127.

    >>> rt_survivors(cap=64).tolist()
    [2]
    >>> rt_survivors(cap=280).tolist()
    [2, 262, 274, 278]
    """
    nums = RT if nums is None else nums
    top = int(cap) + 1 if inclusive else int(cap)
    hit = np.zeros(top, dtype=bool)
    for n in nums:
        n = int(n)
        hit[n::n] = True                       # strike n, 2n, 3n, ... < top
    evens = np.arange(2, top, 2, dtype=np.int64)
    return evens[~hit[evens]]                   # evens never struck


def extend_rt(k):
    """A bigger RT: 2 * (the first k primes). RT itself is extend_rt(31)
    (2 * primes up to 127). Use this to dial the base set larger -- the more
    primes you fold in, the further the survivors stay pure prime
    (see `pure_prime_cap`).

    >>> extend_rt(5).tolist()
    [4, 6, 10, 14, 22]
    >>> (extend_rt(31) == RT).all()
    np.True_
    """
    return 2 * first_n_primes(int(k))


def pure_prime_cap(base=None):
    """The guaranteed composite-free ceiling for a `base` set (default RT).

    Sieving evens by base = 2*{primes up to P} leaves survivors 2*m where m has
    no prime factor <= P. The SMALLEST composite m that dodges every prime <= P
    is (nextprime(P))**2, so every survivor with value below 2*nextprime(P)**2
    is guaranteed 2*prime. This returns that ceiling (exclusive): run
    rt_survivors up to it and every half is prime.

    For the default RT (P=127, nextprime 131): 2*131**2 == 34322, matching the
    first composite survivor we found.

    >>> pure_prime_cap()
    34322
    >>> pure_prime_cap(extend_rt(100))     # primes up to 541 -> nextprime 547
    598418
    """
    base = RT if base is None else np.asarray(base)
    P = int(base.max()) // 2
    try:
        import gmpy2
        q = int(gmpy2.next_prime(P))
    except ImportError:
        q = P + 1
        while not _is_prime(q):
            q += 1
    return 2 * q * q


def linear_sieve(n, return_spf=False):
    """Euler's linear sieve: primes below `n`, striking each composite EXACTLY
    ONCE -- via its smallest prime factor -- so there is zero overlap (the
    redundancy you spotted in the RT cycles). Runs in O(n).

    The trick: for each i, mark i*p for primes p in increasing order, but STOP
    the moment p divides i. That guarantees every composite is reached only from
    its smallest prime factor, never twice.

    return_spf=True also returns the smallest-prime-factor array (spf[x] is the
    least prime dividing x), which is the proof there's no overlap: it's written
    once per composite.

    >>> linear_sieve(20).tolist()
    [2, 3, 5, 7, 11, 13, 17, 19]
    >>> primes, spf = linear_sieve(13, return_spf=True)
    >>> spf.tolist()             # spf[12]==2, spf[9]==3, ...
    [0, 0, 2, 3, 2, 5, 2, 7, 2, 3, 2, 11, 2]
    """
    n = int(n)
    spf = np.zeros(n, dtype=np.int64)          # 0 = untouched (=> prime)
    primes = []
    for i in range(2, n):
        if spf[i] == 0:                        # i survived => prime
            spf[i] = i
            primes.append(i)
        for p in primes:
            ip = i * p
            if ip >= n:
                break
            spf[ip] = p                        # struck once, by smallest factor p
            if i % p == 0:                     # p is i's smallest factor -> stop
                break
    parr = np.array(primes, dtype=np.int64)
    return (parr, spf) if return_spf else parr


def next_lane_prime(n, kmax=1_000_000):
    """Find a prime of the form (n+k)*2 + 1 by walking the 2*3 wheel lane.

    Candidate v = 2*n + 2*k + 1 is always odd, and is divisible by 3 whenever
    k lands v on a multiple of 3 -- those are skipped (a third of all offsets),
    so no doomed test is ever run. Returns (k, v) for the first prime, or
    (None, None) if none within kmax steps.

    This is plain candidate-generate-and-test (gmpy2.is_prime), just with the
    multiples of 2 and 3 structurally skipped -- the same wheel as coprime_wheel,
    driven from the doubling side. Big-int safe.

    >>> k, v = next_lane_prime(20)          # 2*20+1=41 already prime
    >>> k, v
    (0, 41)
    """
    n = int(n)
    for k in range(int(kmax)):
        v = 2 * n + 2 * k + 1
        if v % 3 == 0:                      # off the 2*3 lane -> can't be prime
            continue
        if _is_prime(v):
            return k, v
    return None, None


def prime_chain(seed, length=10, kmax=1_000_000):
    """Generate a chain of primes: from `seed`, find the next lane prime p
    (see next_lane_prime), record it, then feed 2*p back as the new seed and
    repeat `length` times. Returns a list of (k_offset, prime).

    Each prime is confirmed by a primality test, so the list is all-primes by
    construction. Handy for minting verified big primes on demand -- but note
    it is generate-and-test, not a closed-form prime formula (none exist);
    it's reliable because primes stay dense along the 2*3 lane.

    >>> [p for _, p in prime_chain(20, length=3)]
    [41, 167, 673]
    """
    seed = int(seed)
    out = []
    for _ in range(int(length)):
        k, p = next_lane_prime(seed, kmax)
        if p is None:
            break
        out.append((k, p))
        seed = 2 * p
    return out


def xor_prefix(index):
    """Exclusive prefix-XOR scan:  index -> the anchor row `a` (leading 0).

    a[k] = XOR of index[:k].  This is the row-0 / anchor sequence that the whole
    build system hangs off of. Inverse of `recover_index`.

    >>> xor_prefix([0, 1, 2, 4, 6])
    array([0, 0, 1, 3, 7, 1])
    """
    index = np.asarray(index, dtype=np.int64)
    return np.concatenate(([0], np.bitwise_xor.accumulate(index)))


def buildsystem(index):
    """Symmetric XOR Cayley table from an index:  T[i, j] = a_i ^ a_j,
    where a = xor_prefix(index).

    The table is mirror-symmetric with a zero diagonal, and the ENTIRE thing
    rebuilds from its own row 0 (which is `a`). So storage collapses to the
    index alone (itself a formula via structured_index).
    """
    a = xor_prefix(index)
    return a[:, None] ^ a[None, :]


def recover_index(a):
    """Inverse of `xor_prefix`:  recover the generating index from the anchor
    row `a`, or from a full table (uses its row 0).  idx[k] = a[k] ^ a[k+1].

    Round-trips exactly:  recover_index(buildsystem(idx)) == idx.
    """
    a = np.asarray(a, dtype=np.int64)
    if a.ndim == 2:                 # a full table -> its first row is the anchor
        a = a[0]
    return a[:-1] ^ a[1:]


def first_n_primes(n):
    """The first `n` primes as an int64 array (count-based, not value-bound).

    Sieves to an upper bound for the n-th prime (n*(ln n + ln ln n)) and slices,
    which is fast in bulk. Falls back to gmpy2.next_prime for tiny n.

    >>> int(first_n_primes(512)[-1])
    3671
    """
    n = int(n)
    if n <= 0:
        return np.empty(0, dtype=np.int64)
    if n < 6:
        return np.array([2, 3, 5, 7, 11][:n], dtype=np.int64)
    hi = int(n * (np.log(n) + np.log(np.log(n)))) + 10      # > n-th prime (Rosser)
    s = np.ones(hi, dtype=bool)
    s[:2] = False
    for k in range(2, int(hi ** 0.5) + 1):
        if s[k]:
            s[k * k::k] = False
    return np.flatnonzero(s)[:n].astype(np.int64)


def sixes(var, full=False):
    """The sixes-lattice anchor for 2*var:  6 * (var // 3)  -- the multiple of 6
    at or just below 2*var.  This is the closed form you derived; it satisfies

        2*var == sixes(var) + 2*(var % 3)        (the division algorithm)

    so the remainder 2*(var%3) (0, 2, or 4) says how far 2*var sits above the
    anchor: +0 lands on the sixes lattice, +2 on the 2-lattice, +4 between them.

    Default returns just the anchor. full=True returns
    (anchor, offset, target) where target == 2*var == anchor + offset.

    >>> sixes(63)
    126
    >>> sixes(61)
    120
    >>> sixes(65, full=True)
    (126, 4, 130)
    """
    anchor = 6 * (var // 3)
    if full:
        return anchor, 2 * (var % 3), 2 * var
    return anchor


def nearest_pair(var, top=1024):
    """Nearest value to 2*var across BOTH even lattices -- arange(6, top, 6)
    (residue 0 mod 6) and arange(2, top, 6) (residue 2 mod 6) -- by arithmetic,
    no array. Returns (nearest_value, nearest_value ^ 2*var).

    With t = 2*var and r = t % 6:
      r == 0 -> t itself (in the ...6,12,18 lattice)   -> xor 0
      r == 2 -> t itself (in the ...2,8,14 lattice)    -> xor 0
      r == 4 -> t is in neither; nearest is t-2 (t+2 ties) -> off by 2

    For ODD var the xor is therefore only ever 0 (var % 3 in {0,1}) or 2
    (var % 3 == 2).

    >>> nearest_pair(63)     # 126, on the 6-lattice
    (126, 0)
    >>> nearest_pair(61)     # 122, on the 2-lattice
    (122, 0)
    >>> nearest_pair(65)     # 130 -> off by 2
    (128, 2)
    """
    t = 2 * var
    val = t if t % 6 in (0, 2) else t - 2     # r == 4 misses by 2 (t+2 ties)
    val = min(max(val, 2), 2 * ((top - 1) // 2))   # clamp into [2, top)
    return int(val), int(val ^ t)


def zero_start(var, starts=(6, 2, 4), top=None):
    """Which start makes arange(start, top, 3) land EXACTLY on 2*var (so the
    XOR is 0) -- computed from var alone, no array, no search.

    The even members of arange(start, top, 3) all share residue `start % 6`, and
    the target 2*var is even, so the exact hit is the start with
    start % 6 == (2*var) % 6.  Equivalently by var % 3:  0->6, 1->2, 2->4.

    Returns the matching start, or None if 2*var exceeds `top` (when given) so no
    in-range hit exists.

    >>> zero_start(63)        # 126 -> %6 == 0
    6
    >>> zero_start(61)        # 122 -> %6 == 2
    2
    >>> zero_start(62)        # 124 -> %6 == 4
    4
    """
    if top is not None and 2 * var >= top:
        return None
    target_mod = (2 * var) % 6
    for s in starts:
        if s % 6 == target_mod:
            return s
    return None


def nearest_even(top, var, start=6, step=3):
    """The even number in np.arange(start, top, step) closest to the target
    var*2, found by arithmetic -- NO array is built.

    The even members of an odd-step arange form their own progression (first
    even term e0, spacing 2*step), so the nearest one to the target is a single
    round(): n = round((target - e0) / spacing), clamped into range.

    Returns (value, value ^ target). For the original 256/61 case this lands on
    the same numbers xor_twos finds the slow way:

    >>> nearest_even(256, 61, start=6)
    (120, 2)
    >>> nearest_even(256, 61, start=2)
    (122, 0)
    """
    t = var * 2
    e0 = start if start % 2 == 0 else start + step    # first even term
    spacing = step if step % 2 == 0 else 2 * step     # gap between evens
    last = e0 + spacing * ((top - 1 - e0) // spacing)  # largest even in range
    val = min(max(e0 + spacing * round((t - e0) / spacing), e0), last)
    return int(val), int(val ^ t)


def xor_twos(top=256, var=61, step=3, starts=(6, 2), evens_only=True):
    """Build  (np.arange(start, top, step) ^ (var*2))  for each start in `starts`
    and count the 2's several ways. `top` and `var` are the knobs you pass in
    (256 and 61 in the original); var*2 is the XOR constant (122).

    If evens_only, each arange is filtered to just its EVEN members before the
    XOR (the even numbers that live in that np range).

    Returns {start: dict} with, per line:
      base     : the (optionally even-filtered) arange
      line     : base ^ (var*2)
      value2   : how many entries equal the literal value 2
      residue2 : how many entries have (line % 3) == 2   (the [2,0,1] pattern)
      n        : length of the line

    >>> r = xor_twos(256, 61)
    >>> r[6]['n'], r[6]['residue2']
    (42, 21)
    >>> r[2]['n'], r[2]['residue2']
    (43, 10)
    """
    k = var * 2
    out = {}
    for s in starts:
        base = np.arange(s, top, step, dtype=np.int64)
        if evens_only:
            base = base[base % 2 == 0]
        line = base ^ k
        out[s] = {"base": base, "line": line,
                  "value2": int((line == 2).sum()),
                  "residue2": int((line % 3 == 2).sum()),
                  "n": int(line.size)}
    return out


def stripe(top, start=1, step=4, xor=3, mod=3, drop=0, keep_matched=False,
           transform=False):
    """Stripe an arange by an (x ^ xor) % mod residue pattern.

    Builds  a = np.arange(start, top, step)  and  residue = (a ^ xor) % mod.
    The residue cycles every `mod` elements -- i.e. every step*mod in value --
    so the matched positions are themselves a clean arange.

    By default returns `a` with every element whose residue == `drop` STRIPPED
    OUT (the zeros, for the default pattern). Options:
      keep_matched=True -> return just the matched (dropped) elements instead
      transform=True    -> return the XOR'd values (a ^ xor), i.e. the actual
                           sequence you see, not the raw arange

    For your case `^3` is just +1, so transform=True turns the raw 1,9,13,...
    into the cool 2,10,14,...  To begin the output at 26, start one below it
    (start=25, since ^3 adds 1):  stripe(227, start=25, transform=True).

    >>> stripe(227)[:6].tolist()
    [1, 9, 13, 21, 25, 33]
    >>> stripe(227, keep_matched=True).tolist()
    [5, 17, 29, 41, 53, 65, 77, 89, 101, 113, 125, 137, 149, 161, 173, 185, 197, 209, 221]
    >>> stripe(227, start=25, transform=True)[:6].tolist()
    [26, 34, 38, 46, 50, 58]
    """
    a = np.arange(start, top, step, dtype=np.int64)
    matched = (a ^ xor) % mod == drop
    sel = a[matched] if keep_matched else a[~matched]
    return (sel ^ xor) if transform else sel


def _is_square(n):
    """True if non-negative int n is a perfect square. Prefers gmpy2.is_square,
    falls back to math.isqrt -- both exact (no float rounding)."""
    if n < 0:
        return False
    try:
        import gmpy2
        return bool(gmpy2.is_square(int(n)))
    except ImportError:
        import math
        r = math.isqrt(int(n))
        return r * r == n


def even_sieve(top=256, steps=(4, 6, 10, 22, 14), strip=(8, 16, 388),
               drop_two=True, drop_half_squares=True):
    """Even-number sieve over [0, top): start from every even number, then strip
    out a stack of arithmetic progressions -- KEEPING the first term of each.

    The build, in order:
      * all evens 0, 2, 4, ... top-2
      * for each n in `strip`, remove every multiple of n (incl. 0) -- none kept
        (8, 16, 388 by default)
      * for each n in `steps`, remove multiples of n but keep n itself
        (np.arange(n, top, n) minus the leading n)
      * if drop_two, also remove 2
      * if drop_half_squares, remove any n whose half is a perfect square,
        i.e. is_square(n // 2)   (e.g. 338 -> 338//2 == 169 == 13**2)

    `top` is the only thing you usually change -- it's the range you're building
    over (256 here). Pass your own `steps`/`strip` to sieve by a different set,
    or flip drop_two / drop_half_squares off.

    >>> even_sieve(256)[:6]
    array([ 4,  6, 10, 14, 22, 26])
    >>> len(even_sieve(256))
    31
    >>> len(even_sieve(256, drop_two=False, drop_half_squares=False))
    32
    """
    evens = np.arange(0, top, 2, dtype=np.int64)
    excluded = [np.arange(0, top, n, dtype=np.int64) for n in strip]  # none kept
    for n in steps:
        m = np.arange(n, top, n, dtype=np.int64)
        excluded.append(m[m != n])                       # keep the first term (n)
    if drop_two:
        excluded.append(np.array([2], dtype=np.int64))
    remove = np.unique(np.concatenate(excluded))
    out = np.setdiff1d(evens, remove)
    if drop_half_squares:                                # drop n where n//2 is square
        out = out[[not _is_square(int(n) // 2) for n in out]]
    return out


def powers_of_two(lo, hi):
    """All powers of two p with lo <= p <= hi, as a sorted int array.

    >>> powers_of_two(2, 128)
    array([  2,   4,   8,  16,  32,  64, 128])
    """
    return np.array([1 << k for k in range(int(hi).bit_length() + 1)
                     if lo <= (1 << k) <= hi], dtype=np.int64)


def multiples_of(step, lo, hi):
    """Multiples of `step` in the half-open range [lo, hi).

    >>> multiples_of(6, 0, 128)[:5]
    array([ 0,  6, 12, 18, 24])
    """
    first = lo + (-lo) % step                      # smallest multiple of step >= lo
    return np.arange(first, hi, step, dtype=np.int64)


def powers_and_sixes(hi, start6=0, lo=2, step=6, keep_dupes=False):
    """Powers of two in [lo, hi] merged with multiples of `step` in [start6, hi).

    Returns ONE array. By default it is sorted and de-duplicated (union); pass
    keep_dupes=True to keep values that appear in both lists.

    >>> powers_and_sixes(128)
    array([  0,   2,   4,   6,   8,  12,  16,  18,  24,  30,  32,  36,  42,
            48,  54,  60,  64,  66,  72,  78,  84,  90,  96, 102, 108, 114,
           120, 126, 128])
    """
    powers = powers_of_two(lo, hi)
    sixes = multiples_of(step, start6, hi)
    if keep_dupes:
        return np.sort(np.concatenate([powers, sixes]))
    return np.union1d(powers, sixes)


def structured_index(hi, step=6, lo=1, end_offset=4):
    """Column-index set for slicing the gasket table.

    Powers of two and multiples of `step`, both STRICTLY BELOW `hi` (the top
    power-of-two boundary is excluded), plus the end-of-range marker
    `hi - end_offset` that sits just inside that boundary. Sorted, de-duped.

    For hi=128 this reproduces the hand-built `A` (a square 30-wide index that
    pairs with X[1:31] to cut a 30x30 block):

    >>> structured_index(128)
    array([  0,   1,   2,   4,   6,   8,  12,  16,  18,  24,  30,  32,  36,
            42,  48,  54,  60,  64,  66,  72,  78,  84,  90,  96, 102, 108,
           114, 120, 124, 126])
    """
    powers = powers_of_two(lo, hi - 1)          # 2^0..2^(k-1); excludes hi itself
    steps = multiples_of(step, 0, hi)           # 0..hi-step
    end = np.array([hi - end_offset], dtype=np.int64)
    return np.union1d(np.union1d(powers, steps), end)
