"""anchors -- the two ends of an octave, and where the property breaks.

Every octave [2^k, 2^(k+1)) is straddled by two anchors, and both read the
primes inside it by pure arithmetic -- because of XOR:

    FLOOR   2^k + 1        Fermat-form, TWO bits    XOR by even w < 2^k IS addition
    CEIL    2^(k+1) - 1    Mersenne-form, all ones  XOR IS subtraction

They are mirrors.  For every prime p in the octave:

    floor_wrap + ceil_wrap  ==  2^k - 2

The reason both work is the same one: `a ^ b == a + b` exactly when `a & b == 0`.
The floor is two bits so an even wrap below 2^k never collides; the ceiling is
all ones so XOR is subtraction outright.

**Where it breaks.**  The property that BOTH anchors are themselves prime holds
at k = 1, 2, 4, 16 and never again.  k = 16 -- floor 65537 = F_4, ceiling
131071 = M_17 -- is the last one.  Two conditions must coincide:

    floor prime  <=>  k is a power of two          (a^m+1 has factor a+1 for odd m)
    ceil prime   <=>  k+1 is a Mersenne exponent

k = 8 fails the second (511 = 7*73).  Past k = 16 the floor fails every
remaining chance, because the Fermat primes stop -- which is open.

    from sitoolz import anchors

    anchors.octave(8)             # both rulers over one octave, as a frame
    anchors.invariant_check(12)   # the mirror identity, every octave to k=12
    anchors.anchor_frame(40)      # is each anchor prime?  where both are
    anchors.both_prime(64)        # the k where both anchors are prime
    anchors.after_16(1024)        # what the floor does past the break
    anchors.report()              # the whole picture in one call

Companion to numtools.ceiling_floor / gap_ruler / under_p2, which is where this
came from.  Stdlib + gmpy2 if present.
"""
try:
    import gmpy2 as _g
except ImportError:
    _g = None


def _is_prime(n):
    n = int(n)
    if _g is not None:
        return bool(_g.is_prime(n))
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


def _frame(rows, index=None):
    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        return df.set_index(index) if index else df
    except ImportError:
        from sitoolz.dfkit import frame
        return frame(rows, index=index)


def floor_anchor(k):
    """2^k + 1 -- the Fermat-form end. Two bits."""
    return (1 << k) + 1


def ceil_anchor(k):
    """2^(k+1) - 1 -- the Mersenne-form end. All ones."""
    return (1 << (k + 1)) - 1


def octave(k):
    """Every prime in [2^k, 2^(k+1)), measured from both ends.

      prime        the prime
      floor_wrap   p - (2^k + 1), reachable by XOR since the floor is two bits
      ceil_wrap    (2^(k+1) - 1) ^ p, which is plain subtraction
      total        the two summed -- constant, = 2^k - 2
      xor_is_add   does (floor ^ floor_wrap) == p, i.e. no bit collision
    """
    fl, ce, inv = floor_anchor(k), ceil_anchor(k), (1 << k) - 2
    rows = []
    for p in range(fl, 1 << (k + 1), 2):
        if not _is_prime(p):
            continue
        fw, cw = p - fl, ce ^ p
        rows.append({"prime": p, "floor_wrap": fw, "ceil_wrap": cw,
                     "total": fw + cw, "invariant": inv,
                     "xor_is_add": (fl ^ fw) == p,
                     "holds": fw + cw == inv})
    return _frame(rows, index="prime")


def invariant_check(kmax=12):
    """The mirror identity floor_wrap + ceil_wrap == 2^k - 2, octave by octave."""
    rows = []
    for k in range(2, kmax + 1):
        fl, ce, inv = floor_anchor(k), ceil_anchor(k), (1 << k) - 2
        ps = [p for p in range(fl, 1 << (k + 1), 2) if _is_prime(p)]
        bad = [p for p in ps if (p - fl) + (ce ^ p) != inv]
        xr = [p for p in ps if (fl ^ (p - fl)) != p]
        rows.append({"k": k, "primes": len(ps), "invariant": inv,
                     "violations": len(bad), "xor_add_failures": len(xr)})
    return _frame(rows, index="k")


def anchor_frame(kmax=40):
    """Is each anchor itself prime?  This is where the property lives or dies.

      k_is_pow2    the floor can ONLY be prime here (a^m+1 factors for odd m)
      floor_prime  2^k + 1 is prime  -- a Fermat prime
      ceil_prime   2^(k+1) - 1 is prime -- a Mersenne prime
      both         the shared property
    """
    rows = []
    for k in range(1, kmax + 1):
        f, c = floor_anchor(k), ceil_anchor(k)
        fp, cp = _is_prime(f), _is_prime(c)
        rows.append({"k": k, "k_is_pow2": (k & (k - 1)) == 0,
                     "floor": f if f < 10 ** 9 else None,
                     "floor_prime": fp,
                     "ceil": c if c < 10 ** 9 else None,
                     "ceil_prime": cp, "both": fp and cp})
    return _frame(rows, index="k")


def both_prime(kmax=64):
    """The k where BOTH anchors are prime.  Ends at 16, and that is the question."""
    out = []
    for k in range(1, kmax + 1):
        if _is_prime(floor_anchor(k)) and _is_prime(ceil_anchor(k)):
            out.append({"k": k, "floor": floor_anchor(k), "ceil": ceil_anchor(k),
                        "floor_is": f"F_{k.bit_length()-1}" if (k & (k-1)) == 0 else "-",
                        "ceil_is": f"M_{k+1}"})
    return _frame(out)


def after_16(kmax=1024):
    """The floor's remaining chances past k=16 -- every power of two, all failing.

    2^k + 1 is prime only if k is a power of two, so the Fermat side gets one
    shot per doubling.  F_5..F_11 are all known composite; beyond that the
    numbers are too large for a quick test and are listed as unresolved here.
    """
    rows = []
    k = 1
    while k <= kmax:
        n = (1 << k) + 1
        idx = k.bit_length() - 1
        if k <= 16:
            st = "PRIME" if _is_prime(n) else "composite"
        elif k <= 4096:
            st = "composite" if idx <= 11 else "not tested here"
        else:
            st = "not tested here"
        rows.append({"k": k, "fermat_index": f"F_{idx}",
                     "digits": len(str(n)) if k < 4000 else int(k * 0.30103) + 1,
                     "status": st})
        k <<= 1
    return _frame(rows)


def report(kmax=40):
    print("the two anchors, and the mirror identity")
    print(invariant_check(10))
    print("\nwhere both anchors are prime:")
    print(both_prime(kmax))
    print("\n=> k = 16 is the last.  floor 65537 = F_4, ceiling 131071 = M_17.")
    print("   floor prime requires k a power of two; ceil prime requires k+1")
    print("   a Mersenne exponent.  k = 8 fails the second (511 = 7*73).")
    print("\nthe floor's remaining chances:")
    print(after_16(4096))


if __name__ == "__main__":
    report()
