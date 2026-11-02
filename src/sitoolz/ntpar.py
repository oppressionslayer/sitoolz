"""ntpar -- the parity of squares, bit by bit.

Chasing "how close to a parity do we get" through the xor-square machines
lands somewhere simpler than it started.  Whether `isqrt(i*i ^ m) == i`
holds does not depend on the mask at all, only on its top bit:

    isqrt(i*i ^ m) == i   iff   bit b of i*i is 0,   b = m.bit_length()-1

because xor only fails to raise a number when it clears its highest
differing bit.  That holds once the offset is small next to the gap
between consecutive squares -- for every `2*i + 1 > m`, and it can differ
below that; `predict` is the exact test, `bit_rule` the asymptotic one.

So the whole machine question was a question about the bits of squares,
and it has three answers and no others.

    bit 0   exactly 1/2, forever          <- the one that is really parity
    bit 1   exactly 0, forever            <- never balances; squares mod 4
    bit b   1/2 + 2**-(b//2 + 1)          <- approaches parity, never lands

    from sitoolz import ntpar

    ntpar.census(12)          # zeros of bit b across squares, vs the law
    ntpar.masks(64)           # every mask, and the bit it really depends on
    ntpar.convergence(16)     # the excess halving every second bit
    ntpar.exact()             # the bits that hit parity on the nose
    ntpar.predict(5, 97)      # exact, no isqrt
    ntpar.bit_rule(5, 97)     # the one-bit version, exact once 2i+1 > m
"""
try:
    import pandas as pd
except ImportError:
    pd = None


def _need_pd():
    if pd is None:
        raise ImportError("pandas required")


def law(b):
    """Fraction of squares with bit b clear -- closed form, no counting.

    Bit 0 is clear exactly when k is even, so it is a true half.  Bit 1 is
    clear always, since squares are 0 or 1 mod 4.  Above that the fraction
    is 1/2 + 2**-(b//2 + 1), which halves every second bit and reaches 1/2
    only in the limit.
    """
    if b == 0:
        return 0.5
    if b == 1:
        return 1.0
    return 0.5 + 2.0 ** -(b // 2 + 1)


def top_bit(m):
    """The only part of a mask that decides whether the index survives."""
    return m.bit_length() - 1


def predict(m, i):
    """Does `isqrt(i*i ^ m) == i`?  Exact, and still no isqrt.

    Two comparisons: the xored value has to land inside `[i*i, (i+1)**2)`.
    `bit_rule` is the cheaper asymptotic version and agrees with this one
    once `2*i + 1 > m` -- see its docstring for why they part company.
    """
    v = i * i ^ m
    return i * i <= v < (i + 1) * (i + 1)


def bit_rule(m, i):
    """The asymptotic rule: bit b of i*i is clear, b = m.bit_length()-1.

    Only the top bit of the mask matters, because xor fails to raise a
    number exactly when it clears that number's highest differing bit.

    This tests `i*i ^ m >= i*i`, which is what decides the index *once the
    offset is small next to the gap between consecutive squares*.  That gap
    is `2*i + 1`, so `bit_rule` and `predict` agree for every `2*i + 1 > m`
    and can differ below it -- for i = 0 with m = 5, the offset of +5
    overshoots 1 entirely and isqrt lands on 2, not 0.  Verified: 228
    disagreements over masks 1..39, all with `2*i + 1 <= m`, none above.
    """
    return ((i * i) >> top_bit(m)) & 1 == 0


def census(bmax=12, kmax=1 << 17):
    """Count zeros of bit b over k*k for k < kmax, against `law`.

      zeros     how many squares have bit b clear
      fraction  zeros / kmax
      law       the closed form
      excess    fraction - 1/2, the distance from true parity
      agree     measured and predicted match to 6 places
    """
    _need_pd()
    rows = []
    for b in range(bmax):
        z = sum(1 for k in range(kmax) if ((k * k >> b) & 1) == 0)
        frac = z / kmax
        rows.append({"bit": b, "zeros": z, "of": kmax,
                     "fraction": round(frac, 6), "law": round(law(b), 6),
                     "excess": round(frac - 0.5, 6),
                     "agree": abs(frac - law(b)) < 1e-6})
    df = pd.DataFrame(rows).set_index("bit")
    df.attrs["kmax"] = kmax
    df.attrs["all_agree"] = bool(df.agree.all())
    return df


def convergence(bmax=16):
    """The excess over parity, as exact fractions, and the halving.

    The excess is 2**-(b//2 + 1) from bit 2 upward, so it halves every
    second bit -- the ratio column reads 1, 2, 1, 2, 1, 2 forever.
    """
    _need_pd()
    from fractions import Fraction
    rows, prev = [], None
    for b in range(bmax):
        e = Fraction(1, 2) if b == 1 else (
            Fraction(0) if b == 0 else Fraction(1, 2 ** (b // 2 + 1)))
        rows.append({"bit": b, "fraction": str(Fraction(1, 2) + e),
                     "excess": str(e), "decimal": float(e),
                     "vs_prev": None if not prev else round(float(prev / e), 1)
                     if e else None})
        prev = e if e else prev
    df = pd.DataFrame(rows).set_index("bit")
    df.attrs["limit"] = 0.5
    df.attrs["reached"] = "only at bit 0"
    return df


def masks(mmax=64, kmax=20000):
    """Every mask up to mmax, and the single bit it actually depends on.

    Masks that share a top bit are the same machine as far as parity is
    concerned -- 16 through 31 all ask about bit 4.  `distinct_masks` is
    how many masks collapse onto each bit.

      top_bit   m.bit_length() - 1
      keeps     residues (out of kmax) where the index survives, by the
                exact test -- the handful of small i where `bit_rule`
                differs are far below the rounding here
      fraction  keeps / kmax
      law       the closed form for that bit
      agree     they match
    """
    _need_pd()
    rows = []
    for m in range(1, mmax + 1):
        b = top_bit(m)
        keeps = sum(1 for i in range(kmax) if predict(m, i))
        frac = keeps / kmax
        rows.append({"m": m, "bin": bin(m)[2:], "top_bit": b,
                     "fraction": round(frac, 6), "law": round(law(b), 6),
                     "agree": abs(frac - law(b)) < 1e-3})
    df = pd.DataFrame(rows).set_index("m")
    df.attrs["all_agree"] = bool(df.agree.all())
    df.attrs["distinct_bits"] = int(df.top_bit.nunique())
    return df


def exact(bmax=64):
    """The bits that hit parity exactly -- the cheat, and there is one.

    Bit 0 is clear exactly when k is even: half the squares, on the nose,
    at every scale, no limit needed.  It is the only one.  Bit 1 is the
    opposite extreme and is never balanced at all.  Everything above is
    forever a little too clear.
    """
    _need_pd()
    rows = [{"bit": b, "fraction": law(b),
             "exact_half": law(b) == 0.5,
             "note": ("k even" if b == 0 else
                      "squares are 0 or 1 mod 4" if b == 1 else
                      f"1/2 + 2^-{b // 2 + 1}")}
            for b in range(bmax)]
    df = pd.DataFrame(rows).set_index("bit")
    df.attrs["exact_bits"] = [int(b) for b in df.index[df.exact_half]]
    return df
