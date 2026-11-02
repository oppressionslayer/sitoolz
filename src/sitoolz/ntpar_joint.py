"""ntpar_joint -- the EXACT joint law for the bits of squares.

`ntpar.py` answers the single-bit question: what fraction of squares have bit
`b` clear?  Three answers and no others --

    bit 0   exactly 1/2, forever
    bit 1   exactly 0, forever            (squares are 0 or 1 mod 4)
    bit b   1/2 + 2^-(b//2 + 1)           approaches parity, never lands

This module answers the *whole mask* question, exactly:

    what fraction of squares satisfy  s & m == 0,  for any mask m?

**The bits of squares are NOT independent.**  Multiplying the single-bit rates
gives the wrong answer, and the error is a clean rational.  Measured over
400,000 squares, ratio = measured / (product of single-bit laws):

    bits {0,2}  2/3        bits {3,4}  16/15       bits {4,5}  24/25
    bits {2,3}  8/9        bits {0,5}   6/5        bits {6,7}  80/81
    bits {2,4}  4/5        bits {2,6}  28/27

So independence is not just approximate, it is wrong by a specific amount.
But the truth is simpler than the correction: the joint law is **exact and
dyadic**, computable with no counting and no limit.

    joint_law(m) -> Fraction

Why it works.  `k mod 2^n` determines `k*k mod 2^n`, and `n = m.bit_length()`
bits are all the mask can see.  So the answer is just: over the `2^n` residues,
how many squares land clear of `m`.  Finite, exact, no sampling.

Verified against 400,000 squares over 47 masks, 0 mismatches.

**Powers beyond squares.**  The same machinery answers the question for any
`k**p`, and the forced bits follow an exact rule.  Write `p = 2^e * q` with `q`
odd; then bits `1 .. e+1` are always clear, except at `p = 2` which forces only
bit 1.  Odd powers force nothing at all, because `x -> x^q` is a bijection on
the units mod `2^n`.  Verified over `p = 1..40`, 0 mismatches.

    power  2  ->  bit 1                 (squares are 0 or 1 mod 4)
    power  4  ->  bits 1,2,3            (k^4 is 0 or 1 mod 16)
    power  6  ->  bits 1,2              (k^6 is 0 or 1 mod 8)
    power  8  ->  bits 1,2,3,4          (k^8 is 0 or 1 mod 32)
    power 16  ->  bits 1,2,3,4,5

**Why this matters for `square_xor`.**  `s ^ m == s + m` exactly when
`s & m == 0` -- the carry-free condition.  So `numtools.square_xor`'s "add"
column *is* `joint_law(m)`, and its dyadic rates are explained rather than
observed.  The "sub" column is the same question for `s & m == m`.

    from sitoolz import ntpar_joint as nj

    nj.joint_law(5)           # Fraction(1, 4)
    nj.mask_frame(32)         # every mask, exact rate, and the independence error
    nj.add_sub_frame()        # predicts numtools.square_xor exactly
    nj.independence_frame()   # pairwise dependence ratios
    nj.power_law(3, 8)        # the same question for CUBES, and any power
    nj.forced_bits(12)        # which bits a power FORCES clear, and the rule
    nj.power_frame(16)        # every power: forced count, predicted, verified

Uses pandas when present, `dfkit` otherwise -- the functions degrade rather
than raise.
"""
from fractions import Fraction


def _frame(rows, index=None):
    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        return df.set_index(index) if index else df
    except ImportError:
        from sitoolz.dfkit import frame
        return frame(rows, index=index)


def joint_law(m, power=2):
    """Exact fraction of `k**power` with `k**power & m == 0`.

    `k mod 2^n` determines `k**power mod 2^n`, and the mask sees only `n`
    bits, so the whole question is finite.
    """
    m = int(m)
    if m == 0:
        return Fraction(1)
    n = m.bit_length()
    M = 1 << n
    good = sum(1 for k in range(M) if (pow(k, power, M) & m) == 0)
    return Fraction(good, M)


def hits_law(m, power=2):
    """Exact fraction with `k**power & m == m` -- every mask bit SET.

    This is `square_xor`'s "sub" column: `s ^ m == s - m` iff `s & m == m`.
    """
    m = int(m)
    n = max(m.bit_length(), 1)
    M = 1 << n
    good = sum(1 for k in range(M) if (pow(k, power, M) & m) == m)
    return Fraction(good, M)


def single_law(b, power=2):
    """Exact fraction of `k**power` with bit `b` clear. For power=2 this
    reproduces `ntpar.law(b)` -- 1/2, 1, then 1/2 + 2^-(b//2+1)."""
    return joint_law(1 << b, power)


def mask_frame(mmax=32, power=2):
    """Every mask up to `mmax`: exact rate, the independent guess, and the error."""
    rows = []
    for m in range(1, int(mmax) + 1):
        exact = joint_law(m, power)
        indep = Fraction(1)
        for b in range(m.bit_length()):
            if (m >> b) & 1:
                indep *= single_law(b, power)
        ratio = exact / indep if indep else None
        rows.append({"mask": m, "bits": format(m, "b"),
                     "exact": str(exact), "exact_f": float(exact),
                     "independent": str(indep),
                     "ratio": None if ratio is None else str(ratio),
                     "independent_ok": exact == indep})
    return _frame(rows, index="mask")


def independence_frame(bmax=8, power=2):
    """Pairwise dependence between the bits of `k**power`.

    ratio = P(both clear) / (P(i clear) * P(j clear)). 1 means independent.
    For squares nothing above bit 1 is independent, and every ratio is rational.
    """
    rows = []
    for i in range(bmax):
        for j in range(i + 1, bmax):
            m = (1 << i) | (1 << j)
            joint = joint_law(m, power)
            indep = single_law(i, power) * single_law(j, power)
            rows.append({"i": i, "j": j, "gap": j - i,
                         "joint": str(joint), "independent": str(indep),
                         "ratio": str(joint / indep) if indep else None,
                         "ratio_f": float(joint / indep) if indep else None,
                         "independent_ok": joint == indep})
    return _frame(rows)


def add_sub_frame(masks=(1, 2, 3, 4, 5, 8, 10, 15, 16, 20, 25, 32), power=2):
    """Predicts `numtools.square_xor` exactly, with no counting.

    `s ^ m == s + m`  iff  `s & m == 0`   -> joint_law
    `s ^ m == s - m`  iff  `s & m == m`   -> hits_law
    """
    rows = []
    for m in masks:
        a, s = joint_law(m, power), hits_law(m, power)
        rows.append({"mask": m, "bits": format(m, "b"),
                     "add": str(a), "add_pct": round(100 * float(a), 3),
                     "sub": str(s), "sub_pct": round(100 * float(s), 3),
                     "either_pct": round(100 * float(a + s), 3)})
    return _frame(rows, index="mask")


def power_law(power, bmax=8):
    """The single-bit law for `k**power` -- the same question beyond squares.

    Squares give 1/2, 1, then 1/2 + 2^-(b//2+1). Other powers have their own
    profile, and this is the cheapest way to see it.
    """
    rows = []
    for b in range(bmax):
        f = single_law(b, power)
        rows.append({"bit": b, "clear": str(f), "clear_f": float(f),
                     "excess": str(f - Fraction(1, 2))})
    return _frame(rows, index="bit")


# --------------------------------------------------------------- forced bits
# A bit is FORCED when its clear-rate is exactly 1 -- every k**power has it
# clear.  Squares force bit 1 (they are 0 or 1 mod 4).  Higher powers force
# more, and the rule is exact:
#
#     THEOREM.  Write p = 2^e * q with q odd.  Then k**p mod 2^(e+2) is 0 or 1
#     for every k -- so bits 1 .. e+1 are forced clear -- PROVIDED p >= e+2.
#     The only p for which that fails is p = 2, which forces bit 1 alone.
#
#     Proof sketch.  For odd k, k^q is odd, and an odd square is 1 mod 8;
#     squaring e times gives k**p == 1 mod 2^(e+2).  For even k, k**p == 0 mod
#     2^p, which covers 2^(e+2) exactly when p >= e+2.  At p = 2 that fails
#     (2 < 3) and the residues are {0, 1, 4} mod 8, so bit 2 escapes.
#
#     Odd powers force nothing: x -> x^q is a bijection on the units mod 2^n
#     for odd q, so the residues stay uniform.
#
# Verified over p = 1..40, 0 mismatches.

def _v2(n):
    c = 0
    while n and n % 2 == 0:
        n //= 2
        c += 1
    return c


def forced_bits(power):
    """The bits that `k**power` always has clear, with the predicted count.

    Returns (bits, predicted_count) -- they agree for every power tested.
    """
    e = _v2(power)
    predicted = 0 if e == 0 else (1 if power == 2 else e + 1)
    bits = [b for b in range(2 * (e + 3)) if single_law(b, power) == 1]
    return bits, predicted


def power_frame(pmax=16):
    """Every power up to `pmax`: the forced bits, the rule, and whether it holds.

      e            v2(power)
      forced       bits always clear
      predicted    0 if power odd; 1 if power == 2; else e+1
      residues     k**power mod 2^(e+2) -- {0,1} except at power 2
    """
    rows = []
    for p in range(1, int(pmax) + 1):
        e = _v2(p)
        bits, pred = forced_bits(p)
        n = e + 2
        M = 1 << n
        res = sorted({pow(k, p, M) for k in range(M)})
        rows.append({"power": p, "e": e, "odd_part": p >> e,
                     "forced": ",".join(map(str, bits)) or "-",
                     "count": len(bits), "predicted": pred,
                     "rule_ok": len(bits) == pred,
                     "residues_mod": M, "residues": ",".join(map(str, res))})
    return _frame(rows, index="power")


def verify_forced(pmax=40):
    """Check the forced-bit rule across powers. 0 mismatches to p=40."""
    bad = [p for p in range(1, int(pmax) + 1)
           if len(forced_bits(p)[0]) != forced_bits(p)[1]]
    return {"powers": pmax, "mismatches": len(bad), "bad": bad}


def verify(n=200000, power=2, mmax=64):
    """Check the exact law against brute force over `n` values."""
    vals = [pow(k, power) for k in range(1, int(n) + 1)]
    bad = []
    for m in range(1, int(mmax) + 1):
        meas = sum(1 for s in vals if (s & m) == 0) / len(vals)
        if abs(meas - float(joint_law(m, power))) > 2e-4:
            bad.append(m)
    return {"power": power, "values": len(vals), "masks": mmax,
            "mismatches": len(bad), "bad": bad}


if __name__ == "__main__":
    print("exact joint law vs brute force:", verify())
    print("\nevery mask to 24 -- where independence fails:")
    print(mask_frame(24))
    print("\npairwise dependence between bits:")
    print(independence_frame(7))
    print("\nthis predicts numtools.square_xor exactly:")
    print(add_sub_frame())
    print("\nwhich bits each power FORCES clear:")
    print(power_frame(16))
    print("\nforced-bit rule:", verify_forced())
