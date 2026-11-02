"""squarewalk -- the square family, walked by the parity rule instead of isqrt.

`ntfun.square_family` addresses the numbers that are square *and* whose root
xored with `m` is square:  x = (i*i ^ m)**2, indexed by i.  It finds them in
closed form, and asks `isqrt` whether the xor kept the root inside index i.

`ntpar` answered that question without the square root:

    isqrt(i*i ^ m) == i   iff   the xor did not clear the top bit of i*i

so the whole thing is one bit test.  This file is what that buys.  It is new
work built on top of both -- neither `ntfun` nor `ntpar` is touched.

    from sitoolz import squarewalk as sw

    sw.ladder(5)              # the family with a stays/jumps column
    sw.survivors(5, 0, 64)    # the indices the xor keeps
    sw.collapse(64)           # 64 masks are 6 questions -- computed, not asserted
    sw.pattern(32, 64)        # the survivor grid, as a picture
    sw.chain(5, 3)            # iterate i -> i*i ^ m and see where it goes
    sw.bench(5)               # predict vs isqrt, timed
    sw.agree(64, 4096)        # predict against the real isqrt, exhaustively

Frames: `ladder_frame`, `collapse_frame`, `chain_frame`, `mask_frame`.
"""
from math import isqrt
from sitoolz import ntpar

predict = ntpar.predict
bit_rule = ntpar.bit_rule
law = ntpar.law
top_bit = ntpar.top_bit


# ------------------------------------------------------------------- the map
def root(m, i):
    """The root of the family member at index i:  i*i ^ m."""
    return i * i ^ m


def member(m, i):
    """The family member itself:  (i*i ^ m)**2."""
    r = i * i ^ m
    return r * r


def stays(m, i):
    """Did the xor keep the root inside index i?  One bit test, no isqrt."""
    return predict(m, i)


def jumped_to(m, i):
    """Which index the root actually landed in.  isqrt only for the report."""
    return isqrt(i * i ^ m)


# --------------------------------------------------------------- the survivors
def survivors(m, lo=0, hi=64):
    """Indices in [lo, hi) that the xor keeps -- the set `predict` cuts out."""
    return [i for i in range(lo, hi) if predict(m, i)]


def density(m, n=1 << 14):
    """Measured survivor density against the closed form for this mask's bit."""
    keeps = sum(1 for i in range(n) if predict(m, i))
    b = top_bit(m)
    return {"m": m, "top_bit": b, "measured": keeps / n, "law": law(b),
            "agree": abs(keeps / n - law(b)) < 1e-3}


def signature(m, n=64):
    """The survivor pattern as a bitmask -- two masks match iff this matches."""
    s = 0
    for i in range(n):
        if predict(m, i):
            s |= 1 << i
    return s


def collapse(mmax=64, n=64, lo=None):
    """Group masks by identical survivor patterns.

    `ntpar` says only the top bit of the mask matters -- but that is the
    *asymptotic* statement, true once `2*i + 1 > m`.  Below that boundary the
    exact test disagrees with it, and every mask is an individual.

    So the window matters, and getting it wrong is the easy mistake: scan
    i < 64 with masks up to 64 and most of the window sits *below* the
    boundary, so 64 masks look like 23 groups instead of 7.  `lo` defaults to
    just above the boundary for every mask in range, where the collapse is
    real.  Pass lo=0 to see the exceptions instead.
    """
    if lo is None:
        lo = mmax // 2 + 1
    groups = {}
    for m in range(1, mmax + 1):
        sig = 0
        for k, i in enumerate(range(lo, lo + n)):
            if predict(m, i):
                sig |= 1 << k
        groups.setdefault(sig, []).append(m)
    out = []
    for sig, ms in sorted(groups.items(), key=lambda kv: min(kv[1])):
        bits = {top_bit(m) for m in ms}
        out.append({"masks": ms, "n_masks": len(ms), "top_bits": sorted(bits),
                    "one_top_bit": len(bits) == 1,
                    "survivors": bin(sig).count("1"), "of": n,
                    "window": (lo, lo + n)})
    return out


def boundary_effect(mmax=64, n=64):
    """How many groups you get below the boundary versus above it.

    Below: every mask is nearly its own case.  Above: one group per top bit,
    which is the collapse `ntpar` states.
    """
    return {"below_boundary": len(collapse(mmax, n, lo=0)),
            "above_boundary": len(collapse(mmax, n)),
            "distinct_top_bits": len({top_bit(m) for m in range(1, mmax + 1)})}


def pattern(mmax=32, n=64, mark="#", blank="."):
    """The survivor grid as text: a row per mask, a column per index."""
    lines = []
    for m in range(1, mmax + 1):
        row = "".join(mark if predict(m, i) else blank for i in range(n))
        lines.append(f"m={m:3} b={top_bit(m):2} |{row}|")
    return "\n".join(lines)


# ------------------------------------------------------------------- walking
def chain(m, start, steps=12):
    """Iterate i -> i*i ^ m and record what each step did.

    This is the walk the family invites: the root of one index is a number,
    so read it as the next index and go again.  `stays` is the parity rule,
    so the whole walk runs without a square root.
    """
    out, i, seen = [], start, {}
    for k in range(steps):
        if i in seen:
            out.append({"step": k, "i": i, "root": None, "stays": None,
                        "note": f"cycle back to step {seen[i]}"})
            break
        seen[i] = k
        r = root(m, i)
        out.append({"step": k, "i": i, "root": r, "member": r * r,
                    "stays": stays(m, i), "landed_in": jumped_to(m, i),
                    "note": ""})
        if r == i:
            out[-1]["note"] = "fixed point"
            break
        i = r
    return out


def fixed_points(m, hi=1 << 16):
    """Indices with i*i ^ m == i -- the walk stops dead."""
    return [i for i in range(hi) if (i * i ^ m) == i]


def ladder(m=5, lo=0, hi=12):
    """The family in index order, with the parity verdict alongside."""
    out, prev = [], None
    for i in range(lo, hi):
        r = root(m, i)
        x = r * r
        out.append({"i": i, "i_sq": i * i, "root": r, "member": x,
                    "offset": r - i * i, "stays": stays(m, i),
                    "landed_in": jumped_to(m, i),
                    "gap": None if prev is None else x - prev})
        prev = x
    return out


# ------------------------------------------------------- checking and timing
def agree(mmax=64, imax=4096):
    """`predict` against a real isqrt, and `bit_rule` against `predict`.

    The documented boundary is 2*i + 1 > m: above it all three agree, below
    it `bit_rule` can differ.  This checks both claims exhaustively.
    """
    pred_bad, rule_bad, rule_bad_above = 0, 0, 0
    total = 0
    for m in range(1, mmax + 1):
        for i in range(imax):
            total += 1
            truth = isqrt(i * i ^ m) == i
            if predict(m, i) != truth:
                pred_bad += 1
            if bit_rule(m, i) != truth:
                rule_bad += 1
                if 2 * i + 1 > m:
                    rule_bad_above += 1
    return {"pairs": total,
            "predict_mismatches": pred_bad,
            "bit_rule_mismatches": rule_bad,
            "bit_rule_mismatches_above_2i+1>m": rule_bad_above}


def bench(m=5, n=200_000, base=0):
    """The three tests, same answers, timed.  `base` sets the index scale.

    In CPython `math.isqrt` is C and `predict` is big-integer arithmetic, so
    isqrt wins on small indices.  predict overtakes it once the numbers get
    big enough that the square root costs more than a multiply -- around
    10^12 -- and `bit_rule` is cheaper again, being one multiply and a shift.
    """
    import time
    rng = range(base, base + n)
    t0 = time.perf_counter()
    a = sum(1 for i in rng if isqrt(i * i ^ m) == i)
    t1 = time.perf_counter()
    b = sum(1 for i in rng if predict(m, i))
    t2 = time.perf_counter()
    c = sum(1 for i in rng if bit_rule(m, i))
    t3 = time.perf_counter()
    return {"n": n, "m": m, "base": base, "digits": len(str(base + n)),
            "predict_agrees": a == b,
            "isqrt_sec": round(t1 - t0, 5),
            "predict_sec": round(t2 - t1, 5),
            "bit_rule_sec": round(t3 - t2, 5),
            "predict_speedup": round((t1 - t0) / (t2 - t1), 2),
            "bit_rule_speedup": round((t1 - t0) / (t3 - t2), 2),
            "bit_rule_count_differs_by": abs(a - c)}


def scale_table(m=5, n=3000, exponents=(3, 6, 12, 30, 100, 300, 1000)):
    """Where each test wins, across index scales."""
    return [bench(m, n, base=10 ** e) | {"scale": f"10^{e}"} for e in exponents]


# ---------------------------------------------------------------- dataframes
def ladder_frame(m=5, lo=0, hi=16):
    from sitoolz.dfkit import frame
    return frame(ladder(m, lo, hi), index="i")


def collapse_frame(mmax=64, n=64):
    from sitoolz.dfkit import frame
    rows = []
    for g in collapse(mmax, n):
        rows.append({"top_bits": ",".join(map(str, g["top_bits"])),
                     "n_masks": g["n_masks"], "survivors": g["survivors"],
                     "of": g["of"],
                     "density": round(g["survivors"] / g["of"], 4),
                     "law": round(law(g["top_bits"][0]), 4),
                     "masks": ",".join(map(str, g["masks"][:10]))
                              + (" …" if g["n_masks"] > 10 else "")})
    return frame(rows)


def chain_frame(m=5, start=3, steps=12):
    from sitoolz.dfkit import frame
    return frame(chain(m, start, steps))


def scale_frame(m=5, n=3000):
    from sitoolz.dfkit import frame
    rows = [{"scale": r["scale"], "digits": r["digits"],
             "isqrt_sec": r["isqrt_sec"], "predict_sec": r["predict_sec"],
             "bit_rule_sec": r["bit_rule_sec"],
             "predict_speedup": r["predict_speedup"],
             "bit_rule_speedup": r["bit_rule_speedup"],
             "agrees": r["predict_agrees"]} for r in scale_table(m, n)]
    return frame(rows, index="scale")


def mask_frame(mmax=64, n=1 << 12):
    from sitoolz.dfkit import frame
    return frame([density(m, n) for m in range(1, mmax + 1)], index="m")


if __name__ == "__main__":
    print("the family for m=5, with the parity verdict:")
    for r in ladder(5, 0, 8):
        print(f"  i={r['i']:2}  {r['i_sq']:4} ^5 = {r['root']:4}  "
              f"member {r['member']:8}  stays={str(r['stays']):5}  "
              f"landed in {r['landed_in']}")
    print("\n64 masks grouped by survivor pattern, above the 2i+1 > m boundary:")
    for g in collapse(64):
        print(f"  top bit {g['top_bits']}: {g['n_masks']:2} masks, "
              f"{g['survivors']}/{g['of']} survive   one_top_bit={g['one_top_bit']}")
    print("\n  window matters:", boundary_effect(64))
    print("\npredict vs isqrt:", bench(5, 100_000))
    print("exhaustive check:", agree(32, 2048))
