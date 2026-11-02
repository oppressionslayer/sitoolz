"""What a strong answer to 'write a program' actually looks like.

Given ONLY the terms and no rule, don't write a sequence generator -- write
a search over RULES. This cracks the k=5 sequence in milliseconds, and the
same 60 lines crack the whole family.

The moves, in the order a person would try them:
  1. are the terms squares?           -> take roots
  2. are the roots near squares?      -> get offsets
  3. are the offsets periodic?        -> quasi-polynomial, solves everything
  4. do the offsets come from XOR?    -> the WHY, and the only step that
                                         needs binary thinking
"""

from math import isqrt

TERMS = [1, 16, 25, 144, 441, 784, 1089, 2704, 4761, 7056, 9409, 15376]

is_square = lambda v: v >= 0 and isqrt(v) ** 2 == v


def crack(terms, max_k=1024, max_period=32):
    out = {}

    # 1. squares?
    if all(is_square(t) for t in terms):
        roots = [isqrt(t) for t in terms]
        out["all_squares"] = True
    else:
        roots = list(terms)
        out["all_squares"] = False

    # 2. THE FAST PATH: search over XOR constants directly. No need to
    #    recover m at all -- just ask whether root ^ k is a square for all.
    for k in range(max_k):
        if all(is_square(n ^ k) for n in roots):
            out["xor_k"] = k
            break

    # 3. with k known, m falls out and the quasi-polynomial is visible.
    if "xor_k" in out:
        k = out["xor_k"]
        ms = [isqrt(n ^ k) for n in roots]
        offs = [n - m * m for n, m in zip(roots, ms)]
        out["m_consecutive"] = sorted(ms) == list(range(min(ms), max(ms) + 1))
        out["distinct_offsets"] = sorted(set(offs))
        for p in range(1, max_period + 1):
            table, ok = {}, True
            for m, o in zip(ms, offs):
                if table.setdefault(m % p, o) != o:
                    ok = False
                    break
            if ok and len(table) == p:
                out["period"] = p
                out["table"] = dict(sorted(table.items()))
                break

    # 4. THE SLOW PATH, for comparison: find the quasi-polynomial WITHOUT
    #    ever thinking about XOR. Assign m by rank in sorted order. The head
    #    is scrambled (the bijection is non-monotonic for small m), so drop
    #    leading terms until a clean period appears.
    srt = sorted(roots)
    for drop in range(len(srt) // 2):
        tail = srt[drop:]
        ms2 = [isqrt(n) if n - isqrt(n) ** 2 <= (isqrt(n) + 1) ** 2 - n
               else isqrt(n) + 1 for n in tail]
        offs2 = [n - m * m for n, m in zip(tail, ms2)]
        for p in range(2, max_period + 1):
            table, ok = {}, True
            for m, o in zip(ms2, offs2):
                if table.setdefault(m % p, o) != o:
                    ok = False
                    break
            if ok and len(table) == p and len(tail) >= 2 * p:
                out["qp_without_xor"] = {"dropped": drop, "period": p,
                                         "table": dict(sorted(table.items()))}
                return out
    return out


def generate(k, i):
    """Once you have k, everything is one line."""
    return (i * i ^ k) ** 2


if __name__ == "__main__":
    import time

    t = time.time()
    r = crack(TERMS)
    dt = time.time() - t

    print(f"cracked in {dt*1000:.2f} ms\n")
    for key in ("all_squares", "m_consecutive", "distinct_offsets",
                "period", "table", "xor_k"):
        if key in r:
            print(f"  {key:20} {r[key]}")

    k = r["xor_k"]
    print(f"\n  rule: x = (m^2 ^ {k})^2")
    print(f"  reproduces the input: "
          f"{sorted(generate(k, m) for m in range(12)) == sorted(TERMS)}")
    print(f"  next term:            {generate(k, 12)}")
    print(f"  10^25-th term has     {len(str(generate(k, 10**25)))} digits")
