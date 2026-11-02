"""ntprobe -- probe how the square-XOR machines index a value, WITH an honest control.

A tempting headline: "the first 1024 Rule-30 centre bits land EXACTLY on the m2
machine!"  It is true -- and it is not special.  The m2 machine's index of ``v*v``
is ``v`` itself, so the correction ``K == 0`` for *every* input, Rule-30 or random
alike: squaring then indexing round-trips exactly on the plain-square (m2) machine
by construction.  The fixed-mask machines (xor2 / xor5 / m2+2) just carry a constant
tiny correction (K_pop 1-2) set by their mask, the same for any v.  The one genuinely
input-dependent index here is ``nibrev``.

Moral (the Monty-Hall move): before calling a coincidence deep, run the control on a
random value of the same size.  If the "miracle" survives, it was the method, not the
data.  ``report()`` shows the observation and the control side by side.

    from sitoolz import ntprobe
    ntprobe.report(1024)
"""
import random
from sitoolz import ntfun


def rule30_centre_value(gens: int = 1024) -> int:
    """The Rule 30 centre column for ``gens`` generations, read as one integer.

    Seed is a single 1; generation 0 (the seed) is the most-significant bit.
    """
    W = 2 * gens + 5
    row = [0] * W
    row[W // 2] = 1
    v = 0
    for _ in range(gens):
        v = (v << 1) | row[W // 2]
        row = [(row[i - 1] ^ (row[i] | row[i + 1])) if 0 < i < W - 1 else 0
               for i in range(W)]
    return v


def machine_exactness(v: int, machines=None) -> dict:
    """Index ``v*v`` on each square-XOR machine and report the leftover correction K.

    Returns ``{machine: {"K_bits", "K_pop", "exact"}}`` where ``exact`` (K == 0) means
    ``v*v`` sits exactly on that machine's lattice.  Wraps :func:`ntfun.machine_index`.
    """
    res = ntfun.machine_index(v * v, machines)
    out = {}
    for name, (K, i, m) in res.items():
        K = int(K)
        out[name] = {"K_bits": K.bit_length(),
                     "K_pop": bin(K).count("1"),
                     "exact": K == 0}
    return out


def report(gens: int = 1024, seed: int = 7):
    """Print the Rule-30 observation beside a random control of the same bit length.

    The take-away: whatever is ``exact`` for Rule 30 is exact for the random control
    too -- so it is a property of the machine, not of Rule 30.  Returns
    ``(rule30_result, random_result)``.
    """
    rng = random.Random(seed)
    r30 = rule30_centre_value(gens)
    rnd = rng.getrandbits(r30.bit_length())
    a, b = machine_exactness(r30), machine_exactness(rnd)

    print(f"Rule 30 centre ({gens} gens, {r30.bit_length()} bits)  vs  random control:")
    print(f"  {'machine':8} | {'rule30 K_pop':>12} {'exact':>6} | {'random K_pop':>12} {'exact':>6}")
    print("  " + "-" * 54)
    for name in a:
        print(f"  {name:8} | {a[name]['K_pop']:>12} {str(a[name]['exact']):>6} "
              f"| {b[name]['K_pop']:>12} {str(b[name]['exact']):>6}")
    same = all(a[n]['exact'] == b[n]['exact'] for n in a)
    print("\n  Control: exact-flags match Rule-30 vs random for every machine:", same)
    print("  => m2's exactness is the machine (plain-square round-trip), not Rule 30.")
    return a, b


def _clip(v, clip):
    """Decimal string of v, clipped to ~`clip` chars as head…tail so it fits a screen."""
    s = str(v)
    if clip is None or len(s) <= clip:
        return s
    keep = max(4, clip - 3)
    h, t = (keep + 1) // 2, keep // 2
    return f"{s[:h]}…{s[-t:]}"


def walk_down(N, clip=64):
    """Walk a number DOWN its index tower and return the walk as a DataFrame.

    Each rung is one big-integer square root: ``N -> isqrt(N) -> isqrt(that) -> … -> 1``.
    The bit length halves every step, so ANY N reaches its 1-bit seed in about
    ``log2(bit_length)`` rungs -- and the funny part is that each rung IS the index of
    the one above it (the member's index is its square root).  So this is literally
    "index the index, all the way down."

    Big values are clipped to ``clip`` characters (head…tail) so you can read them on
    screen; pass ``clip=None`` for the full digits.  Columns:

        bits      bit length of this rung
        digits    decimal length
        rem       v - isqrt(v)**2  -- the part shed descending one step
        exact_sq  whether this rung is itself a perfect square (rem == 0)
        value     the rung, clipped for display

    Example::

        from sitoolz import ntprobe
        ntprobe.walk_down(ntprobe.rule30_centre_value(1024) ** 2)   # 617 digits -> 1
    """
    import pandas as pd
    from gmpy2 import mpz, isqrt

    rows, v, level = [], mpz(N), 0
    while True:
        r = isqrt(v)
        rem = v - r * r
        rows.append({"level": level, "bits": v.bit_length(), "digits": len(str(v)),
                     "rem": _clip(int(rem), clip), "exact_sq": rem == 0,
                     "value": _clip(v, clip)})
        if v <= 1:
            break
        v, level = r, level + 1
    return pd.DataFrame(rows).set_index("level")


def walk_down_demo(gens=1024, clip=64):
    """The showcase walk: the first `gens` Rule-30 centre bits, SQUARED, walked all the
    way down its index tower.  Row 0 is the 2*gens-bit member; row 1 is its index (the
    Rule-30 centre value itself); each row after is the index of the one above."""
    return walk_down(rule30_centre_value(gens) ** 2, clip=clip)
