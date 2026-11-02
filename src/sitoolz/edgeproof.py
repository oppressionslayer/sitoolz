"""edgeproof -- the period-spectrum account of rule 30's edge law.

Reproduces every claim in EDGE_LAW.md.  Stdlib only.

    python3 edgeproof.py            # run the whole battery

The claims, in order:

    C1  every bit of the edge is exactly periodic in g, period a power of 2
    C2  the period doubles from bit b-1 to bit b iff the driver
        d = bit_{b-1} | bit_{b-2} has odd weight over one of its periods
    C3  bits whose period divides 2^v are all 1 at g = -1 mod 2^v (in phase)
    C4  the first bit whose period exceeds 2^v is 0 there (the run stops)
    C5  edgelaw.TABLE[v] == #{b : period(b) divides 2^v}
"""
from collections import Counter


# ------------------------------------------------------------------ evolve
def edge_rows(G, W=128):
    """Generations 0..G-1 of the edge, in the frame co-moving with the cone.

    E'[j] = E[j] XOR (E[j-1] OR E[j-2]).  The cone edge never consults the
    interior, so this is exact, not an approximation of the full automaton.

    Bit j depends only on bits j-1, j-2, so information flows strictly from
    low bits to high.  Masking at W therefore never corrupts anything below
    W -- a narrow window is exact for the bits it keeps, at any generation.
    """
    m = (1 << W) - 1
    E = 1
    rows = [E]
    for _ in range(G - 1):
        E = (E ^ ((E << 1) | (E << 2))) & m
        rows.append(E)
    return rows


_COLCACHE = {}


def column(rows, b):
    """Bit b of the edge, as a function of generation.  Memoised."""
    key = (id(rows), b)
    if key not in _COLCACHE:
        _COLCACHE[key] = [(r >> b) & 1 for r in rows]
    return _COLCACHE[key]


# ------------------------------------------------------------------ periods
def period_exp(seq, kmax):
    """Least k with seq period 2^k, or None if none <= 2^kmax fits."""
    for k in range(kmax + 1):
        p = 1 << k
        if all(seq[i] == seq[i + p] for i in range(len(seq) - p)):
            return k
    return None


def spectrum(rows, nbits, kmax):
    """p(b) for b < nbits.  None where the sample is too short to trust."""
    return [period_exp(column(rows, b), kmax) for b in range(nbits)]


# ------------------------------------------------------------------- claims
def c1_periodic(P):
    """C1: every resolved bit has a power-of-two period (None = unresolved)."""
    resolved = [(b, p) for b, p in enumerate(P) if p is not None]
    return {"resolved_bits": len(resolved), "max_exp": max(p for _, p in resolved)}


def c2_weight_rule(rows, P):
    """C2: p(b) == p(b-1) + 1  iff  the driver has odd weight over its period.

    bit_b is the running XOR-sum of d = bit_{b-1} | bit_{b-2}.  A running
    XOR-sum of a period-P signal keeps period P when the signal's weight
    over one period is even, and doubles to 2P when it is odd.
    """
    out = []
    for b in range(2, len(P)):
        if P[b] is None or P[b - 1] is None or P[b - 2] is None:
            break
        d = [x | y for x, y in zip(column(rows, b - 1), column(rows, b - 2))]
        Pd = 1 << max(P[b - 1], P[b - 2])
        weight_odd = sum(d[:Pd]) % 2 == 1
        doubled = P[b] == max(P[b - 1], P[b - 2]) + 1
        out.append({"bit": b, "driver_exp": max(P[b - 1], P[b - 2]),
                    "weight_odd": weight_odd, "doubled": doubled,
                    "agree": weight_odd == doubled})
    return out


def c3_in_phase(rows, P):
    """C3: bit b is 1 at every g == -1 mod 2^p(b)."""
    G = len(rows)
    out = []
    for b, p in enumerate(P):
        if p is None:
            break
        step = 1 << p
        gs = range(step - 1, G, step)
        ones = all(column(rows, b)[g] == 1 for g in gs)
        out.append({"bit": b, "p": p, "always_one_at_-1_mod_2^p": ones,
                    "samples": len(gs)})
    return out


def c4_run_stops(rows, P, vmax):
    """C4: at level v the first bit with p(b) > v is 0, so the run ends."""
    G = len(rows)
    out = []
    for v in range(vmax + 1):
        stop = next((b for b, p in enumerate(P) if p is not None and p > v), None)
        if stop is None:
            break
        step = 1 << v
        gs = [g for g in range(step - 1, G, 2 * step)]   # v2(g+1) == v exactly
        zero = all(column(rows, stop)[g] == 0 for g in gs)
        out.append({"v": v, "stop_bit": stop, "zero_there": zero,
                    "samples": len(gs)})
    return out


def c5_table(P, TABLE, vmax):
    """C5: TABLE[v] == #{b : p(b) <= v}."""
    out = []
    for v in range(vmax + 1):
        if v not in TABLE:
            break
        pred = sum(1 for p in P if p is not None and p <= v)
        out.append({"v": v, "table": TABLE[v], "spectrum_count": pred,
                    "agree": pred == TABLE[v]})
    return out


# ------------------------------------------------------------------- driver
def run(G=1 << 16, nbits=48, kmax=12, vmax=12):
    from sitoolz import edgelaw
    rows = edge_rows(G)
    P = spectrum(rows, nbits, kmax)

    print(f"edge evolved for {G:,} generations; periods trusted to 2^{kmax}\n")

    print("C1  every bit is exactly periodic, period a power of two")
    info = c1_periodic(P)
    print(f"    {info['resolved_bits']} bits resolved, "
          f"largest period 2^{info['max_exp']}")
    print("    b   :", " ".join(f"{b:3}" for b, p in enumerate(P) if p is not None))
    print("    p(b):", " ".join(f"{p:3}" for p in P if p is not None))
    print()

    print("C2  period doubles exactly when the driver has odd weight")
    rows2 = c2_weight_rule(rows, P)
    bad = [r for r in rows2 if not r["agree"]]
    print(f"    {len(rows2)} bits tested, {len(bad)} disagreements")
    print()

    print("C3  bits with period 2^p are all 1 at g == -1 mod 2^p")
    rows3 = c3_in_phase(rows, P)
    bad = [r for r in rows3 if not r["always_one_at_-1_mod_2^p"]]
    print(f"    {len(rows3)} bits tested, {len(bad)} failures")
    print()

    print("C4  the first out-of-phase bit is 0 there, so the run stops")
    rows4 = c4_run_stops(rows, P, vmax)
    bad = [r for r in rows4 if not r["zero_there"]]
    print(f"    {len(rows4)} levels tested, {len(bad)} failures")
    print()

    print("C5  TABLE[v] == #{b : period(b) divides 2^v}")
    rows5 = c5_table(P, edgelaw.TABLE, vmax)
    print("     v   TABLE   spectrum   ")
    for r in rows5:
        print(f"    {r['v']:2}   {r['table']:5}   {r['spectrum_count']:8}   "
              f"{'ok' if r['agree'] else 'MISMATCH'}")
    bad = [r for r in rows5 if not r["agree"]]
    print(f"\n    {len(rows5)} levels tested, {len(bad)} mismatches")

    return {"C2": rows2, "C3": rows3, "C4": rows4, "C5": rows5}


if __name__ == "__main__":
    run()


# ---------------------------------------------------------------- deep mode
def bytes_columns(rows, nbits):
    """Columns as bytes objects -- period tests then run at C speed."""
    return [bytes((r >> b) & 1 for r in rows) for b in range(nbits)]


def period_exp_bytes(C, kmax):
    for k in range(kmax + 1):
        p = 1 << k
        if C[p:] == C[:-p]:
            return k
    return None


def deep(G=1 << 19, nbits=48, kmax=15, W=128):
    """C1-C5 at depth.  Slower to build, but resolves many more bits.

    Uses byte-per-generation columns so the period test is a single bytes
    comparison instead of a Python loop.  Defaults take a few minutes.
    """
    from sitoolz import edgelaw
    rows = edge_rows(G, W)
    cols = bytes_columns(rows, nbits)
    P = [period_exp_bytes(C, kmax) for C in cols]
    res = [(b, p) for b, p in enumerate(P) if p is not None]

    print(f"deep: {G:,} generations, {nbits} bits, periods trusted to 2^{kmax}\n")
    print(f"C1  {len(res)} bits resolved, largest period 2^{max(p for _, p in res)}")
    print("    p(b):", " ".join(str(p) for _, p in res))
    print()

    fails3 = 0
    for b, p in res:
        step = 1 << p
        if any(cols[b][g] != 1 for g in range(step - 1, G, step)):
            fails3 += 1
    print(f"C3  {len(res)} bits tested at g = -1 mod 2^p, {fails3} failures")

    levels = 0
    fails4 = 0
    for v in range(kmax):
        stop = next((b for b, p in enumerate(P) if p is not None and p > v), None)
        if stop is None:
            break
        levels += 1
        step = 1 << v
        if any(cols[stop][g] != 0 for g in range(step - 1, G, 2 * step)):
            fails4 += 1
    print(f"C4  {levels} levels tested, {fails4} failures")
    print()

    print("C5   v   TABLE   spectrum")
    n5 = bad5 = 0
    for v in sorted(edgelaw.TABLE):
        if v >= kmax:
            break
        pred = sum(1 for p in P if p is not None and p <= v)
        ok = pred == edgelaw.TABLE[v]
        n5 += 1
        bad5 += (not ok)
        print(f"    {v:2}   {edgelaw.TABLE[v]:5}   {pred:8}   {'ok' if ok else 'MISMATCH'}")
    print(f"\n    {n5} levels tested, {bad5} mismatches")
    return P
