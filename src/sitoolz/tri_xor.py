"""tri_xor -- triangular XOR-driven systems over GF(2), standalone.

One file, stdlib only, no imports from anywhere else.  Everything in
AN_INVITATION.md can be reproduced from here.

A *triangular XOR-driven system* is

    x_j(t+1) = x_j(t) XOR G_j( x_{j-1}(t), ..., x_{j-d}(t) )      x_j = 0 for j < 0

for arbitrary Boolean G_j.  Bit j reads only bits below it, so the system is a
T-function in the sense of Klimov-Shamir (2002) -- equivalently a 1-Lipschitz
map on the 2-adic integers.

    THEOREM.  Every coordinate is exactly periodic in t, with period a power
    of two.

Rule 30's light-cone edge is the instance G = OR of the two below.

    from sitoolz import tri_xor as tx

    tx.spectrum(tx.OR)             # period exponents p(j)
    tx.delta(tx.OR)                # doubling indicator
    tx.edge_run(255)               # rule 30's edge run at generation 255 -> 24
    tx.survey()                    # all sixteen 2-input drivers
    tx.check()                     # the theorem's claims, verified
"""

# ---------------------------------------------------------------- drivers
# A 2-input driver is an int 0..15; bit ((a<<1)|b) of it is G(a, b).
ZERO, NOR, B_AND_NOTA, NOT_A = 0, 1, 2, 3
A_AND_NOTB, NOT_B, XOR, NAND = 4, 5, 6, 7
AND, XNOR, B, NOTA_OR_B     = 8, 9, 10, 11
A, A_OR_NOTB, OR, ONE       = 12, 13, 14, 15

NAMES = {0:"0", 1:"NOR", 2:"b AND NOT a", 3:"NOT a", 4:"a AND NOT b", 5:"NOT b",
         6:"XOR", 7:"NAND", 8:"AND", 9:"XNOR", 10:"b", 11:"(NOT a) OR b",
         12:"a", 13:"a OR (NOT b)", 14:"OR", 15:"1"}


def evolve(driver=OR, nbits=40, steps=1 << 13, seed=None):
    """Run the system and return one bytes column per coordinate.

    `driver` is an int 0..15 (a 2-input G, the same for every j).  `seed`
    defaults to a single 1 in coordinate 0.
    """
    x = bytearray(nbits)
    if seed is None:
        x[0] = 1
    else:
        x[:len(seed)] = bytes(seed)
    cols = [bytearray(steps) for _ in range(nbits)]
    for t in range(steps):
        for j in range(nbits):
            cols[j][t] = x[j]
        nx = bytearray(x)
        for j in range(nbits):
            a = x[j-1] if j >= 1 else 0
            b = x[j-2] if j >= 2 else 0
            nx[j] = x[j] ^ ((driver >> ((a << 1) | b)) & 1)
        x = nx
    return [bytes(c) for c in cols]


def period_exp(column, kmax=12):
    """Least k with `column` of period 2^k, or None if none <= 2^kmax."""
    for k in range(kmax + 1):
        p = 1 << k
        if column[p:] == column[:-p]:
            return k
    return None


def spectrum(driver=OR, nbits=40, steps=1 << 13, kmax=12):
    """p(j) for each coordinate.  None where the run is too short to resolve."""
    return [period_exp(c, kmax) for c in evolve(driver, nbits, steps)]


def delta(driver=OR, **kw):
    """The doubling indicator: 1 where the period doubled, else 0.

    This is the only place the choice of driver shows up.  The theorem fixes
    the periods as powers of two for every driver; delta is what differs.
    """
    P = spectrum(driver, **kw)
    out = []
    for j in range(1, len(P)):
        if P[j] is None or P[j-1] is None:
            out.append(None)
        else:
            out.append(1 if P[j] == P[j-1] + 1 else 0)
    return out


def plateaus(driver=OR, **kw):
    """Run lengths of coordinates sharing a period -- the gaps in delta, +1."""
    P = [p for p in spectrum(driver, **kw) if p is not None]
    out, run = [], 1
    for a, b in zip(P, P[1:]):
        if b == a:
            run += 1
        else:
            out.append(run); run = 1
    out.append(run)
    return out


# ------------------------------------------------------ rule 30's edge law
def v2(n):
    """2-adic valuation: how many times 2 divides n."""
    return (n & -n).bit_length() - 1


def edge_measure(g, width=192):
    """Rule 30's edge run at generation g, by evolving the edge alone.

    E' = E ^ ((E << 1) | (E << 2)).  Bit j reads only j-1 and j-2, so the
    window is exact for every bit it keeps, at any generation.
    """
    m = (1 << width) - 1
    E = 1
    for _ in range(g):
        E = (E ^ ((E << 1) | (E << 2))) & m
    r = 0
    while (E >> r) & 1:
        r += 1
    return r


def edge_table(vmax=20, **kw):
    """run(v) = #{ j : period(x_j) divides 2^v }, from the OR spectrum."""
    P = spectrum(OR, **kw)
    return [sum(1 for p in P if p is not None and p <= v) for v in range(vmax + 1)]


def edge_run(g, table=None):
    """The edge law: a valuation and a table index, no simulation."""
    if table is None:
        table = edge_table()
    return table[v2(g + 1)]


# ------------------------------------------------------------- the claims
def check(nbits=34, steps=1 << 14, kmax=13, verbose=True):
    """Verify the theorem and the two phase lemmas on rule 30's edge."""
    cols = evolve(OR, nbits, steps)
    P = [period_exp(c, kmax) for c in cols]
    res = [(j, p) for j, p in enumerate(P) if p is not None]
    out = {}

    out["C1_power_of_two"] = len(res)

    bad = sum(1 for j, p in res
              if any(cols[j][t] != 1 for t in range((1 << p) - 1, steps, 1 << p)))
    out["C3_in_phase_failures"] = bad

    bad = 0; lv = 0
    for v in range(kmax):
        j = next((i for i, p in enumerate(P) if p is not None and p > v), None)
        if j is None:
            break
        lv += 1
        step = 1 << v
        if any(cols[j][t] != 0 for t in range(step - 1, steps, 2 * step)):
            bad += 1
    out["C4_levels"] = lv
    out["C4_failures"] = bad

    tbl = [sum(1 for p in P if p is not None and p <= v) for v in range(lv)]
    bad = sum(1 for g in range(1, 4096)
              if v2(g + 1) < len(tbl) and edge_measure(g, 192) != tbl[v2(g + 1)])
    out["edge_law_mismatches_to_4096"] = bad

    if verbose:
        print(f"C1  {out['C1_power_of_two']} coordinates resolved, all periods powers of two")
        print(f"C3  in-phase bits show 1        -- {out['C3_in_phase_failures']} failures")
        print(f"C4  out-of-phase bits show 0    -- {out['C4_failures']} failures "
              f"over {out['C4_levels']} levels")
        print(f"law run(g) == table[v2(g+1)]    -- "
              f"{out['edge_law_mismatches_to_4096']} mismatches to generation 4096")
    return out


def survey(nbits=40, steps=1 << 12, kmax=11):
    """Every 2-input driver: how much of the spectrum resolves, and delta."""
    rows = []
    for g in range(16):
        P = spectrum(g, nbits, steps, kmax)
        res = [p for p in P if p is not None]
        d = [v for v in delta(g, nbits=nbits, steps=steps, kmax=kmax) if v is not None]
        rows.append({"driver": NAMES[g], "code": g,
                     "resolved": f"{len(res)}/{nbits}",
                     "max_period": f"2^{max(res)}" if res else "-",
                     "delta_density": round(sum(d)/len(d), 3) if d else None,
                     "opaque": len(res) < nbits})
    return rows


def show_survey(**kw):
    print(f"{'driver':>14}  {'resolved':>9}  {'max period':>10}  {'density':>7}  verdict")
    for r in survey(**kw):
        v = "opaque" if r["opaque"] else "explicit"
        print(f"{r['driver']:>14}  {r['resolved']:>9}  {r['max_period']:>10}  "
              f"{str(r['delta_density']):>7}  {v}")


if __name__ == "__main__":
    print("theorem and phase lemmas, on rule 30's edge:")
    check()
    print()
    print("driver space:")
    show_survey()
    print()
    print("edge law spot checks (all four share v2 = 8, so all give 24):")
    for g in (255, 767, 1279, 1791):
        print(f"  g={g:>5}  v2(g+1)={v2(g+1):2}  run={edge_run(g)}  "
              f"measured={edge_measure(g)}")
