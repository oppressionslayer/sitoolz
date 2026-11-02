"""driver_lab -- the mixture question and the depth law, measured.

Stdlib only.  Answers items 1 and 2 of TOMORROW.md:

    1. Is "linear part tame + nonlinear part -> linear growth" generic?
       128 drivers = 8 linear parts x 16 nonlinear parts, four seeds each.
    2. The closed form for shift depth > 1, per residue class.

A *triangular XOR-driven system* is  x_j(t+1) = x_j(t) XOR G(x_{j-1}, x_{j-2},
x_{j-3}),  x_j = 0 for j < 0.  Written over GF(2), G is a polynomial in
a = x_{j-1}, b = x_{j-2}, c = x_{j-3}:

    G = (linear part L subset of {a,b,c}) + (nonlinear part N subset of {ab,ac,bc,abc})

Rule 30's light-cone edge is L = {a,b}, N = {ab}  --  a OR b.

    python3 driver_lab.py            # both experiments, ~2 minutes

Everything here is bit-parallel: coordinate j lives in bit j of one int, so
a = x << 1 reads the coordinate *below* j.  `x >> 1` would read the one above
and silently run a different, non-triangular system -- so `crossvalidate()`
checks every 2-input driver against the list-indexed evaluator in tri_xor.py
before any result below is trusted.
"""

# ----------------------------------------------------------------- drivers
LIN_NAMES = {0: "-", 1: "a", 2: "b", 3: "a+b", 4: "c", 5: "a+c", 6: "b+c", 7: "a+b+c"}
NON_NAMES = {0: "-", 1: "ab", 2: "ac", 4: "bc", 8: "abc"}


def non_name(N):
    if N == 0:
        return "-"
    return "+".join(NON_NAMES[bit] for bit in (1, 2, 4, 8) if N & bit)


def name(L, N):
    l, n = LIN_NAMES[L], non_name(N)
    if L and N:
        return f"{l}+{n}"
    return l if L else n


def step(x, L, N, m):
    """One generation, bit-parallel.  Bit j of the result is coordinate j."""
    a = (x << 1) & m
    b = (x << 2) & m
    c = (x << 3) & m
    g = 0
    if L & 1: g ^= a
    if L & 2: g ^= b
    if L & 4: g ^= c
    if N & 1: g ^= a & b
    if N & 2: g ^= a & c
    if N & 4: g ^= b & c
    if N & 8: g ^= a & b & c
    return (x ^ g) & m


def rows(L, N, nbits, steps, seed=1):
    m = (1 << nbits) - 1
    x = seed & m
    out = [x]
    for _ in range(steps - 1):
        x = step(x, L, N, m)
        out.append(x)
    return out


# --------------------------------------------------------------- spectrum
def spectrum(R, nbits, kmax):
    """p(j) for every coordinate at once: least k with column j of period 2^k.

    All coordinates share one pass per k.  `bad` accumulates, over t, the
    coordinates where row t and row t+2^k differ; a coordinate still clear
    at the end of the pass has period 2^k.
    """
    n = len(R)
    full = (1 << nbits) - 1
    p = [None] * nbits
    open_ = full
    for k in range(kmax + 1):
        P = 1 << k
        if P >= n:
            break
        bad = 0
        for t in range(n - P):
            bad |= R[t] ^ R[t + P]
            if bad & open_ == open_:
                break
        settled = open_ & ~bad
        for j in range(nbits):
            if (settled >> j) & 1:
                p[j] = k
        open_ &= bad
        if not open_:
            break
    return p


def envelope(p):
    """max p(0..j), the only classifier that survives non-monotone spectra."""
    out, best = [], -1
    for v in p:
        if v is None:
            out.append(None)
            continue
        best = max(best, v)
        out.append(best)
    return out


def profile(L, N, nbits=32, steps=1 << 13, kmax=12, seed=1):
    p = spectrum(rows(L, N, nbits, steps, seed), nbits, kmax)
    env = envelope(p)
    star = next((j for j, v in enumerate(p) if v is None), None)
    last = max((v for v in p if v is not None), default=0)
    return {"p": p, "env": env, "first_unresolved": star, "max_exp": last}


def classify(pr, nbits=32, kmax=12):
    """LOG (closed-form candidate), LINEAR (outruns the window), or FROZEN."""
    if pr["max_exp"] == 0 and pr["first_unresolved"] is None:
        return "FROZEN"
    if pr["first_unresolved"] is not None:
        return "LINEAR"
    # resolved everywhere: logarithmic iff the envelope stays near log2(nbits)
    return "LOG" if pr["max_exp"] <= (nbits - 1).bit_length() + 1 else "LINEAR"


# --------------------------------------------------------- cross-validation
def crossvalidate(nbits=24, steps=1 << 11, kmax=10):
    """Every 2-input driver, bit-parallel vs tri_xor's list-indexed evaluator."""
    from sitoolz import tri_xor as tx
    bad = []
    for code in range(16):
        # tri_xor's code is bit ((a<<1)|b) of `code`; translate to (L, N).
        g00 = code & 1
        if g00:                       # constant term: outside this parametrisation
            continue
        gb = (code >> 1) & 1          # G(0,1) -> coefficient of b
        ga = (code >> 2) & 1          # G(1,0) -> coefficient of a
        gab = (code >> 3) & 1         # G(1,1)
        L = ga | (gb << 1)
        N = (ga ^ gb ^ gab) & 1       # coefficient of ab in the ANF
        mine = spectrum(rows(L, N, nbits, steps), nbits, kmax)
        theirs = tx.spectrum(code, nbits, steps, kmax)
        if mine != theirs:
            bad.append((code, tx.NAMES[code], mine, theirs))
    return bad


# ------------------------------------------------- experiment 1: the mixture
DENSE = [0xA5A5A5A5, 0x6DB7C31F, 0x1F3E5C7B]


def mixture(nbits=32, steps=1 << 13, kmax=12, verbose=True):
    """8 linear parts x 16 nonlinear parts, single-cell seed + 3 dense seeds.

    The question: does adding any nonlinear part to a tame (depth-1 linear)
    driver always destroy the closed form?
    """
    seeds = [1] + DENSE
    table = {}
    for L in range(8):
        for N in range(16):
            verdicts = []
            for s in seeds:
                pr = profile(L, N, nbits, steps, kmax, seed=s)
                verdicts.append(classify(pr, nbits, kmax))
            table[(L, N)] = verdicts
    if verbose:
        print("mixture: rows = linear part, cols = nonlinear part")
        print("verdict per driver over [single-cell, dense1, dense2, dense3]")
        print("L = log-growth (closed form)   X = linear growth   . = frozen\n")
        code = {"LOG": "L", "LINEAR": "X", "FROZEN": "."}
        head = "".join(f"{non_name(N):>10}" for N in range(16))
        print(f"{'':>8}{head}")
        for L in range(8):
            cells = []
            for N in range(16):
                cells.append("".join(code[v] for v in table[(L, N)]))
            print(f"{LIN_NAMES[L]:>8}" + "".join(f"{c:>10}" for c in cells))
    return table


def mixture_report(table):
    """The claim, and every driver that violates it."""
    depth1 = [L for L in range(8) if L & 1]        # linear part contains a
    tame_lin, exceptions, always_tame = [], [], []
    for L in range(8):
        v = table[(L, 0)]
        tame_lin.append((LIN_NAMES[L], v))
    for L in depth1:
        for N in range(1, 16):
            v = table[(L, N)]
            if all(x != "LINEAR" for x in v):
                exceptions.append((name(L, N), v))
    for L in range(8):
        for N in range(16):
            if all(x == "LOG" for x in table[(L, N)]):
                always_tame.append(name(L, N))
    print("\npure linear parts (no nonlinear term) -- should all be LOG:")
    for nm, v in tame_lin:
        print(f"   {nm:>8}  {v}")
    print(f"\ndepth-1 linear part + any nonzero nonlinear part: {len(depth1)*15} drivers")
    if exceptions:
        print(f"   {len(exceptions)} never reach linear growth under any seed:")
        for nm, v in exceptions:
            print(f"   {nm:>12}  {v}")
    else:
        print("   0 exceptions: every one goes linear under at least one seed.")
    print(f"\ndrivers that are LOG under all four seeds: {len(always_tame)}")
    print("   " + ", ".join(always_tame))
    return {"exceptions": exceptions, "always_tame": always_tame}


# ------------------------------------------- experiment 2: the depth law
def shift_rows(S, nbits, steps, seed=1):
    """x' = x XOR sum_{s in S} (x << s) -- a purely linear driver, any depth."""
    m = (1 << nbits) - 1
    x = seed & m
    out = [x]
    for _ in range(steps - 1):
        g = 0
        for s in S:
            g ^= (x << s) & m
        x = (x ^ g) & m
        out.append(x)
    return out


def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def clog2(n):
    """ceil(log2(n)) for n >= 1."""
    return 0 if n <= 1 else (n - 1).bit_length()


def predict(S, j, seeded_classes):
    """Predicted p(j) for the linear shift-set system.

    g = gcd(S) splits the coordinates into g interleaved subsystems on the
    residue classes mod g; class r carries subsystem index m = (j-r)/g and is
    a copy of the S/g system, whose minimum shift is d = min(S)/g.  The
    Frobenius argument (1+Lambda)^(2^k) = 1 + Lambda^(2^k) then gives

        p = ceil(log2((m+1)/d))

    on a seeded class, and 0 on a class that is identically zero.
    """
    g = 0
    for s in S:
        g = gcd(g, s)
    r = j % g
    if r not in seeded_classes:
        return 0
    m = j // g
    d = min(S) // g
    return clog2((m + 1 + d - 1) // d) if m + 1 > d else 0


def depth_law(nbits=40, steps=1 << 13, kmax=12, verbose=True):
    """Every shift set S subset of {1..5}: measured spectrum vs prediction."""
    sets = []
    for mask in range(1, 32):
        S = tuple(s for s in range(1, 6) if mask >> (s - 1) & 1)
        sets.append(S)
    results = []
    for S in sets:
        g = 0
        for s in S:
            g = gcd(g, s)
        for label, seed, classes in (("single", 1, {0}),
                                     ("dense", 0x6DB7C31F, set(range(g)))):
            p = spectrum(shift_rows(S, nbits, steps, seed), nbits, kmax)
            pred = [predict(S, j, classes) for j in range(nbits)]
            miss = [(j, p[j], pred[j]) for j in range(nbits)
                    if p[j] is not None and p[j] != pred[j]]
            results.append({"S": S, "seed": label, "gcd": g, "dmin": min(S),
                            "resolved": sum(1 for v in p if v is not None),
                            "mismatches": miss, "p": p, "pred": pred})
    if verbose:
        print(f"{'shift set':>16} {'seed':>7} {'gcd':>4} {'d':>3} "
              f"{'resolved':>9} {'mismatches':>11}")
        for r in results:
            print(f"{str(r['S']):>16} {r['seed']:>7} {r['gcd']:>4} "
                  f"{r['dmin']:>3} {r['resolved']:>9} {len(r['mismatches']):>11}"
              + ("" if not r["mismatches"] else
                 "   first: j=%d got %s want %s" % r["mismatches"][0]))
        bad = [r for r in results if r["mismatches"]]
        tot = sum(r["resolved"] for r in results)
        print(f"\n{len(results)} systems, {tot} resolved coordinates, "
              f"{len(bad)} systems with any mismatch")
    return results



# ------------------------------------------- experiment 1b: the exceptions
SEEDS8 = [1, 0xA5A5A5A5, 0x6DB7C31F, 0x1F3E5C7B, 0xFFFFFFFF, 0x2, 0x80000001,
          0x13579BDF]


def truth(L, N):
    """The driver as a truth table: bit (a<<2|b<<1|c) of the returned int."""
    tt = 0
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                g = 0
                if L & 1: g ^= a
                if L & 2: g ^= b
                if L & 4: g ^= c
                if N & 1: g ^= a & b
                if N & 2: g ^= a & c
                if N & 4: g ^= b & c
                if N & 8: g ^= a & b & c
                tt |= g << ((a << 2) | (b << 1) | c)
    return tt


def is_subcube(tt):
    """Is the driver a conjunction of literals -- the indicator of a subcube?"""
    pts = [i for i in range(8) if tt >> i & 1]
    if not pts:
        return False
    fixed = [k for k in range(3)
             if len({(q >> k) & 1 for q in pts}) == 1]
    return len(pts) == 1 << (3 - len(fixed))


def fires_on_ones(tt):
    """G(1,1,1) -- does the driver fire when every input coordinate is 1?"""
    return (tt >> 7) & 1


def mixture_deep(nbits=40, steps=1 << 14, kmax=13, seeds=None, verbose=True):
    """All 128 drivers x 8 seeds, and the structural test of the tame set.

    HYPOTHESIS.  A driver is tame (logarithmic period growth under every seed)
    exactly when it is linear, or it is the indicator of a subcube that does
    not contain the all-ones point.  Everything else reaches linear growth
    under at least one seed.
    """
    seeds = seeds or SEEDS8
    rec = []
    for L in range(8):
        for N in range(16):
            v = [classify(profile(L, N, nbits, steps, kmax, seed=s), nbits, kmax)
                 for s in seeds]
            tt = truth(L, N)
            tame = all(x != "LINEAR" for x in v)
            pred = (tt == 0) or (N == 0) or (is_subcube(tt) and not fires_on_ones(tt))
            rec.append({"L": L, "N": N, "name": name(L, N), "tt": tt,
                        "verdicts": v, "tame": tame, "predicted_tame": pred,
                        "subcube": is_subcube(tt), "ones": fires_on_ones(tt)})
    agree = [r for r in rec if r["tame"] == r["predicted_tame"]]
    if verbose:
        tame = [r for r in rec if r["tame"]]
        print(f"{len(rec)} drivers x {len(seeds)} seeds, "
              f"{nbits} coordinates, periods resolved to 2^{kmax}\n")
        print(f"tame under every seed: {len(tame)}")
        print(f"{'driver':>16} {'subcube':>8} {'G(1,1,1)':>9}   verdicts")
        for r in tame:
            print(f"{r['name']:>16} {str(r['subcube']):>8} {r['ones']:>9}   "
                  + "".join({'LOG': 'L', 'LINEAR': 'X', 'FROZEN': '.'}[x]
                            for x in r["verdicts"]))
        print(f"\nhypothesis: tame  <=>  linear, or a subcube avoiding (1,1,1)")
        print(f"   {len(agree)}/{len(rec)} drivers agree")
        for r in rec:
            if r["tame"] != r["predicted_tame"]:
                print(f"   COUNTEREXAMPLE {r['name']:>16} measured "
                      f"{'tame' if r['tame'] else 'wild'}, predicted "
                      f"{'tame' if r['predicted_tame'] else 'wild'}   {r['verdicts']}")
    return rec


# ------------------------------- experiment 2b: the spectrum without simulation
def pmul(a, b, N):
    """Multiply two GF(2) polynomials held as bitmasks, truncated at z^N."""
    r = 0
    m = (1 << (N + 1)) - 1
    while b:
        low = b & -b
        r ^= (a << (low.bit_length() - 1))
        b ^= low
    return r & m


def spectrum_linear_exact(S, jmax=63):
    """p(j) for x' = x XOR sum_{s in S} (x << s), exactly, with no simulation.

    Write P(z) = sum_{s in S} z^s, so one generation is multiplication by
    1 + P and the orbit from a single cell is x(t) = (1+P)^t.  Frobenius gives
    (1+P)^(2^k) = 1 + P(z^(2^k)) =: 1 + Q_k, hence

        x_j(t + 2^k) + x_j(t) = [z^j] ( (1+P)^t * Q_k ).

    By Lucas, (1+P)^t factors over the binary digits of t as a product of
    distinct (1 + Q_i), so the whole orbit is indexed by finite subsets T of
    the indices i with min(S)*2^i <= jmax -- about log2(jmax) of them.  So

        p(j) = least k such that [z^j] (F_T * Q_k) = 0 for every subset T.

    That is a polynomial-time computation: 2^|I| products of truncated
    polynomials, no generations evolved at all.  It is exact for every shift
    set, including the ones where the ceil(log2) closed form is only a bound.
    """
    P = 0
    for s in S:
        P |= 1 << s
    d = min(S)
    I = [i for i in range(jmax.bit_length() + 2) if d << i <= jmax]
    K = len(I) + 1

    Q = []
    for k in range(K + 1):
        q = 0
        for s in S:
            if s << k <= jmax:
                q |= 1 << (s << k)
        Q.append(q)

    A = [0] * (K + 1)                       # A[k] = coordinates that still move

    def walk(idx, F):
        for k in range(K + 1):
            if Q[k]:
                A[k] |= pmul(F, Q[k], jmax)
        if idx == len(I):
            return
        walk(idx + 1, F)
        walk(idx + 1, pmul(F, 1 ^ Q[I[idx]], jmax))

    walk(0, 1)

    out = []
    for j in range(jmax + 1):
        p = None
        for k in range(K + 1):
            if not (A[k] >> j) & 1:
                p = k
                break
        out.append(p)
    return out


def check_exact(nbits=40, steps=1 << 13, kmax=12, smax=5, verbose=True):
    """The simulation-free spectrum against the measured one, every shift set."""
    rows_out = []
    for mask in range(1, 1 << smax):
        S = tuple(s for s in range(1, smax + 1) if mask >> (s - 1) & 1)
        meas = spectrum(shift_rows(S, nbits, steps), nbits, kmax)
        exact = spectrum_linear_exact(S, nbits - 1)
        miss = [(j, meas[j], exact[j]) for j in range(nbits)
                if meas[j] is not None and meas[j] != exact[j]]
        g = 0
        for s in S:
            g = gcd(g, s)
        closed = [predict(S, j, {0}) for j in range(nbits)]
        cmiss = [j for j in range(nbits)
                 if meas[j] is not None and meas[j] != closed[j]]
        rows_out.append({"S": S, "gcd": g, "min": min(S),
                         "exact_mismatches": len(miss),
                         "closed_form_mismatches": len(cmiss)})
    if verbose:
        print(f"{'shift set':>16} {'gcd':>4} {'min':>4} "
              f"{'ceil-log2 form':>15} {'exact rule':>11}")
        for r in rows_out:
            print(f"{str(r['S']):>16} {r['gcd']:>4} {r['min']:>4} "
                  f"{r['closed_form_mismatches']:>15} {r['exact_mismatches']:>11}")
        tot = len(rows_out)
        bad_c = sum(1 for r in rows_out if r["closed_form_mismatches"])
        bad_e = sum(1 for r in rows_out if r["exact_mismatches"])
        print(f"\n{tot} shift sets: ceil-log2 form fails on {bad_c}, "
              f"exact rule fails on {bad_e}")
        print("the closed form is exact exactly when min(S) == gcd(S):",
              all((r["closed_form_mismatches"] == 0) == (r["min"] == r["gcd"])
                  for r in rows_out))
    return rows_out


if __name__ == "__main__":
    print("cross-validation against tri_xor's list-indexed evaluator")
    bad = crossvalidate()
    print(f"   {len(bad)} disagreements over the 2-input drivers\n")

    t = mixture()
    mixture_report(t)
    print("\n" + "=" * 72)
    print("depth law: p(j) = ceil(log2((m+1)/d)),  m = j // gcd(S), d = min(S)/gcd(S)\n")
    depth_law()
    print("\n" + "=" * 72)
    print("the same spectra with nothing evolved -- and where the closed form fails\n")
    check_exact()


# ---------------------------------------------------------------- dataframes
def driver_frame(nbits=32, steps=1 << 13, kmax=12, seeds=None, verbose=False):
    """Every driver as a row: structure, prediction, and measured verdict.

    Columns: L, N, name, truth (table as an int), linear, subcube, ones
    (G(1,1,1)), predicted_tame, tame (measured), n_linear/n_log/n_frozen
    (verdicts across the seeds), and the per-seed verdict string.

    Defaults are the snappy ones (~5 s).  The published 128/128 figure used
    nbits=40, steps=2^14, kmax=13 and the 8 seeds in SEEDS8.
    """
    from sitoolz.dfkit import frame
    seeds = seeds or SEEDS8[:4]
    rows = []
    for L in range(8):
        for N in range(16):
            v = [classify(profile(L, N, nbits, steps, kmax, seed=s), nbits, kmax)
                 for s in seeds]
            tt = truth(L, N)
            sub = is_subcube(tt)
            ones = bool(fires_on_ones(tt))
            rows.append({
                "L": L, "N": N, "name": name(L, N), "truth": tt,
                "linear": N == 0, "subcube": sub, "ones": ones,
                "predicted_tame": (tt == 0) or (N == 0) or (sub and not ones),
                "tame": all(x != "LINEAR" for x in v),
                "n_log": v.count("LOG"), "n_linear": v.count("LINEAR"),
                "n_frozen": v.count("FROZEN"),
                "verdicts": "".join({"LOG": "L", "LINEAR": "X",
                                     "FROZEN": "."}[x] for x in v),
            })
    return frame(rows)


def spectrum_frame(L=3, N=1, seed=1, nbits=32, steps=1 << 13, kmax=12):
    """One driver's spectrum, coordinate by coordinate.

    Columns: j, p (period exponent), period, envelope (max p up to j),
    resolved, frozen.
    """
    from sitoolz.dfkit import frame
    pr = profile(L, N, nbits, steps, kmax, seed=seed)
    rows = []
    for j, (p, e) in enumerate(zip(pr["p"], pr["env"])):
        rows.append({"j": j, "p": p,
                     "period": None if p is None else 1 << p,
                     "envelope": e, "resolved": p is not None,
                     "frozen": p == 0})
    return frame(rows, index="j")


def shiftset_frame(smax=5, nbits=40, steps=1 << 13, kmax=12):
    """Every linear shift set: gcd, min, and how each rule scores on it.

    Columns: S, gcd, min, terms (XORs per generation), closed_form_ok
    (the ceil-log2 form), exact_ok (the algebraic rule), mismatches for each.
    """
    from sitoolz.dfkit import frame
    rows = []
    for mask in range(1, 1 << smax):
        S = tuple(s for s in range(1, smax + 1) if mask >> (s - 1) & 1)
        g = 0
        for s in S:
            g = gcd(g, s)
        meas = spectrum(shift_rows(S, nbits, steps), nbits, kmax)
        exact = spectrum_linear_exact(S, nbits - 1)
        closed = [predict(S, j, {0}) for j in range(nbits)]
        em = sum(1 for j in range(nbits)
                 if meas[j] is not None and meas[j] != exact[j])
        cm = sum(1 for j in range(nbits)
                 if meas[j] is not None and meas[j] != closed[j])
        rows.append({"S": "(" + ",".join(map(str, S)) + ")", "gcd": g,
                     "min": min(S), "terms": len(S),
                     "min_eq_gcd": min(S) == g,
                     "closed_form_ok": cm == 0, "closed_form_miss": cm,
                     "exact_ok": em == 0, "exact_miss": em})
    return frame(rows)
