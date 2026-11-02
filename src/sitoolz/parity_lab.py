"""parity_lab -- what the doubling indicator is actually made of.

Stdlib only.  Item 4 of TOMORROW.md asked whether delta can be computed
without paying 2^p generations per coordinate.  This file does not remove the
wall, but it localises it precisely, and proves two of the three terms away.

Write c_j(t) for coordinate j of the light-cone edge at generation t, and
d_j = c_{j-1} OR c_{j-2} for its driver.  Coordinate j is the running XOR-sum
of d_j, so its period doubles exactly when d_j has odd weight over one period:

    delta_j = parity of  sum_{t < T} d_j(t),   T = 2^max(p_{j-1}, p_{j-2}).

Since  w(A OR B) = w(A) + w(B) - w(A AND B),  that splits into three parity
terms.  Write sigma_j for the weight parity of coordinate j over one of *its
own* periods, and J_j for the joint term:

    delta_j = sigma_{j-1}  XOR  [p_{j-2} == p_{j-1}] sigma_{j-2}  XOR  J_j     (D)

    J_j = parity of  sum_{t < T} c_{j-1}(t) c_{j-2}(t)

The two singleton terms are free, by:

    ANTIPODAL LEMMA.  If coordinate j doubles (delta_j = 1) then
    c_j(t + T_j/2) = NOT c_j(t) for every t.

    Proof.  c_j(t + T_j/2) XOR c_j(t) is the XOR-sum of d_j over a window of
    length T_j/2, which is one full period of d_j; that sum is 1 for every t,
    because odd weight over one period is exactly what made j double.

    COROLLARY.  A doubling coordinate has weight exactly T_j/2 over its
    period, so sigma_j = 0 whenever p_j >= 2.

    FRESH RULE.  If coordinate j-1 doubled and p_{j-1} >= 2, then
    sigma_{j-1} = 0 and p_{j-2} < p_{j-1}, so (D) collapses to

        delta_j = J_j

    -- the doubling of a coordinate that follows a doubling is *entirely*
    the joint AND-parity.  Nothing else survives.

So the exponential wall is not in the OR at all.  It is one term, and only on
the coordinates that follow a doubling.

    python3 parity_lab.py
"""
from sitoolz import driver_lab as dl

EDGE_L, EDGE_N = 3, 1                 # a + b + ab  ==  a OR b


def edge_columns(nbits=26, steps=1 << 15, kmax=14):
    R = dl.rows(EDGE_L, EDGE_N, nbits, steps)
    p = dl.spectrum(R, nbits, kmax)
    cols = [bytes((r >> j) & 1 for r in R) for j in range(nbits)]
    return p, cols


def quantities(p, cols):
    """p, sigma, delta, J for every resolved coordinate."""
    res = [j for j, v in enumerate(p) if v is not None]
    sigma, delta, J, weight = {}, {}, {}, {}
    for j in res:
        T = 1 << p[j]
        weight[j] = sum(cols[j][:T])
        sigma[j] = weight[j] & 1
    for j in res:
        if j >= 1 and p[j - 1] is not None:
            delta[j] = 1 if p[j] == p[j - 1] + 1 else 0
        if j >= 2 and p[j - 1] is not None and p[j - 2] is not None:
            T = 1 << max(p[j - 1], p[j - 2])
            J[j] = sum(x & y for x, y in
                       zip(cols[j - 1][:T], cols[j - 2][:T])) & 1
    return res, sigma, delta, J, weight


def check_antipodal(p, cols, delta):
    """c_j(t + T/2) == NOT c_j(t) for every doubling coordinate, every t."""
    tested = bad = 0
    for j, d in delta.items():
        if not d or p[j] < 1:
            continue
        T = 1 << p[j]
        h = T >> 1
        tested += 1
        if any(cols[j][t] == cols[j][t + h] for t in range(h)):
            bad += 1
    return tested, bad


def check_halfweight(p, delta, weight):
    """Corollary: a doubling coordinate has weight exactly T/2."""
    tested = bad = 0
    for j, d in delta.items():
        if not d:
            continue
        tested += 1
        if weight[j] != (1 << p[j]) >> 1:
            bad += 1
    return tested, bad


def check_decomposition(p, sigma, delta, J):
    """Identity (D)."""
    tested = bad = 0
    for j in sorted(delta):
        if j not in J:
            continue
        tested += 1
        rhs = sigma[j - 1] ^ (sigma[j - 2] if p[j - 2] == p[j - 1] else 0) ^ J[j]
        if rhs != delta[j]:
            bad += 1
    return tested, bad


def check_fresh(p, delta, J):
    """FRESH RULE: delta_{j-1} = 1 and p_{j-1} >= 2  =>  delta_j == J_j."""
    tested = bad = 0
    for j in sorted(delta):
        if j not in J or (j - 1) not in delta:
            continue
        if delta[j - 1] == 1 and p[j - 1] >= 2:
            tested += 1
            if delta[j] != J[j]:
                bad += 1
    return tested, bad


def report(nbits=26, steps=1 << 15, kmax=14):
    p, cols = edge_columns(nbits, steps, kmax)
    res, sigma, delta, J, weight = quantities(p, cols)
    print(f"rule 30's edge, {len(res)} coordinates resolved, "
          f"largest period 2^{max(p[j] for j in res)}\n")
    print("  j   p   period   weight   w/T      sigma  delta  J")
    for j in res:
        T = 1 << p[j]
        print(f"{j:3} {p[j]:3} {T:8} {weight[j]:8}   {weight[j]/T:.4f}"
              f"{sigma[j]:7}{delta.get(j,'-'):>7}{J.get(j,'-'):>3}")
    print()
    for label, (n, bad) in (
            ("antipodal lemma  c(t+T/2) = NOT c(t) on doublers",
             check_antipodal(p, cols, delta)),
            ("corollary        weight == T/2 on doublers",
             check_halfweight(p, delta, weight)),
            ("identity (D)     delta = sigma + sigma' + J",
             check_decomposition(p, sigma, delta, J)),
            ("fresh rule       delta_j = J_j after a doubling",
             check_fresh(p, delta, J))):
        print(f"{label:<52} {n:3} tested, {bad} failures")
    frac = check_fresh(p, delta, J)[0] / max(1, len(J))
    print(f"\nthe joint term is the whole story on {frac:.0%} of coordinates;"
          f"\nthe rest carry a sigma from a plateau coordinate as well.")
    return {"p": p, "sigma": sigma, "delta": delta, "J": J}


if __name__ == "__main__":
    report()


# ---------------------------------------------------------------- dataframes
def parity_frame(nbits=26, steps=1 << 15, kmax=14):
    """The whole parity picture, one row per resolved coordinate.

    Columns: j, p, period, weight, ratio (weight/period), sigma, delta, J,
    doubled, half_weight (the antipodal corollary), fresh (whether the fresh
    rule applies here), and fresh_ok (delta == J where it does).
    """
    from sitoolz.dfkit import frame
    p, cols = edge_columns(nbits, steps, kmax)
    res, sigma, delta, J, weight = quantities(p, cols)
    rows = []
    for j in res:
        T = 1 << p[j]
        d = delta.get(j)
        fresh = (j - 1) in delta and delta.get(j - 1) == 1 and p[j - 1] >= 2 \
            and j in J
        rows.append({
            "j": j, "p": p[j], "period": T, "weight": weight[j],
            "ratio": round(weight[j] / T, 4), "sigma": sigma[j],
            "delta": d, "J": J.get(j),
            "doubled": None if d is None else bool(d),
            "half_weight": weight[j] == T >> 1,
            "fresh": fresh,
            "fresh_ok": (delta[j] == J[j]) if fresh else None,
        })
    return frame(rows, index="j")
