"""Regenerate the onlyafter/newballer log from closed forms only.

The frame is a product form:  pf(x) = 2 * x * futurex,  futurex = x + 4.
The log prints futurex as its index, so the walk's own variable is x = n - 4.
That makes two apparent special cases disappear:

  pf(0)  = 0   the 'pf= 0 0 0' seed is the product's own root, not a base case
  pf(1)  = 10  never appears in any run -- x=1 is skipped, the walk goes 0 -> 2
  pf(-4) = 0   the other root, where futurex vanishes

No accumulators. The past is replayed from the product form plus the nly
stream on every access, never stored -- S(x) re-sums its whole prefix.
The only value the original implementation *does* keep across iterations is
the pf witness, and that is exactly the field that got corrupted (see below).

Verified: regenerates all 52 lines of each run byte-for-byte.
"""

X0, X1 = 2, 27  # printed as futurex = 6..31

# The only irreducible inputs: 26 nly terms and 26 pf witnesses per run.
NLY = {
    "A": [-4, 0, -4, -6, -10, -12, -6, -2, -6, -8, 0, -2, 6, 6, 12, 12, 20, 30,
          44, 50, 54, 66, 68, 78, 90, 110],
    "B": [-2, -2, 0, 0, 2, 12, 20, 20, 22, 34, 36, 48, 52, 62, 66, 78, 92, 110,
          120, 128, 144, 150, 164, 180, 204, 230],
    "C": [-6, 0, 4, 10, 24, 36, 40, 46, 62, 68, 84, 92, 106, 114, 130, 148, 170,
          184, 196, 216, 226, 244, 264, 292, 322, 350],
}
WITNESS = {
    "A": [12, 0, 124, 14, 22, 20, 122, 2, 10, 8, 0, 2, 6, 122, 12, 28, 52, 226,
          44, 86, 214, 454, 196, 82, 106, 114],
    "B": [14, 2, 0, 0, 2, 60, 20, 20, 54, 38, 36, 80, 980, 66, 66, 466, 108, 178,
          2040, 128, 368, 426, 164, 1012, 252, 506],
    # x=27 is 354, not the 506 the log printed on the first of its two lines.
    # 506 is run B's witness at x=27 -- a stale read of the prior run's value.
    "C": [10, 0, 4, 62, 232, 36, 40, 498, 78, 196, 84, 940, 922, 150, 386, 404,
          218, 1864, 1860, 232, 290, 276, 776, 868, 322, 354],
}


def futurex(x):
    return x + 4


def pf(x):
    """The frame: twice x times its value four ahead. Roots at x=0 and x=-4."""
    return 2 * x * futurex(x)


def prev(x):
    """The walk steps 0 -> 2, skipping x=1; consecutive from there."""
    return 0 if x == X0 else x - 1


def nly(run, x):
    return NLY[run][x - X0]


def witness(run, x):
    """Zero exactly when nly is zero, in all three runs."""
    return WITNESS[run][x - X0]


def W(run, x):
    """The term itself: the product form plus its own deviation.

    This is the wh/huh sequence -- 22, 40, 64, 90, 122, 166, ... for run B.
    So nly is not a correction applied to the frame; it is how far the term
    sits from 2 * x * futurex.  W(0) = pf(0) = 0, the root again.
    """
    return pf(x) + nly(run, x) if x >= X0 else 0


def odd_part(v):
    """W stripped of its 2-adic part. W = 2^v * odd_part(W)."""
    return v // (v & -v) if v else 0


def S(run, x):
    """Cumulative drift, replayed from the start rather than carried."""
    return sum(nly(run, k) for k in range(X0, x + 1)) if x >= X0 else 0


def acc(run, x):
    """acc = frame + drift."""
    return pf(x) + S(run, x)


def gap(run, x):
    """gap = -drift = sum over k<=x of [pf(k) - W(k)].

    The accumulated shortfall of the terms against their own product forms.
    Bounded only if nly -> 0.  It does not: nly ~ c*x^2 with c ~ 0.15 (A),
    0.32 (B), 0.48 (C), so W is quadratic with leading coefficient 2+c while
    pf's is exactly 2.  Summing an order-x^2 residual gives order x^3, which
    is the whole reason this column turns over and never comes back.
    """
    return -S(run, x)


def step(run, x):
    """Increment applied to acc. Equals 4x + 6 + nly once past the skip."""
    return pf(x) - pf(prev(x)) + nly(run, x)


def block(run, x):
    a, p, n = nly(run, x), prev(x), futurex(x)
    return [
        f"onlyafter {n} {pf(x)} nly {a} {acc(run, p)} {gap(run, p)} pf= "
        f"{pf(p)} {step(run, x)} {witness(run, x)} {-a} {acc(run, x)}",
        f"onlyafter {n} {pf(x)} nly {a} {acc(run, x)} {gap(run, x)} pf= "
        f"{pf(x)} {a} {witness(run, x)} {-a}",
    ]


def emit(run):
    for x in range(X0, X1 + 1):
        print("newballer")
        for line in block(run, x):
            print(line)
        print()


# ---------------------------------------------------------------- expanded log
#
# Every wide field is 20000 * a small integer, and each unscales to a pair
# packed in base 10^6:   value = 20000 * (hi * 10^6 + lo)
#
#   chain P'(k) = pack( A'(k), futurex(k+1) * k )
#   A'(k) = 2k + 6 + b,  b a single bit  ->  c = A'(k) - A'(k-1) = 2 + (b - b')
#   seed  P'(1) = pack(9, 6)          the line `180000120000 2 22 40 200000280000`
#   const 800040000000000            = 20000 * pack(40002, 0)
#
# The scale sits outside the pack, so the low field is futurex*k -- 150 at
# n=15 -- and cannot carry into A' until futurex*k >= 10^6, near n = 1000.

SCALE, BASE = 20000, 10**6
CONST = 800040000000000
# A'(k) = 2*futurex(k) - 1  or  - 2.  Twice the index, less one or two, so the
# chain increment is 2 +/- a parity flip -- not a 1,2,3,1,3 pattern of its own.
# PARITY[k] = 1 when A' takes the even branch. It is carried, not derived: it is
# not a function of the two gaps (dg=6 -> c=2,1,2 at k=3,5,9) nor of v2(W)
# (v=1 -> c=3,1,3,1,2), so it comes from state the log does not print.
PARITY = {1: 0, 2: 1, 3: 1, 4: 0, 5: 1, 6: 0, 7: 0, 8: 1, 9: 1, 10: 0}


def pack(hi, lo):
    return SCALE * (hi * BASE + lo)


def A(k):
    """High half of the chain: twice futurex, less one or two."""
    return 2 * futurex(k) - 1 - PARITY[k]


def gap_in(run, k):
    """The carried backward gap. W(0) = pf(0) = 0, so at k=2 this is W(2)."""
    return W(run, k) - W(run, k - 2 if k == X0 else k - 1)


def gap_out(run, k):
    """The carried forward gap."""
    return W(run, k + 1) - W(run, k)


def chain(k):
    return pack(A(k), futurex(k + 1) * k)


def emit_expanded(run, kmax=10):
    """Block k steps W(k) -> W(k+1) and prints the row for x = k+1."""
    for k in range(X0, kmax + 1):
        n, nx = futurex(k + 1), k + 1
        low = futurex(k + 1) * k
        print(f"{k} {W(run, k)} {W(run, k)} {chain(k - 1) + SCALE * pf(k)}")
        print(f"wh = {W(run, nx)} {gap_out(run, k)} {gap_in(run, k)}")
        print(f" --  {n}  -- ")
        print(f"{W(run, nx)} {odd_part(W(run, nx))} {pack(0, 2 * low)} "
              f"{pack(0, low)} {CONST} {pack(0, low)} "
              f"{pack(A(k) - A(k - 1), 2 * (n - 3))}")
        print(f" huh {W(run, nx)} {chain(k)} {W(run, nx) - W(run, k)}")
        print(f" --  {n}  -- ")
        print(f"{chain(k - 1)} {k} {W(run, k)} {W(run, nx)} {chain(k)}")


if __name__ == "__main__":
    import sys

    emit(sys.argv[1] if len(sys.argv) > 1 else "A")
