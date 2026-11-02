"""A variant of the sequence whose frame actually tracks its terms.

Same skeleton as the original:

    pf(x) = 2 * x * futurex,  futurex = x + 4
    W(x)  = pf(x) + nly
    gap   = sum over k<=x of [pf(k) - W(k)] = -sum(nly)

One rule change. In the original, nly grew like c*x^2 -- the same order as
pf itself -- so summing it gave an order-x^3 gap that could never come back.
Here nly is chosen as the SMALLEST step off the frame that lands on a W whose
odd part is prime. That keeps |nly| tiny, so the gap is a bounded walk instead
of a cubic runaway.

On a tie (+d and -d both work) take the sign that pulls the gap toward zero.
That feedback is what turns a slow drift into a bounded oscillation.

    unsteered:  gap in    0 .. 30   over x = 2..31
    steered:    gap in   -4 .. 14   over x = 2..59
    original:   gap = -586 / -1970 / -3422 at x = 27 (runs A / B / C)

And the odd part is prime for every term by construction -- what was a
one-directional hint in the original ("A' even => odd part prime") is here
an unconditional property of the sequence.
"""

from sympy import isprime

SEARCH = 400


def futurex(x):
    return x + 4


def pf(x):
    return 2 * x * futurex(x)


def odd_part(v):
    return v // (v & -v) if v else 0


def step(x, gap):
    """Nearest even W to pf(x) with a prime odd part, steering toward gap=0."""
    for d in range(0, SEARCH, 2):
        signs = [0] if d == 0 else ([-1, 1] if gap <= 0 else [1, -1])
        for s in signs:
            W = pf(x) + s * d
            if W > 0 and isprime(odd_part(W)):
                return W, s * d
    raise ValueError(f"no term within {SEARCH} of pf({x})")


def run(lo=2, hi=60):
    gap = 0
    for x in range(lo, hi):
        W, nly = step(x, gap)
        gap -= nly
        yield x, pf(x), W, nly, odd_part(W), gap


if __name__ == "__main__":
    rows = list(run())
    print("  x     pf      W    nly   odd(W)    gap")
    for x, p, W, nly, odd, gap in rows:
        print(f" {x:3} {p:6} {W:6} {nly:5} {odd:8} {gap:6}")
    g = [r[5] for r in rows]
    print(f"\n  |nly| max {max(abs(r[3]) for r in rows)},  "
          f"gap in {min(g)} .. {max(g)},  all odd parts prime: "
          f"{all(isprime(r[4]) for r in rows)}")
