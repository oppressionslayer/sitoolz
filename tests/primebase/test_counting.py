"""li must not lie.

The Ramanujan series ``li(x) = gamma + ln(z) + sum z**k/(k*k!)`` has terms that PEAK
at ``k ~ z = ln(x)`` and only decay after. A fixed iteration cap therefore truncates
mid-climb once ``z`` reaches it and returns a confidently wrong low answer -- no
exception, no warning. That shipped: the cap was 400, so every ``x`` between about
1e173 and 1e308 came back ~45% low, and past 1e308 ``float(x)`` raised OverflowError.

These pin both, plus the accuracy floor, so the estimate can be trusted at the sizes
this library actually chases.
"""
import math

import pytest

import sitoolz.primebase as pb

gmpy2 = pytest.importorskip("gmpy2", reason="ground-truth li needs MPFR")


def true_li(x, bits=12000):
    """li(x) = Ei(ln x), to `bits` of precision. MPFR is the reference."""
    gmpy2.get_context().precision = bits
    return gmpy2.eint(gmpy2.log(gmpy2.mpz(x)))


def agree_digits(a, b, bits=12000):
    """Leading decimal digits shared by two positive values (0 if magnitudes differ)."""
    gmpy2.get_context().precision = bits
    ma, ea, _ = gmpy2.digits(gmpy2.mpfr(a), 10)
    mb, eb, _ = gmpy2.digits(gmpy2.mpfr(b), 10)
    if ea != eb:
        return 0
    n = 0
    for x, y in zip(ma, mb):
        if x != y:
            break
        n += 1
    return n


# (x, worst acceptable relative error against EXACT pi) -- li overshoots pi throughout
PI_ACCURACY = [(10 ** 4, 2e-2), (10 ** 6, 2e-3), (10 ** 8, 2e-4), (10 ** 11, 1e-5)]


@pytest.mark.parametrize("x,tol", PI_ACCURACY)
def test_li_tracks_exact_pi(x, tol):
    exact = pb.pi(x)
    assert abs(pb.li(x) - exact) / exact < tol


def test_li_overshoots_pi():
    """li(x) > pi(x) everywhere in reach -- the sieve upper bound depends on this."""
    for e in range(4, 12):
        x = 10 ** e
        assert pb.li(x) > pb.pi(x), f"li undershot pi at 1e{e}"


@pytest.mark.parametrize("e", [4, 8, 11, 50, 100, 173, 300, 500, 1000, 2441])
def test_li_matches_mpfr(e):
    """The regression that matters: 1e173..1e308 used to be silently ~45% low."""
    x = 10 ** e
    assert agree_digits(pb.li(x), true_li(x)) >= 12


def test_li_survives_past_float_range():
    """float(x) used to raise OverflowError here; pi of a 2441-digit x is a 2438-digit int."""
    v = pb.li(10 ** 2441)
    assert isinstance(v, int)
    assert len(str(v)) == 2438


def test_li_returns_float_while_one_fits():
    assert isinstance(pb.li(10 ** 11), float)
    assert isinstance(pb.li(1e11), float)


def test_li_big_carries_16_digits():
    """_li_big computes C/z in float, so ~16 significant digits -- know the ceiling."""
    x = 10 ** 2441
    assert agree_digits(pb.li(x), true_li(x)) >= 15


def test_li_edge_values():
    assert pb.li(1) == 0.0
    assert pb.li(0) == 0.0
    assert pb.li(-5) == 0.0
    assert pb.li(1.0) == 0.0


def test_wheel_law_pi_over_4():
    """The 2,10 wheel holds pi(top)/4 -- the predictor, and it leans on pi/li."""
    import numpy as np
    r = np.arange(2, 2 ** 16 * 10, 10)
    X = r ^ 5
    top = int(X.max())
    measured = sum(pb.is_prime(int(v)) for v in X)
    assert abs(measured - pb.pi(top) / 4) / (pb.pi(top) / 4) < 0.01
