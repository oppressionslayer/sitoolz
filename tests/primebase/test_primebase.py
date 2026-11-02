import numpy as np
import pytest

import sitoolz.primebase as pb


def test_first():
    assert pb.first(10).tolist() == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert pb.first(1)[0] == 2


def test_nth():
    assert pb.nth(1) == 2
    assert pb.nth(6) == 13
    assert pb.nth(1_000_000) == 15_485_863   # known value


def test_primes_range_matches_reference():
    got = pb.primes(1000, 1100)
    ref = [n for n in range(1000, 1100)
           if all(n % d for d in range(2, int(n ** 0.5) + 1))]
    assert got.tolist() == ref


def test_primes_window_far_out():
    # a window past a trillion -- seeded by the base, no full table
    got = pb.primes(10 ** 12, 10 ** 12 + 100)
    assert got.tolist() == [1000000000039, 1000000000061,
                            1000000000063, 1000000000091]


def test_twins():
    assert pb.twins(0, 20).tolist() == [[3, 5], [5, 7], [11, 13], [17, 19]]


def test_is_prime_and_factor():
    assert pb.is_prime(9_999_999_999_999_917) is True
    assert pb.is_prime(9_999_999_999_999_919) is False
    assert pb.factor(72_092_246_240_999_997) == 3
    assert pb.factor(97) is None


def test_reach_guard():
    with pytest.raises(ValueError):
        pb.is_prime(pb.REACH + 1)


def test_base_upto_prefix():
    up = pb.base_upto(100)
    assert up.tolist() == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41,
                           43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
