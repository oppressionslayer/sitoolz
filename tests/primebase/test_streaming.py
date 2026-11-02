"""The decoders must not scale with the base's size.

A decoded prime costs 8 bytes, so an 11.8-billion-prime base is 94 GB -- decoding it
whole OOM-kills the process. These tests pin the two properties that prevent that:
every library path streams, and the two array-returning functions refuse rather than
try. They also pin ``mod``, since a big base's reach runs past int64 and plain
``n % chunk`` raises ``OverflowError`` there.
"""
import numpy as np
import pytest

import sitoolz.primebase as pb
from sitoolz.primebase import _base

REF_LIMIT = 10 ** 7          # 620k primes: a reference small enough to hold outright


@pytest.fixture(scope="module")
def ref():
    return _base.base_upto(REF_LIMIT)


def test_stream_reproduces_base_upto(ref):
    assert np.array_equal(np.concatenate(list(_base.stream(REF_LIMIT))), ref)


@pytest.mark.parametrize("limit", [0, 1, 2, 3, 4, 5, 6, 7, 100, 10 ** 6])
def test_stream_edge_limits(ref, limit):
    """The 2/3 head is spliced in by hand, so the small limits are the risky ones."""
    chunks = list(_base.stream(limit))
    got = np.concatenate([np.empty(0, dtype=np.int64)] + chunks)
    assert np.array_equal(got, ref[ref <= limit])


def test_stream_chunks_stay_small():
    """The point of the whole exercise: chunk size is bounded, not base-proportional."""
    for i, c in zip(range(4), _base.stream(chunk=1 << 20)):
        assert c.size <= (1 << 20)


def test_count_upto_matches_searchsorted(ref):
    for x in list(range(0, 3000)) + [10 ** 5, 10 ** 6, REF_LIMIT - 1]:
        assert _base.count_upto(x) == int(np.searchsorted(ref, x, side="right")), x


def test_count_upto_spans_whole_base():
    assert _base.count_upto(pb.largest()) == pb.count()
    assert _base.count_upto(pb.largest() + 1) == pb.count()


def test_pi_uses_metadata_at_the_top():
    assert pb.pi(pb.largest()) == pb.count()


@pytest.mark.parametrize("bits", [40, 62, 63, 64, 80, 200, 700])
def test_mod_matches_python_past_int64(ref, bits):
    """``n % chunk`` raises OverflowError past 2**63; ``mod`` must not."""
    n = 0xDEADBEEF * (1 << bits) + 12345
    chunk = ref[1000:6000]
    assert np.array_equal(_base.mod(n, chunk),
                          np.array([n % int(p) for p in chunk], dtype=np.int64))


def test_mod_fast_path_still_exact(ref):
    n = 2 ** 61 + 7                       # below the 2**62 switch: plain numpy path
    chunk = ref[:5000]
    assert np.array_equal(_base.mod(n, chunk), n % chunk)


def test_base_refuses_instead_of_oom(monkeypatch):
    """base() must raise, not allocate, when the base is too big to decode."""
    monkeypatch.setattr(_base, "_MAX_DECODE", 10)
    pb.base.cache_clear()
    try:
        with pytest.raises(MemoryError, match="stream"):
            pb.base()
    finally:
        pb.base.cache_clear()


def test_base_upto_refuses_instead_of_oom(monkeypatch):
    monkeypatch.setattr(_base, "_MAX_DECODE", 10)
    with pytest.raises(MemoryError, match="stream"):
        _base.base_upto(REF_LIMIT)


def test_base_upto_under_the_cap_still_works(monkeypatch):
    monkeypatch.setattr(_base, "_MAX_DECODE", 10)
    assert _base.base_upto(10).tolist() == [2, 3, 5, 7]     # 4 primes, under the cap


def test_hot_paths_never_call_base(monkeypatch):
    """factor/first/nth/pi must not route through the un-streamed decoder."""
    def boom():
        raise AssertionError("called base() -- this is the 94 GB path")

    monkeypatch.setattr(_base, "base", boom)
    monkeypatch.setattr(pb, "base", boom)
    assert pb.factor(10 ** 12 + 39) is None
    assert pb.first(5).tolist() == [2, 3, 5, 7, 11]
    assert pb.nth(100) == 541
    assert pb.pi(10 ** 6) == 78_498
    assert pb.largest_factor(2 * 3 * 5 * 7 * 11 * 13) == 13
    assert pb.factorize(360) == {2: 3, 3: 2, 5: 1}
