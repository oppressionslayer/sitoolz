import doctest
import importlib

import numpy as np
import pytest

import sitoolz
from sitoolz import MODULES


@pytest.mark.parametrize("name", MODULES)
def test_every_module_imports(name):
    mod = importlib.import_module(f"sitoolz.{name}")
    assert getattr(sitoolz, name) is mod
    if hasattr(mod, "__path__"):  # a subpackage: its modules are listed
        assert mod.__all__, f"{name} lists no modules"
        return
    public = [n for n in dir(mod) if not n.startswith("_") and callable(getattr(mod, n))]
    assert public, f"{name} has no public functions"


def test_doctests_in_numtoolz_and_xortri():
    for name in ("numtoolz", "xortri"):
        r = doctest.testmod(importlib.import_module(f"sitoolz.{name}"), verbose=False)
        assert r.failed == 0, f"{name}: {r.failed} doctest failures"


def test_the_session_aliases_still_work():
    from sitoolz import numtools as nt, numtoolz as ntz, numtoolsd as ntd, numtoolsm as ntm
    idx = ntz.structured_index(128)
    T = ntz.buildsystem(idx)
    assert np.array_equal(T, T[0][:, None] ^ T[0][None, :]) and np.array_equal(ntz.recover_index(T), idx)
    assert nt.verdict(91)[0] == "composite" and nt.is_prime(2**61 - 1)
    assert nt.ruler(8).tolist() == [0, 1, 0, 2, 0, 1, 0, 3]
    assert ntd and ntm


def test_the_xor_triangle_is_the_session():
    from sitoolz import xortri
    a = [2, 6, 8, 11, 12, 5, 4]
    assert xortri.square(a)[:, -1].tolist() == [1, 8, 6, 12, 15, 2]
    assert xortri.apex(a) == 2 and xortri.from_left_edge(xortri.left_edge(a)).tolist() == a


def test_nt_and_ntz_share_no_names():
    from sitoolz import numtools as nt, numtoolz as ntz
    a = {n for n in dir(nt) if not n.startswith("_") and callable(getattr(nt, n))}
    b = {n for n in dir(ntz) if not n.startswith("_") and callable(getattr(ntz, n))}
    assert not (a & b)
