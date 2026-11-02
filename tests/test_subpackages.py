import importlib
import pytest


def test_primebase_proves_and_rolls():
    from sitoolz import primebase as pb
    assert [int(x) for x in pb.first(5)] == [2, 3, 5, 7, 11]
    assert int(pb.nth(1_000_000)) == 15485863
    assert pb.is_prime(9_999_999_999_999_917) and not pb.is_prime(72_092_246_240_999_997)
    steps = pb.prove.build_up(2**31 - 1, 512)
    assert pb.prove.verify_certificate(steps)


def test_sequence_imports_when_sympy_is_there():
    pytest.importorskip("sympy")
    for m in ("better", "frame", "regen", "rulefinder", "twotables", "xorsquares"):
        importlib.import_module(f"sitoolz.sequence.{m}")


def test_certificates_cross_the_language_line():
    from pathlib import Path
    from sitoolz.primebase import prove
    chain = prove.build_up(2**31 - 1, 256)
    for conv in ("numtools", "slip"):
        text = prove.to_json(chain, conv)
        back = prove.from_json(text)
        assert [s["N"] for s in back] == [s["N"] for s in chain]
        assert prove.verify_certificate(back)
    data = Path(__file__).parent / "data"
    for name in ("rust_chain_512.json", "rust_chain_512_slip.json"):
        text = (data / name).read_text()
        steps = prove.from_json(text)
        assert steps[0]["N"] == 2**31 - 1 and steps[-1]["N"].bit_length() >= 512
        assert prove.verify_certificate(steps), f"{name}: every rung Rust rolled re-proves in Python"
    # a tampered witness fails
    steps = prove.from_json((data / "rust_chain_512.json").read_text())
    steps[1]["a"] = 1
    assert not prove.verify_certificate(steps)
