"""The docs stay in the public voice, point at shipped code, and link to files that exist."""
import importlib
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = sorted((ROOT / "docs").rglob("*.md")) + [ROOT / "README.md"]

# the banned strings are assembled at run time so this file never contains them itself
BANNED = [
    "honest", "Honest", "Lars Lars", "Wil" + "liam", "Wil" + "iam", "wil" + "iam", "~/",
    "/ho" + "me/", "/Us" + "ers/", "lars_" + "approves", "In [", "❯", "✻", "not a record",
    "Why not just", "not for circulation", "TOMO" + "RROW", "ton" + "ight", "Ton" + "ight",
]

MODULE_RE = re.compile(r"\bsitoolz\.([a-zA-Z_][a-zA-Z0-9_]*)(?:\.([a-zA-Z_][a-zA-Z0-9_]*))?")
IMPORT_RE = re.compile(r"^\s*from sitoolz(?:\.([a-zA-Z_][a-zA-Z0-9_.]*))? import (.+)$", re.M)
LINK_RE = re.compile(r"\]\(([^)#]+)(?:#[^)]*)?\)")


def _fences(text):
    return re.findall(r"```(?:python)?\n(.*?)```", text, flags=re.S)


@pytest.mark.parametrize("doc", DOCS, ids=[str(d.relative_to(ROOT)) for d in DOCS])
def test_no_banned_phrasing(doc):
    text = doc.read_text(encoding="utf-8")
    hits = [b for b in BANNED if b in text]
    assert not hits, f"{doc.name}: {hits}"


@pytest.mark.parametrize("doc", DOCS, ids=[str(d.relative_to(ROOT)) for d in DOCS])
def test_code_points_at_shipped_modules(doc):
    import sitoolz
    text = doc.read_text(encoding="utf-8")
    for fence in _fences(text):
        for m in IMPORT_RE.finditer(fence):
            sub, names = m.group(1), m.group(2).split("#")[0]
            if sub:
                importlib.import_module(f"sitoolz.{sub}")
                mod = importlib.import_module(f"sitoolz.{sub}")
                for n in names.split(","):
                    n = n.strip().split(" as ")[0].strip().strip("()")
                    if n and n != "*":
                        assert hasattr(mod, n), f"{doc.name}: sitoolz.{sub} has no {n}"
            else:
                for n in names.split(","):
                    n = n.strip().split(" as ")[0].strip().strip("()")
                    if n:
                        assert n in sitoolz.MODULES, f"{doc.name}: sitoolz has no module {n}"
        for m in MODULE_RE.finditer(fence):
            top, attr = m.group(1), m.group(2)
            assert top in sitoolz.MODULES, f"{doc.name}: sitoolz.{top} is not a module"
            if attr:
                mod = importlib.import_module(f"sitoolz.{top}")
                assert hasattr(mod, attr), f"{doc.name}: sitoolz.{top} has no {attr}"


@pytest.mark.parametrize("doc", DOCS, ids=[str(d.relative_to(ROOT)) for d in DOCS])
def test_links_resolve(doc):
    text = doc.read_text(encoding="utf-8")
    for target in LINK_RE.findall(text):
        if target.startswith(("http://", "https://", "mailto:")):
            continue
        if " " in target or "·" in target or not re.search(r"\.(md|dat|json|py|txt)$|/", target):
            continue  # math written as [a](b), not a link
        assert (doc.parent / target).exists(), f"{doc.name}: broken link {target}"


def test_twin_k2n256_verifies():
    from sitoolz import twin_k2n256
    r = twin_k2n256.verify()
    assert r["digits"] == 82 and r["difference"] == 2
    assert r["lucas_selftest"] and r["upper_proven"] and r["lower_proven"]
    assert r["witnesses"] == {2: 5, 3: 2, 31: 2, 107: 2}


def test_showcase_loads():
    from sitoolz import showcase
    p = showcase.big_prime()
    assert len(str(p)) == 19538 and p.bit_length() == 64903
    assert str(p).startswith("488722292701599022162589797030") and str(p).endswith("178852201085316611321153499213")
    pairs = {n: showcase.twin_pairs(n) for n in (2048, 4096, 8192)}
    assert [(r["k"], r["digits"]) for r in pairs[2048]] == [(442533, 623)]
    assert [(r["k"], r["digits"]) for r in pairs[4096]] == [(2130693, 1240)]
    assert pairs[8192] == []
    assert len(showcase.twins(4096)) == 121 and all(r["upper_proven_pocklington"] for r in showcase.twins(4096))
    c = showcase.chain_4096()
    assert c["bits"] == 4096 and [r[0].bit_length() for r in c["chain"]][-1] == 4096
