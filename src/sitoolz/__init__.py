"""sitoolz: Lars Rocha's number-theory playground, every module in one package.

    from sitoolz import numtools as nt, numtoolz as ntz, xortri
    from sitoolz import ntfun, ntmc, edgelaw, tri_xor, square_frame, ...

Forty-two modules and two subpackages: `primebase`, the prime engine that proves, factors
and generates up to 72 quadrillion from a 15 MB base and rolls the Pocklington ladder, and
`sequence`, the 2026-08-01 sequence project. The ones that were IPython session companions keep their names and their
aliases (`nt`, `ntz`, `ntd`, `ntm`), so a session that did `import numtools as nt` now does
`from sitoolz import numtools as nt` and nothing else changes. Submodules load on first use,
so `import sitoolz` is cheap; `pandas` is needed by most of them, `gmpy2`, `numba` and `scipy`
speed some up when present and are never required.

The notes that went with the code are in `docs/` of the source tree.
"""

import importlib

__version__ = "0.1.0"

MODULES = [
    "anchors", "antispoof", "backwalk", "centerrandom", "dfkit", "driver_lab", "edgeframe",
    "edgelaw", "edgelaw_deep", "edgeproof", "edgerandom", "edgereach", "fast30", "gap_forge",
    "ntcomb", "ntfun", "ntfunextended", "ntmc", "ntpar", "ntpar_joint", "ntprobe", "numtools",
    "numtoolsd", "numtoolsm", "numtoolz", "oeis_crosscheck", "parity_lab", "prime_check", "qnr2",
    "randbattery", "roguelike", "shopworld", "square_frame", "square_rng", "square_store",
    "squareverse", "squarewalk", "tri_frame", "tri_xor", "twin_forge", "txs", "xortri",
    "primebase", "sequence", "showcase", "twin_k2n256", "city", "motion", "contact", "district", "convoy", "zipper",
]

__all__ = list(MODULES)


def __getattr__(name):
    if name in MODULES:
        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module 'sitoolz' has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(MODULES))
