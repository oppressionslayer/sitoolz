"""sitoolz.sequence: the onlyafter/newballer logs, the A' branch, and the XOR squares.

The sequence project of 2026-08-01: `regen` regenerates both log formats from closed forms,
`twotables` prints the parity-bit bridge, `xorsquares` scans n ^ k for squares (a bijection,
n = m^2 ^ 5), `rulefinder` searches over rules instead of writing a generator, `frame` and
`better` carry the product form pf = 2·x·futurex. Needs `sympy`
(`pip install "sitoolz[sequence]"`). The notes are in `docs/sequence/`.
"""
__all__ = ["better", "frame", "regen", "rulefinder", "twotables", "xorsquares"]
