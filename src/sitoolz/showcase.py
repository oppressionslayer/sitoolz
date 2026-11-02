"""showcase -- the numbers the sessions built, shipped with the package so you can hold them.

    from sitoolz import showcase
    p = showcase.big_prime()          # the 19,538-digit prime built from M127 in nine rolls
    showcase.twins(4096)              # the k*2^4096 +/- 1 hunt: 121 proven +1 sides, one pair
    showcase.chain_4096()             # a 4096-bit prime with its three-rung Pocklington chain

Every value here was built, not found: each one carries the construction that proves it.
"""

import json
import sys
from pathlib import Path

# a 19,538-digit integer is past Python's default string-conversion guard
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

_DATA = Path(__file__).parent / "data"


def big_prime():
    """The 19,538-digit prime (64,903 bits) rolled up from M127, as an int."""
    return int((_DATA / "rollball_M127_19538.dat").read_text().strip())


def twins(n):
    """The hunt records for k*2^n +/- 1, n in (2048, 4096, 8192): a list of dicts with
    k, n, digits, upper, upper_proven_pocklington, lower, lower_status."""
    rows = json.loads((_DATA / f"twins_{n}.json").read_text())
    for r in rows:
        r["upper"] = int(r["upper"])
        r["lower"] = int(r["lower"])
    return rows


def twin_pairs(n):
    """The records of `twins(n)` where the -1 side also passed the probable-prime screen."""
    return [r for r in twins(n) if r["lower_status"] != "composite"]


def chain_4096():
    """A 4096-bit prime and the Pocklington chain that proves it, from the seed 112543."""
    d = json.loads((_DATA / "prime_4096_chain.json").read_text())
    chain = d["chain"]
    if isinstance(chain, str):
        import ast
        chain = ast.literal_eval(chain)
    return {"bits": d["bits"], "prime": int(d["prime"]), "chain": [[int(x) if x is not None else None for x in rung] for rung in chain]}
