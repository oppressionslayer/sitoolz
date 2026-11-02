"""edgerandom -- which parts of rule 30 are predictable, and how fast.

Rule 30 has been used as a randomness source.  Wolfram proposed its **centre
column** for that purpose and Mathematica shipped it as a default RNG for
years; the column's non-periodicity is still an open problem, and Rowland's
Conjecture 3 is about exactly this kind of reducibility.

The **right light-cone edge** is a different story, and this module is the
caveat made runnable:

    edge(g) = TABLE[v2(g+1)]

A table index.  No simulation, at any generation, in under a microsecond.  The
edge of rule 30 is not merely compressible -- it is a closed form.

    from sitoolz import edgerandom as er

    er.edge_at(10**9)            # 51 -- constant time
    er.is_predictable(10**9)     # what is known, and how cheaply
    er.speed_compare(20000)      # lookup vs actually simulating
    er.caveat_frame()            # part by part: predictable or not
    er.edge_stream(0, 32)        # bulk, still O(1) per generation

**What this does and does not say.**  It says the edge is decided by a lookup,
so any sampling scheme that touches it is drawing from a known sequence.  It
says nothing about the centre column, which remains open -- and the contrast is
the point: two features of the same automaton, one settled in microseconds and
one unsettled for forty years.

Proof of the law: `PROOF_V2_LAW.md`.  Table provenance: `edgelaw`, extended in
`edgelaw_deep`.  Stdlib only.
"""
import time

from sitoolz import edgelaw as _el

try:
    from sitoolz import edgelaw_deep as _deep
except ImportError:                                    # pragma: no cover
    _deep = None


def _frame(rows, index=None):
    try:
        import pandas as pd
        df = pd.DataFrame(rows)
        return df.set_index(index) if index else df
    except ImportError:
        from sitoolz.dfkit import frame
        return frame(rows, index=index)


def v2(n):
    """The 2-adic valuation -- how many times 2 divides n."""
    n = int(n)
    if n == 0:
        return 0
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    return c


def _table():
    """The deepest verified table available."""
    if _deep is not None and len(_deep.TABLE) > len(_el.TABLE):
        return _deep.TABLE
    return _el.TABLE


def edge_at(g):
    """Edge run at generation `g`. One lookup. Raises past the verified table.

    This is the whole caveat in one call: whatever `g` is, the answer costs a
    valuation and a dict access.
    """
    T = _table()
    v = v2(int(g) + 1)
    if v not in T:
        raise KeyError(
            f"generation {int(g):,} has v2(g+1) = {v}, past the verified table "
            f"(v2 <= {max(T)}). Extend it rather than extrapolating.")
    return T[v]


def is_predictable(g):
    """Is the edge at `g` decided without simulating? What does it cost?

    Returns a dict: the rung, the value, the cost, and whether the table
    reaches that far. Nothing here estimates or extrapolates.
    """
    T = _table()
    v = v2(int(g) + 1)
    known = v in T
    t0 = time.perf_counter()
    val = T.get(v)
    dt = time.perf_counter() - t0
    return {"generation": int(g), "rung_v2": v, "in_table": known,
            "edge_run": val, "seconds": dt,
            "method": "table lookup" if known else "not verified to this rung",
            "simulation_needed": False if known else True}


def speed_compare(g, width=None):
    """Lookup versus actually running the automaton to generation `g`.

    Simulation is O(g) in time; the lookup is O(1). The ratio is the caveat.
    """
    g = int(g)
    t0 = time.perf_counter()
    look = edge_at(g)
    t_look = time.perf_counter() - t0

    W = width or (g + 8)
    mask = (1 << W) - 1
    t1 = time.perf_counter()
    E = 1
    for _ in range(g):
        E = (E ^ ((E << 1) | (E << 2))) & mask
    r = 0
    while (E >> r) & 1:
        r += 1
    t_sim = time.perf_counter() - t1

    return {"generation": g, "lookup": look, "simulated": r,
            "agree": look == r,
            "lookup_us": round(t_look * 1e6, 3),
            "simulate_us": round(t_sim * 1e6, 1),
            "speedup": round(t_sim / t_look) if t_look else None}


def edge_stream(lo, hi):
    """Edge run for every generation in [lo, hi). O(1) each, no state carried."""
    return [edge_at(g) for g in range(int(lo), int(hi))]


def caveat_frame():
    """Rule 30 part by part: what is settled and what is not.

    The contrast is the content -- the centre column is the part proposed as a
    randomness source, and it is the part still open.
    """
    T = _table()
    rows = [
        {"feature": "right light-cone edge (run length)",
         "status": "CLOSED FORM",
         "how": "edge(g) = TABLE[v2(g+1)]",
         "cost": "one lookup, any g",
         "note": f"proved; table verified to rung {max(T)}"},
        {"feature": "edge, as a bit vector",
         "status": "exactly periodic",
         "how": "each bit has period 2^p(b), from generation 0",
         "cost": "O(width)",
         "note": "no transient -- not merely eventual"},
        {"feature": "period spectrum p(b)",
         "status": "not 2-automatic",
         "how": "parity of the driver weight",
         "cost": "O(2^p) per coordinate",
         "note": "2-kernel grows 3, 7, 15 -- no digit-based closed form"},
        {"feature": "centre column",
         "status": "OPEN",
         "how": "no known formula",
         "cost": "simulation",
         "note": "proposed as an RNG; non-periodicity unproven"},
        {"feature": "left side (column periods)",
         "status": "partly settled",
         "how": "Rowland 2006, Lemma 2 and Proposition 2",
         "cost": "-",
         "note": "doubling criterion needs a white stripe"},
    ]
    return _frame(rows, index="feature")


def rng_warning(samples):
    """How many of these generations have an edge run already determined?

    Pass the generations a sampling scheme would touch. Any that come back
    covered are drawn from a sequence that is one lookup away.
    """
    T = _table()
    covered = [g for g in samples if v2(int(g) + 1) in T]
    return {"sampled": len(samples), "covered_by_table": len(covered),
            "fraction": round(len(covered) / len(samples), 4) if samples else 0.0,
            "values": [edge_at(g) for g in covered[:12]]}


if __name__ == "__main__":
    print("rule 30, part by part:")
    print(caveat_frame())
    print()
    for g in (255, 4095, 10 ** 6, 10 ** 9):
        d = is_predictable(g)
        print(f"  g={g:>12,}  v2={d['rung_v2']:>2}  edge={d['edge_run']}"
              f"  {d['seconds']*1e6:.2f} us  ({d['method']})")
    print()
    print("lookup vs simulation:")
    for g in (2000, 20000):
        c = speed_compare(g)
        print(f"  g={g:>6}: lookup {c['lookup_us']:>7.3f} us   "
              f"simulate {c['simulate_us']:>10.1f} us   "
              f"speedup {c['speedup']:>8,}x   agree={c['agree']}")
