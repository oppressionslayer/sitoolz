"""contact -- first contact: where a fixed-orientation footprint first enters a solid along a
straight segment, as an exact fraction. The browser's First Contact contract, version 1, and the
engine's ``wcraft::contact``, in Python's own exact rationals. A float enters as the binary value
it is (``Fraction(x)``), so a live position off the quarter grid is compared exactly, the way the
studio decodes binary64 into dyadic integers.

    from sitoolz import contact
    contact.sweep({"x": 0, "z": 0, "dx": 80, "dz": 0, "hx": 2, "hz": 1},
                  [{"id": "2", "x0": 20, "x1": 20.25, "z0": -6, "z1": 6}])
    # kind 'hit', fraction 9/40, ids ['2'], normals [(-1, 0)], stop (18.0, 0.0)
"""
from fractions import Fraction


def _axis(p, d, lo, hi):
    """The open interval of times p + t*d is strictly inside (lo, hi): (t0, t1), None when never;
    (None, None) means always, for a still axis already inside."""
    if d == 0:
        return (None, None) if lo < p < hi else None
    a, b = Fraction(lo - p) / Fraction(d), Fraction(hi - p) / Fraction(d)
    return (a, b) if a < b else (b, a)


def sweep(motion, solids):
    m = {k: Fraction(v) for k, v in motion.items()}
    desired = (float(m["x"] + m["dx"]), float(m["z"] + m["dz"]))
    overlap, best, hits = [], None, []
    for s in solids:
        x0, x1 = Fraction(s["x0"]) - m["hx"], Fraction(s["x1"]) + m["hx"]
        z0, z1 = Fraction(s["z0"]) - m["hz"], Fraction(s["z1"]) + m["hz"]
        ix = _axis(m["x"], m["dx"], x0, x1)
        iz = _axis(m["z"], m["dz"], z0, z1)
        if ix is None or iz is None:
            continue
        starts = [t for t in (ix[0], iz[0]) if t is not None]
        ends = [t for t in (ix[1], iz[1]) if t is not None]
        entry = max(starts) if starts else Fraction(-1)
        exit_ = min(ends) if ends else Fraction(2)
        if not entry < exit_ or not exit_ > 0:
            continue
        if entry < 0:
            overlap.append(int(s["id"]))
            continue
        if entry > 1:
            continue
        normals = []
        if ix[0] is not None and ix[0] == entry and m["dx"] != 0:
            normals.append((-1 if m["dx"] > 0 else 1, 0))
        if iz[0] is not None and iz[0] == entry and m["dz"] != 0:
            normals.append((0, -1 if m["dz"] > 0 else 1))
        if best is None or entry < best:
            best, hits = entry, [(int(s["id"]), normals)]
        elif entry == best:
            hits.append((int(s["id"]), normals))
    if overlap:
        return {"kind": "overlap", "fraction": Fraction(0), "ids": sorted(set(overlap)), "normals": [],
                "stop": (float(m["x"]), float(m["z"])), "desired": desired}
    if best is None:
        return {"kind": "clear", "fraction": Fraction(1), "ids": [], "normals": [], "stop": desired, "desired": desired}
    ids = sorted({h[0] for h in hits})
    normals = sorted({n for h in hits for n in h[1]})
    stop = (float(m["x"] + m["dx"] * best), float(m["z"] + m["dz"] * best))
    return {"kind": "hit", "fraction": best, "ids": ids, "normals": normals, "stop": stop, "desired": desired}


def pair(a, b):
    """**Two moving cars as one sweep**, the studio's ``moving-pair-exact-v1``. Each car is a dict
    ``id, x, z, dx, dz, hx, hz``. The relative motion of A against B, with B as a solid of its own
    half-sizes at the origin, goes through ``sweep``; the exact fraction then places both cars on
    their own world paths. Normals point from B toward A. Adding one displacement to both cars
    leaves the fraction and normals alone; swapping the cars swaps the positions and negates the
    normals.

    >>> A = {"id": 45, "x": -22, "z": 0, "dx": 40, "dz": 0, "hx": 2.5, "hz": 1.25}
    >>> B = {"id": 7, "x": 0, "z": -22, "dx": 0, "dz": 40, "hx": 1.25, "hz": 2.75}
    >>> r = pair(A, B); r["kind"], r["fraction"], r["contact"][0], r["normals_on_a"]
    ('hit', Fraction(73, 160), (Fraction(-15, 4), Fraction(0, 1)), [(-1, 0)])
    """
    relative = {"x": a["x"] - b["x"], "z": a["z"] - b["z"], "dx": a["dx"] - b["dx"], "dz": a["dz"] - b["dz"],
                "hx": a["hx"], "hz": a["hz"]}
    obstacle = {"id": str(b["id"]), "x0": -b["hx"], "x1": b["hx"], "z0": -b["hz"], "z1": b["hz"]}
    hit = sweep(relative, [obstacle])
    t = hit["fraction"]
    along = lambda c, t: (Fraction(c["x"]) + t * Fraction(c["dx"]), Fraction(c["z"]) + t * Fraction(c["dz"]))
    return {"kind": hit["kind"], "fraction": t, "ids": sorted([int(a["id"]), int(b["id"])]),
            "normals_on_a": hit["normals"], "relative": relative, "obstacle": obstacle,
            "contact": [along(a, t), along(b, t)], "desired": [along(a, 1), along(b, 1)]}
