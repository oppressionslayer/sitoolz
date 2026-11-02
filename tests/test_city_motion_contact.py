"""The three laws shared with the site and the engine, held to their reference files."""
import json
from fractions import Fraction
from pathlib import Path

from sitoolz import city, contact, motion

DATA = Path(__file__).parent / "data"


def test_the_city_matches_the_rust_districts():
    for name in ("district-45.rust.json", "district-0.rust.json"):
        ref = json.loads((DATA / name).read_text())
        d = city.district(int(ref["district"]))
        for axis in ("x", "z"):
            for mine, theirs in zip(d["streets"][axis], ref["streets"][axis]):
                assert (mine["lo"], mine["hi"], mine["center"], mine["width"]) == (theirs["lo"], theirs["hi"], theirs["center"], theirs["width"])
        assert (d["streets"]["width"], d["streets"]["depth"]) == (ref["streets"]["width"], ref["streets"]["depth"])
        for l, r in zip(d["lots"], ref["lots"]):
            assert l["cell"] == r["cell"] and str(l["address"]) == r["address"]
            for k in ("floors", "width", "depth", "roof", "color", "kind"):
                assert l["fields"][k] == r[k], (name, l["cell"], k)
            assert l["payload"] == r["payload"] and str(l["N"]) == r["N"]
            assert l["source"] == (r["sourceIndex"], r["sourcePrime"])
            assert l["geometry"]["bounds"] == r["bounds"]
            assert city.open_card(l["N"]) == (int(ref["district"]), l["cell"], l["fields"])
            row = city.prefix_row(l)
            assert [row[i] ^ row[i + 1] for i in range(len(row) - 1)][0] == l["address"]
    # the hash is the engine's own: `at` and `stream_of` agree with wlaw::hash on these values
    assert city.mix(0) == 0 and city.mix(1) == 1753845952
    assert city.at(45, 27, 6673, city.stream_of("city.floors")) % 12 + 1 == 4
    f = dict(city.generated(45, 27)); f["floors"] = 12
    n = city.card(45, 27, f)
    assert city.open_card(n) == (45, 27, f)
    assert city.unpack(n + 1) is None


def test_the_one_tick_law_matches_the_studio_trajectories():
    ref = json.loads((DATA / "driving-v1.javascript.json").read_text())
    for case in ref["cases"]:
        h = case["handling"]
        if h["grip"] < 1:
            continue  # the sliding cases need the full kernel; the speed law alone is checked here
        v = case["initial"]["speed"]
        tick = 0
        cps = {c["tick"]: c for c in case["checkpoints"]}
        for run in case["runs"]:
            t, s, brake = motion.decode_input(run["code"])
            for _ in range(run["ticks"]):
                upper = h["speed"] * (0.45 if t < 0 else 1)
                v = motion.longitudinal_step(v, t, h["accel"], motion.drag_of(brake, h["wet"]), motion.DT, 0.45 * h["speed"], upper)
                tick += 1
                if tick in cps:
                    assert abs(v - cps[tick]["speed"]) < 1e-9, (case["name"], tick, v, cps[tick]["speed"])
    assert abs(motion.settles_at(1, 12, 0.7) - 17.042857142857) < 1e-9
    assert motion.cap_tick(20, 16) == 363
    assert abs(motion.closed_form(0, 1, 12, 0.7, 120) - motion.run(0, 1, 12, 0.7, 120)) < 1e-9
    assert abs(motion.distance(0, 1, 12, 0.7, 120) - 4.8331361692129) < 1e-9


def test_first_contact_matches_the_studio_scenes():
    ref = json.loads((DATA / "contact-v1.reference.json").read_text())
    assert len(ref["cases"]) == 56
    for c in ref["cases"]:
        got = contact.sweep(c["motion"], c["obstacles"])
        e = c["expected"]
        assert got["kind"] == e["kind"], c["name"]
        assert got["fraction"] == Fraction(int(e["fraction"]["n"]), int(e["fraction"]["d"])), c["name"]
        assert got["ids"] == [int(i) for i in e["ids"]], c["name"]
        assert got["normals"] == [tuple(n) for n in e["normals"]], c["name"]
        assert abs(got["stop"][0] - e["stop"]["x"]) < 1e-12 and abs(got["stop"][1] - e["stop"]["z"]) < 1e-12, c["name"]


def test_moving_pairs_meet_where_the_studio_says():
    ref = json.loads((DATA / "moving-pair-v1.javascript.json").read_text())
    assert ref["profile"] == "moving-pair-exact-v1"
    fr = lambda v: Fraction(int(v["n"]), int(v["d"]))
    n = 0
    for case in ref["cases"]:
        r = case["result"]
        ids = [int(i) for i in r["ids"]]
        b_id = int(r["obstacle"]["id"]); a_id = next(i for i in ids if i != b_id)
        cars = [dict(c, id=i, hx=f["hx"], hz=f["hz"]) for c, i, f in zip(case["cars"], (a_id, b_id), r["footprints"])]
        p = contact.pair(cars[0], cars[1])
        assert p["kind"] == r["kind"] and p["fraction"] == fr(r["at"]) and p["ids"] == ids, case["name"]
        assert p["normals_on_a"] == [tuple(v) for v in r["normalsOnA"]]
        assert p["contact"] == [(fr(c["x"]), fr(c["z"])) for c in r["contact"]], case["name"]
        assert p["desired"] == [(fr(c["x"]), fr(c["z"])) for c in r["desired"]]
        # the invariants: a common displacement, and a swap
        shift = lambda c: dict(c, x=c["x"] + 7.25, z=c["z"] - 3.5, dx=c["dx"] + 2, dz=c["dz"] - 1)
        q = contact.pair(shift(cars[0]), shift(cars[1]))
        assert (q["kind"], q["fraction"], q["normals_on_a"]) == (p["kind"], p["fraction"], p["normals_on_a"])
        s = contact.pair(cars[1], cars[0])
        assert (s["kind"], s["fraction"], s["contact"]) == (p["kind"], p["fraction"], p["contact"][::-1])
        assert s["normals_on_a"] == sorted((-x, -z) for x, z in p["normals_on_a"])
        n += 1
    assert n == 9
