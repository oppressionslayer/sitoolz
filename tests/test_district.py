"""District drive, held to the studio's six portable cases."""
import json
from fractions import Fraction
from pathlib import Path

from sitoolz import district

DATA = Path(__file__).parent / "data"
FIELDS = ["x", "z", "yaw", "speed", "vx", "vz", "steer_vis", "lean", "squat"]


def diff(a, b):
    return max([abs(a[k] - b[k]) for k in FIELDS] + [0.0 if a["braking"] == b["braking"] else 1.0])


def test_the_six_portable_cases_replay_and_seek():
    ref = json.loads((DATA / "district-drive-v1.javascript.json").read_text())
    assert ref["profile"] == district.PROFILE
    worst, checked, contacts = 0.0, 0, 0
    for case in ref["cases"]:
        sc = case["scene"]
        edits = {e["cell"]: int(e["N"]) for e in sc["edits"]}
        s = district.scene(int(sc["district"]), edits, sc["start"], case["half"], sc["handling"])
        assert diff(s["spawn"], case["initial"]) == 0.0, case["name"]
        inputs = sc["inputs"]
        if sc["active"] == "branch":
            b = sc["branch"]
            assert district.shares_prefix(inputs, b["inputs"], b["fork"])
            assert district.fork(inputs, b["fork"], district.encode(1, 1, False), 160) == b["inputs"]
            a, t = district.trace(s, inputs), district.trace(s, b["inputs"])
            assert district.first_difference(a, t) == b["fork"] + 1
            assert a["contacts"] and not t["contacts"]
            inputs = b["inputs"]
        t = district.trace(s, inputs)
        for cp in case["checkpoints"]:
            tick, want = cp["tick"], cp["pose"]
            got = t["frames"][tick]
            d = diff(got, want)
            worst = max(worst, d)
            assert d <= 1e-9, (case["name"], tick, got, want)
            sk = district.seek(s, t, tick)
            assert sk["pose"] == got
            assert sk["card"] == district.pack_block(inputs[sk["checkpoint"]:sk["checkpoint"] + district.EVERY])
            assert sk["steps"] < district.EVERY
            checked += 1
        fc = case["firstContact"]
        if not t["contacts"]:
            assert fc is None, case["name"]
            continue
        r = t["contacts"][0]
        assert fc is not None, (case["name"], r["tick"], r["ids"])
        assert (r["tick"], r["code"], r["ids"]) == (fc["tick"], fc["code"], [int(i) for i in fc["ids"]])
        assert r["normals"] == [tuple(n) for n in fc["normals"]]
        assert r["fraction"] == Fraction(int(fc["fraction"]["n"]), int(fc["fraction"]["d"]))
        assert abs(r["used"] - fc["used"]) <= 1e-12
        assert abs(r["stop"][0] - fc["stop"]["x"]) <= 1e-9 and abs(r["stop"][1] - fc["stop"]["z"]) <= 1e-9
        contacts += 1
    assert checked == 36 and contacts == 3
    print(f"six portable cases, {checked} checkpoints, {contacts} contacts exact, worst {worst:e}")


def test_seek_lands_on_every_frame_and_the_edge_rebuilds_the_block():
    s = district.scene(45)
    t = district.trace(s, [5] * 300)
    assert t["contacts"][0]["tick"] == 208 and t["contacts"][0]["ids"] == [2907]
    for tick in range(0, 301, 7):
        assert district.seek(s, t, tick)["pose"] == t["frames"][tick]
    b = district.block(list(range(18)) + [17] * 14, 40)
    assert b["start"] == 32 and b["recovered"] == b["values"]


def test_the_lookup_fixture_cards_and_candidates():
    from sitoolz import contact, xortri
    ref = json.loads((DATA / "district-lookup-v1.javascript.json").read_text())
    assert ref["profile"] == district.PROFILE
    for b in ref["blocks"]:
        card = district.pack_block(b["input"])
        assert card == b["card"]
        assert district.card_from_json(district.card_json(card)) == card
        assert [district.command_at(card, j) for j in range(card["length"])] == b["commands"] == b["input"]
        opened = district.open_block(card)
        assert opened["edge"] == [int(e) for e in b["edge"]]
        if b["input"]:
            assert opened["edge"] == [int(v) for v in xortri.left_edge(b["input"])]
        for t in range(len(b["input"]) + 1):
            assert district.pack_block(b["input"][:t])["length"] == t
            assert int(district.pack_block(b["input"][:t])["N"]) == int(card["N"]) & ((1 << (5 * t)) - 1)
    cases = 0
    for n in ref["neighborhoods"]:
        edits = {e["cell"]: int(e["N"]) for e in n["edits"]}
        s = district.scene(int(n["district"]), edits)
        for i, sol in enumerate(n["solids"]):
            assert s["solids"][i]["id"] == sol["id"] and [s["solids"][i][k] for k in "x0 x1 z0 z1".split()] == [sol[k] for k in "x0 x1 z0 z1".split()]
        for c in n["cases"]:
            got = sorted(int(s["solids"][i]["id"]) for i in district.lots_near(s, c["motion"]))
            assert got == sorted(int(i) for i in c["ids"]), (n["district"], c["motion"])
            cases += 1
    assert cases == 20
    w = ref["wallContact"]
    s = district.scene(int(w["district"]))
    idx = district.lots_near(s, w["motion"])
    assert sorted(int(s["solids"][i]["id"]) for i in idx) == sorted(int(i) for i in w["ids"])
    r = contact.sweep(w["motion"], [s["solids"][i] for i in idx])
    assert r["kind"] == "hit" and [str(i) for i in r["ids"]] == w["result"]["ids"]
    assert r["fraction"] == Fraction(int(w["result"]["fraction"]["n"]), int(w["result"]["fraction"]["d"]))
