"""The convoy reserve, held to the studio's eight receipts."""
import json
from fractions import Fraction
from pathlib import Path

from sitoolz import convoy

DATA = Path(__file__).parent / "data"


def test_the_eight_receipts_replay():
    ref = json.loads((DATA / "convoy-brake-v1.javascript.json").read_text())
    assert ref["profile"] == convoy.PROFILE
    worst = 0.0
    for case in ref["cases"]:
        r = case["receipt"]; rc = r["recipe"]
        lead, foll = rc["cars"]
        p = convoy.plan(rc["mode"], rc["leadTick"], rc["delay"], rc["gap"], rc["margin"], rc["top"], rc["accel"],
                        lead["speed"], foll["speed"], lead["wet"], foll["wet"])
        got = convoy.run(p)
        for key, want in (("required", r["required"]), ("tail", r["remainingFollowerTravel"]), ("peak", r["peakClosure"]), ("minimum_gap", r["minimumGap"])):
            d = abs(got[key] - want); worst = max(worst, d)
            assert d <= 1e-9, (case["name"], key, got[key], want)
        assert got["fitted"] == r["fittedGap"] and got["closest_tick"] == r["closestTick"], case["name"]
        for cp in r["checkpoints"]:
            k = cp["tick"]
            for car, want in zip(("leader", "follower"), cp["cars"]):
                assert abs(got[car]["x"][k] - want["distance"]) <= 1e-9 and abs(got[car]["speed"][k] - want["speed"]) <= 1e-9, (case["name"], k, car)
            assert abs(convoy.gap_at(got, k) - cp["gap"]) <= 1e-9 and abs(convoy.reserve_at(got, k) - cp["reserve"]) <= 1e-9, (case["name"], k)
        c, want = got["contact"], r["contact"]
        if want is None:
            assert c is None, case["name"]
        else:
            assert c is not None and c["tick"] == want["tick"] and c["kind"] == want["kind"], case["name"]
            assert c["fraction"] == Fraction(int(want["fraction"]["n"]), int(want["fraction"]["d"])), case["name"]
            assert abs(c["seconds"] - want["seconds"]) <= 1e-12 and abs(c["motion"]["x"] - want["motion"]["x"]) <= 1e-9
            assert c["normals"] == [tuple(n) for n in want["normals"]]
    print(f"eight receipts, worst difference {worst:e}")
