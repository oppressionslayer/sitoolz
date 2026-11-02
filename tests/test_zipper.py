"""The zipper gate, held to the studio's eight receipts."""
import json
from fractions import Fraction
from pathlib import Path

from sitoolz import zipper

DATA = Path(__file__).parent / "data"
KEYS = ["order", "actor", "label", "lane", "position", "ready", "tick", "exit", "wait", "eligible", "reason"]


def test_the_eight_receipts_replay():
    ref = json.loads((DATA / "zipper-gate-v1.javascript.json").read_text())
    assert ref["profile"] == zipper.PROFILE
    for case in ref["cases"]:
        r = case["receipt"]; rc = r["recipe"]
        ready = [[c["ready"] for c in lane] for lane in rc["lanes"]]
        cards = [[int(c["card"]["N"]) for c in lane] for lane in rc["lanes"]]
        p = zipper.plan(ready, cards, rc["startLane"], rc["policy"], rc["speed"], rc["gateTicks"], rc["gapMode"],
                        rc["gap"], rc["delay"], rc["wet"], rc["margin"])
        got = zipper.run(p)
        assert (got["gap"], got["half_max"], got["clearance_ticks"], got["stride"], got["end"]) == \
               (r["gap"], r["halfMax"], r["clearanceTicks"], r["stride"], r["end"]), case["name"]
        if r["brake"] is None:
            assert got["brake"] is None
        else:
            assert abs(got["brake"]["required"] - r["brake"]["required"]) <= 1e-9 and got["brake"]["fit"] == r["brake"]["fit"]
            assert abs(got["brake"]["tail"] - r["brake"]["tail"]) <= 1e-12
        assert len(got["events"]) == len(r["events"]), case["name"]
        for e, w in zip(got["events"], r["events"]):
            for k in KEYS:
                assert e[k] == w[k], (case["name"], e["label"], k, e[k], w[k])
            assert str(e["address"]) == w["address"]
            assert (e["preferred_before"], e["slot_open"], e["arrival_idle"]) == (w["preferredBefore"], w["slotOpen"], w["arrivalIdle"]), (case["name"], e["label"])
            if w["pairGap"] is None:
                assert e["pair_gap"] is None
            else:
                assert e["pair_gap"] == Fraction(int(w["pairGap"]["n"]), int(w["pairGap"]["d"])), (case["name"], e["label"])
                assert e["pair_gap"] >= got["gap"]


def test_the_cards_open_to_the_studio_addresses():
    assert [zipper.open_car(n)[0] for n in zipper.CARDS[0]] == [45, 3, 15, 127]
    assert [zipper.open_car(n)[0] for n in zipper.CARDS[1]] == [7, 31, 63, 255]
    r = zipper.run(zipper.plan([[0, 30, 60, 90], [0, 30, 60, 90]]))
    assert [e["tick"] for e in r["events"]] == [0, 57, 114, 171, 228, 285, 342, 399]
    assert zipper.exit_x(18, 36, 156) == 18.0
