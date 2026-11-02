"""zipper -- two queues through one gate, the studio's ``zipper-gate-v1``.

Two lanes of cars wait at a gate, each car a packed rulerpack car card with a ready tick. Cars
leave the gate east at one common speed, spaced by a stride of ticks that keeps every adjacent
pair at least the chosen gap apart on the quarter grid:

    clearance_ticks = ceil(120 * (2*halfMax + gap) / speed)       (all on the quarter grid)
    stride          = max(gateTicks, clearance_ticks)

At each admission only the two live heads are read. The time is ``max(nextAllowed, earliest
ready head)``. One eligible head is admitted as it is (``only``); two under the zipper policy
go to the preferred lane (``zipper``), which flips after every admission; two under the arrival
policy go to the earlier head (``oldest``, a tie to the start lane). Once both cars of an adjacent
pair are on the exit their gap is exact: ``speed*(t_follower - t_leader)/120 - h_leader -
h_follower``. In braking mode the gap is the convoy reserve's fitted gap (``convoy``). The
engine's ``wcraft::zipper``, held to the studio's eight receipts.

    from sitoolz import zipper
    r = zipper.run(zipper.plan([[0, 30, 60, 90], [0, 30, 60, 90]]))
    [(e["label"], e["tick"]) for e in r["events"]]   # A1 0, B1 57, A2 114, ... B4 399
"""
import math
from fractions import Fraction

from sitoolz import convoy

HZ = 120
PROFILE = "zipper-gate-v1"
CAR_L = 262144
END_BOUND = 4800
# the studio's eight cars at heading east: address -> longitudinal half-extent (C013 enclosure)
HALF = {45: 2.5, 7: 2.75, 3: 2.5, 31: 2.5, 15: 2.75, 63: 2.25, 127: 2.25, 255: 2.25}
# the cards the studio's cases carry, by lane position: A then B
CARDS = [[9676571645257672083898648, 29020048841878392412045628, 145077591607669446730842464, 1228275790065163840840794500],
         [67701174063992540873097288, 299815669299764290805039216, 609301416990872307920273560, 2466212213788705668701356224]]


def open_car(n, L=CAR_L):
    """A version 1 rulerpack car card back to ``(address, paint)``: ``N = K<<64 | X`` with
    ``X = 2p ^ L*i`` and ``K = 2p ^ L*(i^(i+2))``, so ``K ^ X = L*(i+2)``."""
    x, k = n & ((1 << 64) - 1), n >> 64
    z = k ^ x
    if z % L or z // L < 2:
        raise ValueError("not a car card")
    i = z // L - 2
    v = x ^ (L * i)
    if v & 1:
        raise ValueError("not a car card")
    return i, v >> 1


def plan(ready, cards=None, start_lane=0, policy="zipper", speed=18.0, gate_ticks=36, gap_mode="fixed",
         gap=3.0, delay=24, wet=0.0, margin=2.0, half=None):
    """A merge: two lists of ready ticks (nondecreasing, FIFO), the cars' card words (default the
    studio's), and the gate's settings. ``half`` maps address to half-length for cars outside the
    studio's eight."""
    cards = cards or [CARDS[0][:len(ready[0])], CARDS[1][:len(ready[1])]]
    if len(ready) != 2 or any(len(r) > 4 for r in ready) or not 2 <= len(ready[0]) + len(ready[1]) <= 8:
        raise ValueError("two lanes, up to four cars each, two to eight in all")
    if policy not in ("zipper", "arrival") or gap_mode not in ("fixed", "braking") or start_lane not in (0, 1):
        raise ValueError("policy zipper or arrival, gap mode fixed or braking, start lane 0 or 1")
    if not 12 <= gate_ticks <= 90 or not 5 <= speed <= 30 or not 0 <= gap <= 12 or not 0 <= delay <= 120:
        raise ValueError("gate 12 to 90 ticks, speed 5 to 30, gap 0 to 12, delay 0 to 120")
    h = dict(HALF, **(half or {}))
    lanes, seen = [], set()
    for lane, (ticks, words) in enumerate(zip(ready, cards)):
        if any(b < a for a, b in zip(ticks, ticks[1:])) or any(not 0 <= t <= 240 for t in ticks):
            raise ValueError("ready ticks are 0 to 240 and nondecreasing within a lane")
        cars = []
        for pos, (t, n) in enumerate(zip(ticks, words)):
            address, paint = open_car(n)
            if address in seen:
                raise ValueError("addresses are distinct across the fleet")
            seen.add(address)
            cars.append({"card": n, "address": address, "paint": paint, "ready": int(t), "lane": lane, "position": pos,
                         "label": "AB"[lane] + str(pos + 1), "hx": h[address]})
        lanes.append(cars)
    return {"lanes": lanes, "start_lane": start_lane, "policy": policy, "speed": float(speed), "gate_ticks": int(gate_ticks),
            "gap_mode": gap_mode, "gap": float(gap), "delay": int(delay), "wet": float(wet), "margin": float(margin)}


def run(p):
    """**The admission schedule**: every event with its reason, and the exact gap of each pair."""
    cars = [c for lane in p["lanes"] for c in lane]
    for i, c in enumerate(cars):
        c["actor"] = i
    brake = None
    if p["gap_mode"] == "braking":
        r = convoy.run(convoy.plan("coast", 0, p["delay"], p["gap"], p["margin"], leader_speed=p["speed"],
                                   follower_speed=p["speed"], leader_wet=p["wet"], follower_wet=p["wet"]))
        brake = {"required": r["required"], "fit": r["fitted"], "tail": r["tail"]}
    gap = brake["fit"] if brake else p["gap"]
    half_max = max(c["hx"] for c in cars)
    h4, g4, v4 = round(4 * half_max), round(4 * gap), round(4 * p["speed"])
    clearance = -(-HZ * (2 * h4 + g4) // v4)
    stride = max(p["gate_ticks"], clearance)
    heads = [0, 0]
    next_allowed, preferred = 0, p["start_lane"]
    events, prev = [], None
    while any(heads[l] < len(p["lanes"][l]) for l in (0, 1)):
        live = [(l, p["lanes"][l][heads[l]]) for l in (0, 1) if heads[l] < len(p["lanes"][l])]
        earliest = min(c["ready"] for _, c in live)
        tick = max(next_allowed, earliest)
        eligible = [(l, c) for l, c in live if c["ready"] <= tick]
        if len(eligible) == 1:
            lane, car = eligible[0]
            reason = "only"
        elif p["policy"] == "zipper":
            lane, car = next((l, c) for l, c in eligible if l == preferred)
            reason = "zipper"
        else:
            a, b = eligible
            if a[1]["ready"] == b[1]["ready"]:
                lane, car = next((l, c) for l, c in eligible if l == p["start_lane"])
                reason = "tie"
            else:
                lane, car = min(eligible, key=lambda lc: lc[1]["ready"])
                reason = "oldest"
        pair_gap = None
        if prev is not None:
            pair_gap = Fraction(p["speed"]) * (tick - prev["tick"]) / HZ - Fraction(prev["hx"]) - Fraction(car["hx"])
        events.append({"order": len(events), "actor": car["actor"], "label": car["label"], "address": car["address"],
                       "lane": lane, "position": car["position"], "ready": car["ready"], "tick": tick,
                       "exit": tick + p["gate_ticks"], "wait": tick - car["ready"], "eligible": [c["actor"] for _, c in eligible],
                       "preferred_before": preferred, "slot_open": next_allowed, "arrival_idle": max(0, tick - next_allowed),
                       "reason": reason, "pair_gap": pair_gap})
        heads[lane] += 1
        preferred = 1 - lane
        next_allowed = tick + stride
        prev = {"tick": tick, "hx": car["hx"]}
    end = min(END_BOUND, max(240, events[-1]["exit"] + 120))
    return {"profile": PROFILE, "plan": p, "gap": gap, "half_max": half_max, "clearance_ticks": clearance, "stride": stride,
            "end": end, "brake": brake, "events": events}


def exit_x(speed, exit_tick, now):
    """Where a car is on the exit at tick ``now``: centre zero at its exit, then ``speed`` east."""
    return speed * (now - exit_tick) / HZ
