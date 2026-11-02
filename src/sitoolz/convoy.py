"""convoy -- the braking reserve between two cars on one lane, the studio's ``convoy-brake-v1``.

Two cars face east on a straight lane under the one tick law (``motion.step``), each with its own
starting speed and wetness, sharing a top speed and an acceleration. Before its brake tick a car
coasts (code 4) or accelerates (code 5); from its brake tick on it holds the brake (code 13). The
follower's brake comes ``delay`` ticks after the leader's. From the two free travel traces the
closure ``C[k] = F[k] - L[k]`` says how much of the starting gap the follower has eaten, and the
reserve is the gap that keeps them apart for good:

    tail     = follower speed at 600 * q / (120 * (1 - q)),  q = 1 - (3 - 1.3*wet) / 120
    required = margin + max(0, max(C[0..600]), C[600] + tail)
    fitted   = ceil((required + 1e-9) * 4) / 4

The tail is every metre the follower still travels after the five visible seconds, braking all the
way. First contact is the pair sweep (``contact.sweep``) of the follower's clearance against the
leader's, tick by tick, with ``offset = gap + hx_leader + hx_follower``. The engine's
``wcraft::convoy``, held to the studio's eight receipts.

    from sitoolz import convoy
    r = convoy.run(convoy.plan(mode="coast", lead_tick=60, delay=60, gap=4, margin=2))
    r["required"], r["fitted"], r["contact"]["tick"]     # 6.108417666118766, 6.25, 235
"""
import math

from sitoolz import contact, motion

HZ, DT = motion.HZ, motion.DT
TICKS = 600
PROFILE = "convoy-brake-v1"
COAST, THROTTLE, BRAKE = 4, 5, 13
# the two studio cars at heading east: leader address 7, follower 45, C013 clearances
LEADER = {"address": 7, "hx": 2.75, "hz": 1.25}
FOLLOWER = {"address": 45, "hx": 2.5, "hz": 1.25}


def plan(mode="coast", lead_tick=60, delay=60, gap=4.0, margin=2.0, top=30.0, accel=12.0,
         leader_speed=18.0, follower_speed=18.0, leader_wet=0.0, follower_wet=0.0,
         leader=LEADER, follower=FOLLOWER):
    """A braking plan: the recipe's fields, with the two cars' clearances."""
    if mode not in ("coast", "throttle"):
        raise ValueError("mode is 'coast' or 'throttle'")
    if not 0 <= lead_tick <= 240 or not 0 <= delay <= 120:
        raise ValueError("the leader brakes at 0 to 240, the follower 0 to 120 ticks later")
    return {"mode": mode, "lead_tick": int(lead_tick), "delay": int(delay), "gap": float(gap), "margin": float(margin),
            "top": float(top), "accel": float(accel),
            "cars": [dict(leader, speed=float(leader_speed), wet=float(leader_wet)),
                     dict(follower, speed=float(follower_speed), wet=float(follower_wet))]}


def travel(car, top, accel, brake_tick, mode, ticks=TICKS):
    """One car's free travel from x = 0: the distances and speeds at every tick, 0 to ``ticks``."""
    h = {"speed": top, "accel": accel, "turn": 1.2, "grip": 1.0, "wet": car["wet"]}
    p = motion.initial(0.0, 0.0, 0.0)
    p["speed"] = p["vx"] = car["speed"]
    xs, vs = [0.0], [car["speed"]]
    free = THROTTLE if mode == "throttle" else COAST
    for k in range(ticks):
        p = motion.step(p, BRAKE if k >= brake_tick else free, h)
        xs.append(p["x"])
        vs.append(p["speed"])
    return xs, vs


def run(plan):
    """**The reserve and the first contact** for a plan: the receipt."""
    lead, foll = plan["cars"]
    lx, lv = travel(lead, plan["top"], plan["accel"], plan["lead_tick"], plan["mode"])
    fx, fv = travel(foll, plan["top"], plan["accel"], plan["lead_tick"] + plan["delay"], plan["mode"])
    closure = [f - l for f, l in zip(fx, lx)]
    b = 3 - 1.3 * foll["wet"]
    q = 1 - b / HZ
    tail = fv[-1] * q / (HZ * (1 - q))
    peak = max(closure)
    closest = closure.index(peak)
    required = plan["margin"] + max(0.0, peak, closure[-1] + tail)
    fitted = math.ceil((required + 1e-9) * 4) / 4
    offset = plan["gap"] + lead["hx"] + foll["hx"]
    first = None
    solid = [{"id": str(lead["address"]), "x0": -lead["hx"], "x1": lead["hx"], "z0": -lead["hz"], "z1": lead["hz"]}]
    for k in range(TICKS):
        m = {"x": fx[k] - lx[k] - offset, "z": 0.0, "dx": (fx[k + 1] - fx[k]) - (lx[k + 1] - lx[k]), "dz": 0.0,
             "hx": foll["hx"], "hz": foll["hz"]}
        hit = contact.sweep(m, solid)
        if hit["kind"] != "clear":
            part = float(hit["fraction"])
            first = {"tick": k, "kind": hit["kind"], "fraction": hit["fraction"], "part": part, "at_tick": k + part,
                     "seconds": (k + part) / HZ, "normals": hit["normals"], "motion": m}
            break
    return {"profile": PROFILE, "plan": plan, "leader": {"x": lx, "speed": lv}, "follower": {"x": fx, "speed": fv},
            "closure": closure, "tail": tail, "peak": peak, "closest_tick": closest, "minimum_gap": plan["gap"] - peak,
            "required": required, "fitted": fitted, "contact": first}


def gap_at(receipt, k):
    """The front-to-back gap at tick ``k``."""
    return receipt["plan"]["gap"] - receipt["closure"][k]


def reserve_at(receipt, k):
    """The reserve still needed at tick ``k``: the larger of the closures still to come and the
    terminal closure plus the tail, less the closure so far, never below zero, plus the margin."""
    c = receipt["closure"]
    ahead = max(max(c[k:]), c[-1] + receipt["tail"])
    return receipt["plan"]["margin"] + max(0.0, ahead - c[k])
