"""district -- district drive: the city, the car and the wall as one adapter.

The Spindexin studio's ``district-clearance-v1`` profile, composed from three laws this package
already holds: the addressed city (``city``), the one tick law (``motion.step``) and first contact
(``contact.sweep``). A scene is a district scaled to metres once, its 64 footprints and four edges
as solids, and the car's clearance square. One tick proposes the pose, sweeps the square, and at a
contact takes ``max(0, t - 1e-6/distance)`` of the displacement, zeroes speed and velocity and
keeps the yaw. A trace keeps a checkpoint every 32 ticks; the block between two checkpoints is
kept as one packed word of five-bit lanes holding the XOR triangle's left edge, a block card, and
any command is read out of it with at most five whole-word shift/XORs; ``seek`` returns to any tick
from its checkpoint through the commands it reads. The sweep tests only the lots the street layout
puts under the car's box, found by binary search, plus the four edges. The engine's
``wcraft::district``, held to the same six portable cases and the lookup fixture.

    from sitoolz import district
    s = district.scene(45, start="wall")        # (167, 148) facing west, before building 27
    t = district.trace(s, [5] * 600)            # full throttle into the wall
    t["contacts"][0]["tick"], t["contacts"][0]["ids"]   # 208, [2907]
    district.seek(s, t, 208)["pose"] == t["frames"][208]   # True
"""
import math

from sitoolz import city, contact, motion, xortri

HZ, DT = motion.HZ, motion.DT
LIMIT = 600
EVERY = 32
SCALE = 2.0
SKIN = 1e-6
PROFILE = "district-clearance-v1"
EDGE_BASE = (1 << 40) - 4
BLOCK = "spindexin-drive-block"
HANDLING = {"speed": 30.0, "accel": 12.0, "turn": 1.2, "grip": 1.0, "wet": 0.0}


def encode(throttle, steer, brake):
    """The studio's input code: (throttle+1) + 3*(steer+1) + 9*brake."""
    return (throttle + 1) + 3 * (steer + 1) + 9 * int(bool(brake))


def word_at(word, j):
    """Command ``j`` out of a packed edge word: for each set bit ``b`` of ``j`` fold ``word >> 5b``
    in, and the low lane is ``XOR(edge[m] for m submask of j)``, the triangle's closed form.

    >>> word_at(27441, 2)
    11
    """
    b = 1
    while b <= j:
        if j & b:
            word ^= word >> (5 * b)
        b *= 2
    return word & 31


def blank_card():
    return {"format": BLOCK, "version": 1, "length": 0, "N": "0"}


def card_word(card):
    """The card's word, validated: version 1, length 0..32, ``N`` inside its five-bit lanes."""
    if card.get("format") != BLOCK or card.get("version") != 1:
        raise ValueError("open a version 1 drive-block card")
    length = card["length"]
    if not isinstance(length, int) or not 0 <= length <= EVERY:
        raise ValueError("a block holds 0 to 32 commands")
    n = int(card["N"])
    if n < 0 or n >= 1 << (5 * length):
        raise ValueError("keep the word inside its five-bit lanes")
    return n


def extend(card, code):
    """One more command on a card: lane ``j = length`` becomes ``code ^ word_at(N, j)``, the
    triangle transform's own involution, so no earlier command is needed."""
    if not 0 <= int(code) <= 17:
        raise ValueError("an input code is 0 to 17")
    word, j = card_word(card), card["length"]
    if j >= EVERY:
        raise ValueError("a block holds 32 commands")
    e = word_at(word, j) ^ int(code)
    return {"format": BLOCK, "version": 1, "length": j + 1, "N": str(word | (e << (5 * j)))}


def pack_block(codes):
    """Up to 32 commands as one card.

    >>> pack_block([17, 8, 11])["N"]
    '27441'
    """
    card = blank_card()
    for c in codes:
        card = extend(card, c)
    return card


def command_at(card, j):
    """Command ``j`` straight from the card, validated as a code."""
    word = card_word(card)
    if not 0 <= j < card["length"]:
        raise ValueError("the command index is inside the block")
    code = word_at(word, j)
    if code > 17:
        raise ValueError("the card holds a command outside 0 to 17")
    return code


def open_block(card, start=0):
    """The card opened: its edge lanes and every command."""
    word = card_word(card)
    edge = [(word >> (5 * i)) & 31 for i in range(card["length"])]
    recovered = [command_at(card, i) for i in range(card["length"])]
    return {"start": start, "card": dict(card), "edge": edge, "recovered": recovered, "values": list(recovered)}


def card_json(card):
    return '{"format":"%s","version":1,"length":%d,"N":"%s"}' % (BLOCK, card["length"], card["N"])


def card_from_json(text):
    import json
    card = json.loads(text)
    card_word(card)
    return card


def lots_near(s, m):
    """The solids the segment's clearance box can touch: the lot rows and columns under it by
    binary search over the nine street strips, as ``cell = 8*row + col``, then the four edges.
    Returns indices into ``s["solids"]``."""
    def axis(strips, p, d, h):
        low, high = min(p, p + d), max(p, p + d)
        a, b = 0, 8
        while a < b:
            k = (a + b) >> 1
            if strips[k + 1]["lo"] * SCALE + h < low:
                a = k + 1
            else:
                b = k
        first = a
        a, b = 0, 8
        while a < b:
            k = (a + b) >> 1
            if strips[k]["hi"] * SCALE - h <= high:
                a = k + 1
            else:
                b = k
        return first, a
    st = s["district"]["streets"]
    c0, c1 = axis(st["x"], m["x"], m["dx"], m["hx"])
    r0, r1 = axis(st["z"], m["z"], m["dz"], m["hz"])
    return [8 * r + c for r in range(r0, r1) for c in range(c0, c1)] + [64, 65, 66, 67]


def scene(d, edits=None, start="wall", half=2.5, handling=None):
    """District ``d`` in metres with its edited cards, the four edges, the spawn and the clearance.
    ``start`` is ``"avenue"`` (middle avenue, heading north) or ``"wall"`` (before building 27,
    heading west). ``half`` is the car's clearance half-size, a quarter-metre multiple."""
    edits = dict(edits or {})
    c = city.district(d, edits)
    width, depth = c["streets"]["width"] * SCALE, c["streets"]["depth"] * SCALE
    solids = []
    for l in c["lots"]:
        x0, x1, z0, z1 = l["geometry"]["bounds"][:4]
        solids.append({"id": str(l["address"]), "x0": x0 * SCALE, "x1": x1 * SCALE, "z0": z0 * SCALE, "z1": z1 * SCALE})
    for i, (x0, x1, z0, z1) in enumerate([(-16, 0, -16, depth + 16), (width, width + 16, -16, depth + 16),
                                          (0, width, -16, 0), (0, width, depth, depth + 16)]):
        solids.append({"id": str(EDGE_BASE + i), "x0": x0, "x1": x1, "z0": z0, "z1": z1})
    x = c["streets"]["x"][4]["center"] * SCALE
    if start == "wall":
        z, yaw = c["lots"][27]["geometry"]["center"][1] * SCALE, math.pi
    elif start == "avenue":
        z, yaw = c["streets"]["z"][3]["center"] * SCALE, math.pi / 2
    else:
        raise ValueError("start is 'avenue' or 'wall'")
    spawn = motion.initial(x, z, yaw)
    s = {"district": c, "solids": solids, "width": width, "depth": depth, "spawn": spawn, "half": half,
         "handling": dict(handling or HANDLING), "start": start}
    if contact.sweep({"x": x, "z": z, "dx": 0, "dz": 0, "hx": half, "hz": half}, solids)["kind"] == "overlap":
        raise ValueError("this starting clearance overlaps the city")
    return s


def step(s, pose, code, tick=0):
    """One tick under the profile: (next pose, receipt or None)."""
    proposed = motion.step(pose, code, s["handling"])
    m = {"x": pose["x"], "z": pose["z"], "dx": proposed["x"] - pose["x"], "dz": proposed["z"] - pose["z"],
         "hx": s["half"], "hz": s["half"]}
    hit = contact.sweep(m, [s["solids"][i] for i in lots_near(s, m)])
    if hit["kind"] == "clear":
        return proposed, None
    distance = math.hypot(m["dx"], m["dz"])
    used = 0.0 if hit["kind"] == "overlap" else max(0.0, float(hit["fraction"]) - (SKIN / distance if distance else 0.0))
    proposed["x"] = pose["x"] + m["dx"] * used
    proposed["z"] = pose["z"] + m["dz"] * used
    proposed["speed"] = proposed["vx"] = proposed["vz"] = 0.0
    receipt = {"tick": tick, "code": code, "motion": m, "kind": hit["kind"], "fraction": hit["fraction"],
               "ids": [int(i) for i in hit["ids"]], "normals": hit["normals"], "used": used,
               "stop": (proposed["x"], proposed["z"])}
    return proposed, receipt


def empty(s):
    return {"inputs": [], "frames": [dict(s["spawn"])], "contacts": [], "checkpoints": [{"tick": 0, "pose": dict(s["spawn"])}],
            "cards": [blank_card()]}


def append(s, t, code):
    """One more tick on a trace: the frame, any receipt, a checkpoint on every 32nd tick."""
    if len(t["inputs"]) >= LIMIT:
        raise ValueError("this drive has reached 600 ticks")
    if not 0 <= int(code) <= 17:
        raise ValueError("an input code is 0 to 17")
    tick = len(t["inputs"]) + 1
    pose, receipt = step(s, t["frames"][-1], int(code), tick)
    t["inputs"].append(int(code))
    t["cards"][-1] = extend(t["cards"][-1], int(code))
    if tick % EVERY == 0 and tick < LIMIT:
        t["cards"].append(blank_card())
    t["frames"].append(pose)
    if receipt:
        t["contacts"].append(receipt)
    if tick % EVERY == 0:
        t["checkpoints"].append({"tick": tick, "pose": dict(pose)})
    return t


def trace(s, inputs):
    """Drive every code from the spawn."""
    t = empty(s)
    for c in inputs:
        append(s, t, c)
    return t


def block(inputs, tick):
    """The 32-command block holding ``tick``: start, values, the card, its edge lanes, and the
    commands read back from the card."""
    start = tick // EVERY * EVERY
    values = [int(c) for c in inputs[start:start + EVERY]]
    b = open_block(pack_block(values), start)
    b["values"] = values
    return b


def seek(s, t, tick):
    """Return to ``tick`` from its checkpoint, reading each command it needs straight from the
    block's card: at most 31 steps. The pose equals the frame the trace recorded."""
    if not 0 <= tick <= len(t["inputs"]):
        raise ValueError("the replay tick is past the drive")
    mark = t["checkpoints"][tick // EVERY]
    card = t["cards"][tick // EVERY]
    pose, contacts = dict(mark["pose"]), []
    for i in range(tick - mark["tick"]):
        pose, r = step(s, pose, command_at(card, i), mark["tick"] + i + 1)
        if r:
            contacts.append(r)
    return {"pose": pose, "checkpoint": mark["tick"], "steps": tick - mark["tick"], "card": card, "contacts": contacts}


def fork(inputs, tick, code, count):
    """A branch: the original's codes with ``count`` of them replaced by ``code`` from ``tick`` on."""
    if not 0 <= tick <= len(inputs):
        raise ValueError("the fork tick is past the drive")
    if not 1 <= count <= 240:
        raise ValueError("change 1 to 240 ticks")
    if not 0 <= code <= 17:
        raise ValueError("an input code is 0 to 17")
    out = list(inputs)
    for i in range(tick, min(len(out), tick + count)):
        out[i] = code
    return out


def shares_prefix(original, branch, fork_at):
    """Both histories agree before the fork."""
    return fork_at <= len(original) and fork_at <= len(branch) and list(original[:fork_at]) == list(branch[:fork_at])


def first_difference(a, b):
    """The first state where two traces differ, or None."""
    for i, (p, q) in enumerate(zip(a["frames"], b["frames"])):
        if p != q:
            return i
    return None
