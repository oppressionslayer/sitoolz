"""city -- Lars Rocha's addressed city, the law behind the Spindexin district builder.

A district number lays nine street strips on each axis by reading its own hash channels; each
of the 64 lots between them raises a building whose six fields come from the lot's source prime
and the same channels; and every building is one packed integer in rulerpack's ``L = 256, c = 2``
family, the lot card, which opens back to its address and fields. Nothing is stored: the city is
a formula over its address. The same law runs in JavaScript on the site and in Rust in the
engine; the three agree cell for cell.

    from sitoolz import city
    d = city.district(45)               # streets and 64 lots
    d["lots"][27]["fields"]             # floors 4, width 8, depth 9, gable, colour 4, home
    city.open_card(d["lots"][27]["N"])  # (45, 27, fields) back from the packed integer
"""

LOT = 16
MAX_DISTRICT = 4294967295
ROOFS = ["Flat", "Gable", "Terrace", "Spire"]
KINDS = ["Home", "Workshop", "Library", "Landmark"]
COLORS = ["#d3ba92", "#a6bca1", "#acb9c7", "#c79c86", "#c4b1ca", "#e1c887", "#8ea9a1", "#b6c5cd",
          "#c59478", "#d3c9b7", "#92afc4", "#b4b679", "#cb9aa0", "#7eaaa5", "#a49bbd", "#d1aa70"]
_M = 0xFFFFFFFF
_L, _C, _SHIFT = 256, 2, 64


def _mul(a, b):
    return (a * b) & _M


def mix(h):
    """The engine's 32-bit avalanche."""
    h &= _M
    h = _mul(h ^ (h >> 16), 0x7FEB352D)
    h = _mul(h ^ (h >> 15), 0x846CA68B)
    return h ^ (h >> 16)


def at(seed, x, z, stream):
    """The engine's hash field: one 32-bit value per (seed, x, z, stream), stored nowhere."""
    h = seed & _M
    h ^= _mul(x & _M, 0x9E3779B1)
    h ^= _mul(z & _M, 0x85EBCA77)
    h ^= _mul(stream & _M, 0xC2B2AE3D)
    return mix(h)


def stream_of(name):
    """A named stream: FNV-1a over UTF-16 code units, mixed."""
    h = 2166136261
    for unit in name.encode("utf-16-le")[::2] if name.isascii() else _utf16_units(name):
        h ^= unit
        h = _mul(h, 16777619)
    return mix(h)


def _utf16_units(name):
    b = name.encode("utf-16-le")
    return [b[i] | (b[i + 1] << 8) for i in range(0, len(b), 2)]


def channel(d, cell, p, name):
    return at(d, cell, p, stream_of("city." + name))


def layout(d):
    """The streets: nine strips on each axis, 3 or 5 wide from the district's channel, the middle
    one two wider, a sixteen-unit lot between neighbours."""
    def roads(axis):
        cursor, out = 0, []
        for k in range(9):
            width = 3 + 2 * (channel(d, k, 0, axis + "-street") % 2) + (2 if k == 4 else 0)
            out.append({"index": k, "lo": cursor, "hi": cursor + width, "center": cursor + width / 2, "width": width})
            cursor += width + LOT
        return out
    x, z = roads("x"), roads("z")
    return {"x": x, "z": z, "width": x[8]["hi"], "depth": z[8]["hi"]}


_PRIMES = None


def source(d, cell):
    """The lot's source: index 64*(d mod 32) + cell into the first 2048 primes, and that prime."""
    global _PRIMES
    if _PRIMES is None:
        s = bytearray([1]) * 20000
        s[0] = s[1] = 0
        for i in range(2, 142):
            if s[i]:
                s[i * i::i] = bytearray(len(s[i * i::i]))
        _PRIMES = [i for i in range(20000) if s[i]][:2048]
    index = (d % 32) * 64 + cell
    return index, _PRIMES[index]


def generated(d, cell):
    """The building a lot raises on its own."""
    _, p = source(d, cell)
    c = lambda name: channel(d, cell, p, name)
    return {"floors": 1 + c("floors") % 12, "width": 6 + c("width") % 6, "depth": 6 + c("depth") % 6,
            "roof": c("roof") % 4, "color": c("color") % 16, "kind": c("kind") % 4}


def valid(f):
    return (1 <= f["floors"] <= 12 and 4 <= f["width"] <= 11 and 4 <= f["depth"] <= 11
            and 0 <= f["roof"] <= 3 and 0 <= f["color"] <= 15 and 0 <= f["kind"] <= 3)


def payload(f):
    """The 18-bit payload: floors-1 in bits 0-3, width-4 in 4-6, depth-4 in 7-9, roof 10-11,
    colour 12-15, kind 16-17."""
    return (f["floors"] - 1) | ((f["width"] - 4) << 4) | ((f["depth"] - 4) << 7) | (f["roof"] << 10) | (f["color"] << 12) | (f["kind"] << 16)


def from_payload(p):
    f = {"floors": (p & 15) + 1, "width": ((p >> 4) & 7) + 4, "depth": ((p >> 7) & 7) + 4,
         "roof": (p >> 10) & 3, "color": (p >> 12) & 15, "kind": (p >> 16) & 3}
    return f if 0 <= p <= 262143 and valid(f) else None


def address(d, cell):
    return d * 64 + cell


def pack(i, p):
    """rulerpack's pack in the lot-card family: N = K << 64 | X."""
    x = (2 * p) ^ (_L * i)
    k = (2 * p) ^ (_L * (i ^ (i + _C)))
    return (k << _SHIFT) | x


def unpack(n):
    x = n & ((1 << _SHIFT) - 1)
    k = n >> _SHIFT
    z = k ^ x
    if z % _L or z // _L < _C:
        return None
    i = z // _L - _C
    v = x ^ (_L * i)
    if v & 1:
        return None
    p = v >> 1
    if i > MAX_DISTRICT * 64 + 63 or pack(i, p) != n:
        return None
    return i, p


def card(d, cell, f):
    return pack(address(d, cell), payload(f))


def open_card(n):
    """(district, cell, fields) back from a card, or None."""
    u = unpack(n)
    if u is None:
        return None
    i, p = u
    f = from_payload(p)
    return None if f is None else (i // 64, i % 64, f)


def geometry(f, cell, streets):
    """The building's box on its lot and the roof's extra height; the faces as the site draws them."""
    row, col = divmod(cell, 8)
    x = streets["x"][col]["hi"] + LOT / 2
    z = streets["z"][row]["hi"] + LOT / 2
    x0, x1 = x - f["width"] / 2, x + f["width"] / 2
    z0, z1 = z - f["depth"] / 2, z + f["depth"] / 2
    h = f["floors"] * 3
    top = h + [0, 3, 2, 5][f["roof"]]
    return {"bounds": [x0, x1, z0, z1, 0.0, float(top)], "center": [x, z], "wall_height": h,
            "lot": [streets["x"][col]["hi"], streets["x"][col]["hi"] + LOT, streets["z"][row]["hi"], streets["z"][row]["hi"] + LOT]}


def lot(d, cell, edited=None, streets=None):
    streets = streets or layout(d)
    if edited is None:
        f = generated(d, cell)
        n, was = card(d, cell, f), False
    else:
        o = open_card(edited)
        if o is None or o[0] != d or o[1] != cell:
            raise ValueError("keep the card at its own district and lot")
        f, n, was = o[2], edited, True
    return {"cell": cell, "address": address(d, cell), "fields": f, "payload": payload(f), "N": n,
            "source": source(d, cell), "edited": was, "geometry": geometry(f, cell, streets)}


def district(d, edits=None):
    """The whole district: streets and 64 lots, edited cards in place."""
    if not 0 <= d <= MAX_DISTRICT:
        raise ValueError("a district is a whole number from 0 to 4294967295")
    streets = layout(d)
    edits = dict(edits or {})
    return {"district": d, "streets": streets, "lots": [lot(d, c, edits.get(c), streets) for c in range(64)]}


def prefix_row(l):
    """Lars's construction row: prefix XOR over [address, floors, width, depth, roof, color, kind]."""
    f = l["fields"]
    out = [0]
    for v in (l["address"], f["floors"], f["width"], f["depth"], f["roof"], f["color"], f["kind"]):
        out.append(out[-1] ^ v)
    return out
