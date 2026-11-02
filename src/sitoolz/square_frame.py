#!/usr/bin/env python3
"""
square_frame.py -- store anything as a DataFrame of its SQUARE INDEX, and let
                   it re-create itself from that frame.  One self-contained file.

Lars's square-index tree:  an integer is not stored, it is *addressed*.
Split  n = r**2 + rem  with  r = isqrt(n), recurse on both halves down to small
leaves.  The leaves + the tree's shape are a seed; the value regrows bottom-up
by  root*root + rem.  The data creates itself.

This wraps that as a DataFrame you can look at, save compressed, reload, and
rebuild the exact original from -- with three FLAVORS to choose from:

    "node"     one row per tree node -- the whole index, inspectable
    "leaf"     one row per leaf only  -- lean, near the raw size
    "address"  leaves re-expressed as k**2 XOR m, with a carry_free column
               (the ntpar_joint "adds": XOR == + exactly when k**2 & m == 0)

Quick start:
    from sitoolz import square_frame as sf
    df  = sf.to_frame({"lars": "rocha"}, flavor="node")   # object -> square index
    sf.from_frame(df)                                      # index -> object (rebuilds)
    sf.save(df, "thing.sqf")                               # compressed parquet
    sf.from_frame(sf.load("thing.sqf"))                    # reload + rebuild

Files (e.g. RAND's A Million Random Digits):
    df = sf.from_file("AMillionRandomDigits.bin", flavor="leaf")
    sf.to_file(df, "rebuilt.bin")        # byte-for-byte identical

Dependencies are all OPTIONAL: gmpy2 (fast isqrt), pandas (nicer frame),
pyarrow (compressed save).  Core pack/unpack is pure stdlib.

CLI:
    python3 square_frame.py                       # self-test round-trips
    python3 square_frame.py path/to/file [flavor] # frame a real file, all sizes
"""
import json
import os
import pickle

# ------------------------------------------------------------- optional speedups
try:
    from gmpy2 import isqrt as _gsqrt, mpz as _mpz
    def _isqrt(n): return int(_gsqrt(_mpz(n)))
except Exception:
    from math import isqrt as _isqrt

try:
    import pandas as _pd
except Exception:
    _pd = None

# dfkit is Lars's house frame (pandas when present, its own Frame otherwise)
try:
    from sitoolz import dfkit as _dfkit
    def _make_frame(rows):
        return _dfkit.frame(rows)
except Exception:
    def _make_frame(rows):
        if _pd is not None:
            return _pd.DataFrame(rows)
        return rows            # last resort: a list of dicts is still usable


# ------------------------------------------------------------- object <-> integer
def _obj_to_blob(obj):
    """Return (raw_bytes, kind).  Raw bytes stored as-is; anything else pickled."""
    if isinstance(obj, (bytes, bytearray)):
        return bytes(obj), "bytes"
    return pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL), "pickle"

def _blob_to_obj(raw, kind):
    return raw if kind == "bytes" else pickle.loads(raw)


# ------------------------------------------------------------- the square-index tree
def _walk(value, floor):
    """Pre-order walk.  Returns aligned per-node lists: shape, depth, kind, and
    the leaf values (only for shape==0 nodes, in order)."""
    shape, depth, kind, leaves = [], [], [], []
    def go(n, d, k):
        if n.bit_length() <= floor:            # a leaf -- an atom we actually store
            shape.append(0); depth.append(d); kind.append(k); leaves.append(int(n))
            return
        shape.append(1); depth.append(d); kind.append(k)   # internal: n = r*r + rem
        r = _isqrt(n)
        go(r, d + 1, "root")
        go(n - r * r, d + 1, "rem")
    go(int(value), 0, "value")
    return shape, depth, kind, leaves

def _rebuild(shape, leaves):
    """Regrow the integer -- squares and adds, leaves upward."""
    pos = [0]; li = [0]
    def build():
        internal = shape[pos[0]]; pos[0] += 1
        if not internal:
            v = leaves[li[0]]; li[0] += 1
            return v
        root = build(); rem = build()
        return root * root + rem
    return build()


# ------------------------------------------------------------- shape <-> compact hex
def _shape_to_hex(shape):
    n = 0
    for b in shape:
        n = (n << 1) | b
    length = len(shape)
    return f"{length}:" + (n.to_bytes((length + 7) // 8, "big").hex() if length else "")

def _shape_from_hex(s):
    length_s, _, hexs = s.partition(":")
    length = int(length_s)
    n = int.from_bytes(bytes.fromhex(hexs), "big") if hexs else 0
    return [(n >> (length - 1 - i)) & 1 for i in range(length)]


# ------------------------------------------------------------- little bit helpers
def _pop(x):  return bin(x).count("1")
def _v2(x):   return (x & -x).bit_length() - 1 if x else -1


# ------------------------------------------------------------- build a frame
FLAVORS = ("node", "leaf", "address")

def to_frame(obj, flavor="node", floor=64):
    """Object (or raw bytes) -> a square-index DataFrame in the chosen flavor.

    Rebuild metadata rides in df.attrs (also persisted by save()).
    """
    if flavor not in FLAVORS:
        raise ValueError(f"flavor must be one of {FLAVORS}")
    raw, kind_blob = _obj_to_blob(obj)
    value = int.from_bytes(raw, "big")
    shape, depth, kind, leaves = _walk(value, floor)

    if flavor == "node":
        rows, li = [], 0
        for i, (s, d, k) in enumerate(zip(shape, depth, kind)):
            if s:                                   # internal node
                rows.append(dict(node=i, depth=d, kind=k, internal=1, leaf=-1))
            else:
                lv = leaves[li]; li += 1
                rows.append(dict(node=i, depth=d, kind=k, internal=0, leaf=lv))
    elif flavor == "leaf":
        # leaves only, in pre-order; shape rides in attrs
        depth_of_leaf = [d for s, d in zip(shape, depth) if s == 0]
        rows = [dict(leaf_no=i, depth=d, leaf=lv, bits=lv.bit_length(),
                     popcount=_pop(lv), v2=_v2(lv))
                for i, (lv, d) in enumerate(zip(leaves, depth_of_leaf))]
    else:  # address -- the ntpar_joint "adds" flavor
        rows = []
        for i, lv in enumerate(leaves):
            k = _isqrt(lv); ksq = k * k; m = lv ^ ksq   # lv == k**2 XOR m, always
            cf = (ksq & m) == 0                          # XOR == +  (carry-free) ?
            rows.append(dict(leaf_no=i, k=k, m=m, carry_free=int(cf),
                             adds=int(cf), bits=lv.bit_length(), popcount=_pop(lv)))

    df = _make_frame(rows)
    attrs = dict(flavor=flavor, floor=floor, kind_blob=kind_blob,
                 nbytes=len(raw), value_bits=value.bit_length(),
                 shape=_shape_to_hex(shape), n_leaves=len(leaves))
    if hasattr(df, "attrs"):
        df.attrs.update(attrs)
    else:                                    # dfkit.Frame / list fallback carries a sidecar
        try: df._sqf = attrs
        except Exception: pass
    return df


def _attrs_of(df):
    if hasattr(df, "attrs") and df.attrs.get("flavor"):
        return dict(df.attrs)
    if hasattr(df, "_sqf"):
        return dict(df._sqf)
    raise ValueError("this frame carries no square_frame metadata")


def _leaves_of(df, attrs):
    """Recover the ordered leaf integers from any flavor's frame."""
    flavor = attrs["flavor"]
    def col(name):
        c = df[name]
        return list(c) if not hasattr(c, "tolist") else c.tolist()
    if flavor == "node":
        internal = col("internal"); leaf = col("leaf")
        return [lv for s, lv in zip(internal, leaf) if s == 0]
    if flavor == "leaf":
        return [int(x) for x in col("leaf")]
    # address: leaf == k**2 XOR m
    ks = col("k"); ms = col("m")
    return [int(k) * int(k) ^ int(m) for k, m in zip(ks, ms)]


def from_frame(df):
    """DataFrame -> the original object (it re-creates itself)."""
    attrs = _attrs_of(df)
    shape = _shape_from_hex(attrs["shape"])
    leaves = _leaves_of(df, attrs)
    value = _rebuild(shape, leaves)
    raw = int(value).to_bytes(attrs["nbytes"], "big")
    return _blob_to_obj(raw, attrs["kind_blob"])


# ------------------------------------------------------------- file convenience
def from_file(path, flavor="leaf", floor=64):
    with open(path, "rb") as f:
        return to_frame(f.read(), flavor=flavor, floor=floor)

def to_file(df, path):
    obj = from_frame(df)
    if not isinstance(obj, (bytes, bytearray)):
        raise TypeError("frame does not hold raw bytes; use from_frame() instead")
    with open(path, "wb") as f:
        f.write(obj)
    return path


# ------------------------------------------------------------- compressed save/load
def save(df, path, codec="zstd"):
    """Write the frame compressed (parquet), metadata folded into the schema.

    The scaffolding columns (depth/kind/internal) are hugely repetitive, so a
    real codec crushes them -- leaving essentially just the random leaves.
    """
    if _pd is None:
        raise RuntimeError("save() needs pandas + pyarrow")
    import pyarrow as pa, pyarrow.parquet as pq
    pdf = df if isinstance(df, _pd.DataFrame) else _pd.DataFrame(list(df))
    table = pa.Table.from_pandas(pdf, preserve_index=False)
    meta = dict(table.schema.metadata or {})
    meta[b"square_frame"] = json.dumps(_attrs_of(df)).encode()
    table = table.replace_schema_metadata(meta)
    pq.write_table(table, path, compression=codec)
    return path

def load(path):
    if _pd is None:
        raise RuntimeError("load() needs pandas + pyarrow")
    import pyarrow.parquet as pq
    table = pq.read_table(path)
    attrs = json.loads(table.schema.metadata[b"square_frame"])
    df = table.to_pandas()
    df.attrs.update(attrs)
    return df


# ------------------------------------------------------------- tight bit-packed form
def _concat_bits(pairs):
    """Concatenate (value, width) pairs into one integer, divide-and-conquer
    (left-folding would be quadratic once the accumulator is a megabit)."""
    if not pairs:
        return 0, 0
    if len(pairs) == 1:
        return pairs[0]
    mid = len(pairs) // 2
    a, aw = _concat_bits(pairs[:mid])
    b, bw = _concat_bits(pairs[mid:])
    return (a << bw) | b, aw + bw

def _bit_reader(body, total):
    pos = [0]
    def take(k):
        p = pos[0]; pos[0] += k
        if k == 0:
            return 0
        return (body >> (total - (p + k))) & ((1 << k) - 1)
    return take

_WBITS = 7   # a leaf width is 0..floor(<=64) -> fits in 7 bits

def save_compact(obj_or_df, path, floor=64):
    """The 'ship it' form: shape + bit-packed leaves, no int64 padding, no
    scaffolding.  Lands near the raw size even on incompressible data -- the
    leaves' own bits plus a hair of overhead is the information floor.
    Accepts a raw object/bytes, or an existing square_frame (any flavor)."""
    if hasattr(obj_or_df, "__len__") and (hasattr(obj_or_df, "attrs") or hasattr(obj_or_df, "_sqf")):
        attrs = _attrs_of(obj_or_df)
        shape = _shape_from_hex(attrs["shape"])
        leaves = _leaves_of(obj_or_df, attrs)
        kind_blob, nbytes, floor = attrs["kind_blob"], attrs["nbytes"], attrs["floor"]
        value_bits = attrs["value_bits"]
    else:
        raw, kind_blob = _obj_to_blob(obj_or_df)
        value = int.from_bytes(raw, "big")
        shape, _d, _k, leaves = _walk(value, floor)
        nbytes, value_bits = len(raw), value.bit_length()

    pairs = []
    for lv in leaves:
        w = lv.bit_length()
        pairs.append((w, _WBITS)); pairs.append((lv, w))
    body, body_bits = _concat_bits(pairs)
    body_bytes = body.to_bytes((body_bits + 7) // 8, "big") if body_bits else b""
    meta = dict(magic="SQZ1", flavor="compact", floor=floor, kind_blob=kind_blob,
                nbytes=nbytes, value_bits=value_bits, n_leaves=len(leaves),
                wbits=_WBITS, body_bits=body_bits, shape=_shape_to_hex(shape))
    mj = json.dumps(meta).encode()
    with open(path, "wb") as f:
        f.write(b"SQZ1" + len(mj).to_bytes(4, "big") + mj + body_bytes)
    return path

def _load_compact(path):
    with open(path, "rb") as f:
        blob = f.read()
    assert blob[:4] == b"SQZ1"
    mlen = int.from_bytes(blob[4:8], "big")
    meta = json.loads(blob[8:8 + mlen])
    body = int.from_bytes(blob[8 + mlen:], "big")
    take = _bit_reader(body, meta["body_bits"])
    leaves = []
    for _ in range(meta["n_leaves"]):
        w = take(meta["wbits"])
        leaves.append(take(w))
    return meta, leaves

def load_compact(path):
    """SQZ1 -> the original object (recreates itself, tightest form)."""
    meta, leaves = _load_compact(path)
    shape = _shape_from_hex(meta["shape"])
    value = _rebuild(shape, leaves)
    raw = int(value).to_bytes(meta["nbytes"], "big")
    return _blob_to_obj(raw, meta["kind_blob"])


# ------------------------------------------------------------- the structure meter
def meter(data, floor=64, cap=None):
    """How structured is this blob?  Reads the square-index leaves and reports
    the ntpar_joint carry-free rate alongside entropy + gzip as ground truth.

    carry_free_pct = fraction of leaves where  leaf = k**2 XOR m  is secretly
    k**2 + m  (k**2 & m == 0).  Random atoms almost never are; structured data
    yields smaller / sparser leaves that collide less -> the rate climbs.
    """
    import gzip, math, collections
    if not isinstance(data, (bytes, bytearray, str)):
        data, _ = _obj_to_blob(data)          # arbitrary object -> its bytes
    elif isinstance(data, str):
        data = data.encode()
    raw = bytes(data[:cap] if cap else data)
    n = len(raw)
    c = collections.Counter(raw)
    ent = -sum(v / n * math.log2(v / n) for v in c.values()) if n else 0.0
    gz = len(gzip.compress(raw, 9)) / n if n else 0.0
    value = int.from_bytes(raw, "big") if raw else 0
    _shape, _d, _k, leaves = _walk(value, floor)
    cf = tot = 0
    for lv in leaves:
        k = _isqrt(lv); ksq = k * k
        if (ksq & (lv ^ ksq)) == 0:
            cf += 1
        tot += lv.bit_length()
    L = len(leaves) or 1
    return dict(bytes=n, entropy=round(ent, 3), gzip_ratio=round(gz, 3),
                leaves=len(leaves), mean_leaf_bits=round(tot / L, 1),
                carry_free_pct=round(100 * cf / L, 3))


# ------------------------------------------------------------- leaves as random words
def _is_complete(shape, depth, leaves, nbits):
    L = len(leaves)
    if L < 2 or (L & (L - 1)) != 0:                 # leaf count must be a power of two
        return False
    dl = [d for s, d in zip(shape, depth) if s == 0]
    if len(set(dl)) != 1:                            # all leaves at one depth
        return False
    return min(lv.bit_length() for lv in leaves) > nbits   # every leaf wider than a word

def _auto_complete_leaves(value, nbits, floor_cap=1024):
    """Search floor upward for the first COMPLETE tree with all leaves wider than
    nbits -- your tight_pack trick: a higher floor collapses the deep narrow
    leaves into their parents until the tree balances."""
    f = max(64, ((nbits + 8 + 7) // 8) * 8)
    while f <= floor_cap:
        s, d, k, leaves = _walk(value, f)
        if _is_complete(s, d, leaves, nbits):
            return leaves, f
        f += 8
    return None, None

def leaf_words(obj_or_df, nbits=32, floor=64, auto=True, strict=True):
    """The whitened output as fixed-width unsigned words: the low `nbits` of each
    leaf.  Accepts a raw object/bytes or an existing square_frame.

    A clean, PRNG-quality stream requires a COMPLETE tree (all leaves at one depth,
    leaf count a power of two, every leaf wider than a word).  On an incomplete
    tree the deep leaves are narrow (< nbits) and sibling-correlated, which a
    randomness battery detects for ANY input (even os.urandom).

    With `auto=True` (raw input only) this raises the floor automatically until the
    tree is complete -- so any object yields a clean stream with no hand-tuning.
    `strict=True` raises if a clean tree cannot be reached (or if a passed-in frame
    is itself incomplete)."""
    if hasattr(obj_or_df, "__len__") and (hasattr(obj_or_df, "attrs") or hasattr(obj_or_df, "_sqf")):
        leaves = _leaves_of(obj_or_df, _attrs_of(obj_or_df))          # a frame: fixed already
    else:
        raw, _ = _obj_to_blob(obj_or_df)
        value = int.from_bytes(raw, "big")
        if auto:
            got, used = _auto_complete_leaves(value, nbits)
            if got is not None:
                leaves = got
            elif strict:
                raise ValueError(f"could not reach a complete tree up to floor 1024 "
                                 f"for {len(raw)} bytes; try padding the input")
            else:
                _s, _d, _k, leaves = _walk(value, floor)
        else:
            _s, _d, _k, leaves = _walk(value, floor)

    n = len(leaves)
    if strict and (n & (n - 1)) != 0:
        raise ValueError(
            f"incomplete tree ({n:,} leaves, not 2^k): the leaf stream will not be "
            f"clean.  Use auto=True on raw input, or resize, or strict=False to sample anyway.")
    mask = (1 << nbits) - 1
    return [lv & mask for lv in leaves]

def randomness_report(obj_or_df, nbits=32):
    """Run the randbattery on this input's whitened leaves.  One-liner:
        sf.randomness_report(b'ABCD'*100000)   # structured -> is it now random?"""
    from sitoolz import randbattery as rb
    return rb.report(leaf_words(obj_or_df, nbits), nbits=nbits, name="leaves")


# ------------------------------------------------------------- a look at the index
def summary(df):
    """A tiny stats frame about the leaves -- feed it to your plotter."""
    attrs = _attrs_of(df)
    leaves = _leaves_of(df, attrs)
    n = len(leaves)
    pops = [_pop(x) for x in leaves]
    # carry-free rate: how often leaf == k**2 XOR m is secretly k**2 + m
    cf = 0
    for lv in leaves:
        k = _isqrt(lv); ksq = k * k
        if (ksq & (lv ^ ksq)) == 0:
            cf += 1
    rows = [
        dict(stat="flavor",         value=attrs["flavor"]),
        dict(stat="value_bits",     value=attrs["value_bits"]),
        dict(stat="leaves",         value=n),
        dict(stat="tree_depth",     value=len(_shape_from_hex(attrs["shape"])).bit_length()),
        dict(stat="mean_popcount",  value=round(sum(pops) / n, 2) if n else 0),
        dict(stat="carry_free_%",   value=round(100 * cf / n, 2) if n else 0),
    ]
    return _make_frame(rows)


# ------------------------------------------------------------- CLI / self-test
def _selftest():
    trials = ["hello, square index!",
              {"lars": "rocha", "primes": [2, 3, 5, 7], "n": {"a": 1}},
              list(range(300)),
              b"\x00\x00raw bytes with leading zeros\xff\xff"]
    for flavor in FLAVORS:
        print(f"\n=== flavor: {flavor} ===")
        for obj in trials:
            df = to_frame(obj, flavor=flavor)
            back = from_frame(df)
            ok = "OK  " if back == obj else "FAIL"
            rows = len(df) if hasattr(df, "__len__") else "?"
            print(f"  {ok} {rows:>5} rows  <-  {type(obj).__name__}: {str(obj)[:34]}")

def _framefile(path, flavor):
    import time
    raw = open(path, "rb").read()
    print(f"file: {path}\n{len(raw):,} bytes ({len(raw)*8:,} bits)\n")
    for fl in ([flavor] if flavor else FLAVORS):
        t = time.time(); df = to_frame(raw, flavor=fl); tp = time.time() - t
        t = time.time(); back = from_frame(df);          tu = time.time() - t
        rows = len(df)
        line = f"[{fl:7}] {rows:>8,} rows  pack {tp:5.2f}s  unpack {tu:5.2f}s  " \
               f"rebuild {'OK' if back == raw else 'FAIL'}"
        if _pd is not None:
            p = path + f".{fl}.sqf"
            save(df, p, codec="zstd")
            ratio = os.path.getsize(p) / len(raw)
            line += f"  |  parquet {os.path.getsize(p)/1e3:7.0f} KB ({ratio:4.2f}x)"
            assert from_frame(load(p)) == raw            # on-disk form still recreates
        print(line)
    # the tight, ship-it form (flavor-independent)
    if _pd is not None or True:
        import time
        pc = path + ".sqz"
        t = time.time(); save_compact(raw, pc); ts = time.time() - t
        assert load_compact(pc) == raw
        print(f"\n[compact] bit-packed SQZ1: {os.path.getsize(pc)/1e3:.0f} KB "
              f"({os.path.getsize(pc)/len(raw):.2f}x raw) in {ts:.2f}s  "
              f"-- recreates the file exactly")
    print("\n(zstd crushes the repetitive tree scaffolding; the leaves are the "
          "irreducible part -- random data cannot compress below itself.\n"
          " the visual frames are for LOOKING; SQZ1 is for SHIPPING.)")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        _framefile(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    else:
        _selftest()
        print("\nthe square index, as a frame (node flavor):")
        d = to_frame("the data isn't stored, it's addressed")
        print(d.head(8) if hasattr(d, "head") else d[:8])
