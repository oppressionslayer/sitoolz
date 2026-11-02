#!/usr/bin/env python3
"""
square_store.py -- store any object as a DataFrame of its *square index*,
                   and let it re-create itself from that frame.

The idea (Lars's square-index tree):  an integer is not stored, it is
*addressed*.  Split  n = r**2 + rem  with  r = isqrt(n), recurse on both
halves down to small leaves.  The leaves + the tree shape are the seed;
the value regrows bottom-up by  root*root + rem.  The data creates itself.

This wraps that as a pandas DataFrame where every row is one node of the
square-index tree -- so you can *look at* the index, and hand the frame
back to rebuild the original object, exactly.

    from sitoolz import square_store as ss
    df = ss.to_frame({"hello": [1, 2, 3], "lars": "rocha"})
    ss.from_frame(df)                 # -> {'hello': [1, 2, 3], 'lars': 'rocha'}
    df.head()                         # <- the square index, as a table

CLI demo:  python3 square_store.py
"""
import pickle
import pandas as pd

try:                                   # gmpy2 isqrt is much faster on big ints
    from gmpy2 import isqrt as _isqrt, mpz as _mpz
    def _sq(n): return int(_isqrt(_mpz(n)))
except Exception:
    from math import isqrt as _sq      # stdlib fallback


# --------------------------------------------------------------- object <-> int
def _obj_to_int(obj):
    """Serialize any picklable object to a big integer (+ its byte length)."""
    b = pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)
    return int.from_bytes(b, "big"), len(b)

def _int_to_obj(value, nbytes):
    b = int(value).to_bytes(nbytes, "big")   # nbytes restores any leading zeros
    return pickle.loads(b)


# --------------------------------------------------------------- the index tree
def _walk(value, floor):
    """Pre-order square-index tree as a list of node rows."""
    rows = []
    def go(n, depth, kind):
        if n.bit_length() <= floor:          # a leaf: the atom we actually store
            rows.append((depth, kind, 0, int(n)))
            return
        rows.append((depth, kind, 1, -1))    # internal: n = r*r + rem
        r = _sq(n)
        go(r, depth + 1, "root")             # the square root half
        go(n - r * r, depth + 1, "rem")      # the remainder half
    go(int(value), 0, "value")
    return rows

def _rebuild(rows):
    """Regrow the integer from the node rows -- squares and adds, leaves up."""
    pos = [0]
    def build():
        _depth, _kind, internal, leaf = rows[pos[0]]
        pos[0] += 1
        if not internal:
            return leaf
        root = build()
        rem = build()
        return root * root + rem
    return build()


# --------------------------------------------------------------- public API
def to_frame(obj, floor=64):
    """Object -> DataFrame of its square-index tree (self-describing).

    Columns: depth, kind (value/root/rem), internal (1=node,0=leaf), leaf.
    Metadata needed to rebuild rides in df.attrs.
    """
    value, nbytes = _obj_to_int(obj)
    rows = _walk(value, floor)
    df = pd.DataFrame(rows, columns=["depth", "kind", "internal", "leaf"])
    df.attrs.update(nbytes=nbytes, floor=floor, value_bits=int(value).bit_length())
    return df

def from_frame(df):
    """DataFrame -> the original object (it re-creates itself)."""
    rows = list(df[["depth", "kind", "internal", "leaf"]].itertuples(index=False, name=None))
    value = _rebuild(rows)
    return _int_to_obj(value, df.attrs["nbytes"])


def save_frame(df, path):
    """Persist the index frame + its rebuild metadata to parquet."""
    out = df.copy()
    for k, v in df.attrs.items():        # parquet drops .attrs; fold into columns
        out[f"_{k}"] = v
    out.to_parquet(path)
    return path

def load_frame(path):
    df = pd.read_parquet(path)
    meta = {c[1:]: int(df[c].iloc[0]) for c in df.columns if c.startswith("_")}
    df = df[["depth", "kind", "internal", "leaf"]].copy()
    df.attrs.update(meta)
    return df


# --------------------------------------------------------------- demo
if __name__ == "__main__":
    for obj in ["hello, square index!",
                {"lars": "rocha", "primes": [2, 3, 5, 7], "nested": {"a": 1}},
                list(range(200))]:
        df = to_frame(obj, floor=64)
        back = from_frame(df)
        leaves = int((df.internal == 0).sum())
        ok = "OK " if back == obj else "FAIL"
        print(f"{ok} {leaves:4d} leaves, depth {df.depth.max():2d}, "
              f"{df.attrs['value_bits']:6d} value-bits  <-  {type(obj).__name__}: {str(obj)[:38]}")
    # show the index as a table
    df = to_frame("the data isn't stored, it's addressed")
    print("\nthe square index, as a DataFrame:")
    print(df.head(9).to_string(index=False))
