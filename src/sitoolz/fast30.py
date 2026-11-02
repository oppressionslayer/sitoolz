"""fast30 -- rule 30 centre column at compiled speed.

The evolution packed 64 cells to a word and JIT-compiled, with the working
window tracking the light cone rather than the full ring.

    from sitoolz import fast30
    fast30.centre_column(1_000_000)     # numpy uint8 array, 12 s
    fast30.centre_bits(1_000_000)       # the same thing as a string

Verified bit-identical to ntmc.center_column, and 30x faster than the
big-integer version (12.1 s vs 364 s for a million generations).
"""
import numpy as np

try:
    from numba import njit
except ImportError:                                   # pragma: no cover
    def njit(*a, **k):
        def wrap(f):
            return f
        return wrap if not a else a[0]


@njit(cache=True)
def _evolve(G):
    W = 2 * G + 130
    nw = (W + 63) // 64
    a = np.zeros(nw, dtype=np.uint64)
    b = np.zeros(nw, dtype=np.uint64)
    c = W // 2
    a[c // 64] = np.uint64(1) << np.uint64(c & 63)

    out = np.zeros(G + 1, dtype=np.uint8)
    out[0] = 1
    cw = c // 64
    cb = np.uint64(c & 63)

    for g in range(1, G + 1):
        w0 = (c - g) // 64 - 1                        # only the live cone
        w1 = (c + g) // 64 + 2
        if w0 < 0:
            w0 = 0
        if w1 > nw:
            w1 = nw
        carry = np.uint64(0)
        for j in range(w0, w1):
            cur = a[j]
            nxt = a[j + 1] if j + 1 < nw else np.uint64(0)
            lo = (cur << np.uint64(1)) | carry                  # a[i-1]
            hi = (cur >> np.uint64(1)) | (nxt << np.uint64(63))  # a[i+1]
            b[j] = lo ^ (cur | hi)
            carry = cur >> np.uint64(63)
        for j in range(w0, w1):
            a[j] = b[j]
        out[g] = np.uint8((a[cw] >> cb) & np.uint64(1))
    return out


def centre_column(gens):
    """Centre cell at each generation 0..gens, as a uint8 array."""
    return _evolve(int(gens))


def centre_bits(gens):
    """Centre column as a string of '0' and '1'."""
    return "".join(map(str, _evolve(int(gens))))


def centre_value(gens):
    """Centre column read as a single integer, generation 0 leading."""
    return int(centre_bits(gens), 2)
