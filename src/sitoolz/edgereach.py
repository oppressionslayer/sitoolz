"""edgereach -- push the edge law to higher rungs, compiled.

The run at level v can be read at a single generation, g = 2^v - 1, so one
sequential pass records every rung it passes.  The edge only reads bits
below itself, so a fixed 192-bit window is exact at any generation.

    from sitoolz import edgereach
    edgereach.reach(38)        # rungs 0..38, one pass, numba-compiled

Bare python fallback works but is ~100x slower; use `reach(28)` at most.
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
def _sweep(vmax):
    """Evolve the edge to 2^vmax - 1, recording the run at each 2^v - 1.

    Window is three uint64 words, little-endian: word 0 holds bits 0..63.
    E' = E ^ ((E << 1) | (E << 2)), which never reads above itself.
    """
    w0 = np.uint64(1)
    w1 = np.uint64(0)
    w2 = np.uint64(0)
    out = np.zeros(vmax + 1, dtype=np.int64)

    one = np.uint64(1)
    two = np.uint64(2)
    s62 = np.uint64(62)
    s63 = np.uint64(63)

    g = 0
    nxt = 0                      # next generation of interest: 2^v - 1
    v = 0
    while v <= vmax:
        if g == nxt:
            # trailing ones of the window
            r = 0
            x = w0
            while r < 64 and (x >> np.uint64(r)) & one == one:
                r += 1
            if r == 64:
                x = w1
                while r < 128 and (x >> np.uint64(r - 64)) & one == one:
                    r += 1
            out[v] = r
            v += 1
            if v > vmax:
                break
            nxt = (1 << v) - 1

        a0 = w0 << one
        a1 = (w1 << one) | (w0 >> s63)
        a2 = (w2 << one) | (w1 >> s63)

        b0 = w0 << two
        b1 = (w1 << two) | (w0 >> s62)
        b2 = (w2 << two) | (w1 >> s62)

        w0 = w0 ^ (a0 | b0)
        w1 = w1 ^ (a1 | b1)
        w2 = w2 ^ (a2 | b2)
        g += 1

    return out


def reach(vmax=36):
    """Run lengths for rungs 0..vmax.  Cost is dominated by 2^vmax steps."""
    return [int(x) for x in _sweep(vmax)]


def compare(vmax=36):
    """New rungs against edgelaw.TABLE, and the plateau lengths."""
    from sitoolz import edgelaw
    r = reach(vmax)
    print(" v    run   plateau   status")
    prev = 0
    for v, run in enumerate(r):
        known = edgelaw.TABLE.get(v)
        if known is None:
            status = "NEW"
        elif known == run:
            status = "ok"
        else:
            status = f"MISMATCH (table says {known})"
        print(f"{v:2}   {run:4}   {run - prev:7}   {status}")
        prev = run
    return r
