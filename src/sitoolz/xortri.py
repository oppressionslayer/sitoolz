"""xortri — the XOR difference triangle, and the square it packs into.

Start from a row of integers and XOR each neighbour pair; do it again on the result, down to
one number. That is Pascal's triangle with XOR for addition, and because a binomial
coefficient is odd exactly when the lower index is a submask of the upper (Lucas), row `k`
has a closed form:

    row_k[j] = XOR of arr[j + m]  over every m with  m & k == m

which is the gasket, `(a & b) == 0`, read as arithmetic. Lars Rocha solved this closed form
years before this package; it is written here as he had it. Two things follow. Any single row,
the apex included, comes straight from the input in `popcount(k)` shift-XORs (`row_at`),
without building the rows above it. And either edge of the triangle rebuilds the whole
thing (`from_left_edge`), which is the rank-1 property of the build system again.

The session this came from (`holder[:, -1]` of `[2, 6, 8, 11, 12, 5, 4]` is
`[1, 8, 6, 12, 15, 2]`) is `square` and `right_edge`.
"""

import numpy as np

__all__ = ["triangle", "triangle_fast", "row_at", "apex", "square", "left_edge", "right_edge",
           "from_left_edge", "from_right_edge", "feeds", "taps", "explain"]


def _arr(a):
    return np.asarray(a, dtype=np.int64)


def triangle(arr):
    """The rows under `arr`, each the XOR of neighbours in the one above, down to one number.

    The session's loop, as it was written. Returns the list of rows, not including `arr`.

    >>> [r.tolist() for r in triangle([2, 6, 8, 11, 12, 5, 4])]
    [[4, 14, 3, 7, 9, 1], [10, 13, 4, 14, 8], [7, 9, 10, 6], [14, 3, 12], [13, 15], [2]]
    """
    a = _arr(arr)
    rows = []
    while a.size > 1:
        a = a[:-1] ^ a[1:]
        rows.append(a)
    return rows


def row_at(arr, k):
    """Row `k` of the triangle (row 0 is `arr` itself) in `popcount(k)` shift-XORs.

    Row 2^b is `a[:-2^b] ^ a[2^b:]` in one step, and the shift-XOR operators commute, so row
    `k` is the product over the set bits of `k`. This never touches the rows in between, which
    is what makes the apex of a long row cheap.

    >>> row_at([2, 6, 8, 11, 12, 5, 4], 3).tolist()
    [7, 9, 10, 6]
    >>> row_at([2, 6, 8, 11, 12, 5, 4], 4).tolist()
    [14, 3, 12]
    >>> row_at([2, 6, 8, 11, 12, 5, 4], 0).tolist()
    [2, 6, 8, 11, 12, 5, 4]
    """
    a = _arr(arr)
    k = int(k)
    if k < 0 or k >= a.size:
        raise ValueError(f"row {k} of a triangle with {a.size} rows")
    b = 0
    while k:
        if k & 1:
            s = 1 << b
            a = a[:-s] ^ a[s:]
        k >>= 1
        b += 1
    return a


def triangle_fast(arr):
    """The same rows as `triangle`, each straight from the input by `row_at`.

    Rows cost `n · popcount(k)` instead of chaining through every row above, and any subset of
    rows can be asked for alone. Same answer, bit for bit.

    >>> a = [2, 6, 8, 11, 12, 5, 4]
    >>> all((x == y).all() for x, y in zip(triangle(a), triangle_fast(a)))
    True
    """
    a = _arr(arr)
    return [row_at(a, k) for k in range(1, a.size)]


def apex(arr):
    """The bottom of the triangle: XOR of `arr[m]` over the submasks `m` of `n - 1`.

    For seven inputs that is the even positions, 2 ^ 8 ^ 12 ^ 4.

    >>> apex([2, 6, 8, 11, 12, 5, 4])
    2
    """
    a = _arr(arr)
    return int(row_at(a, a.size - 1)[0])


def feeds(k, n):
    """Which positions of an `n`-long input feed `row_k[0]`: the submasks of `k`.

    `feeds(k, n)[j]` is True when `j & k == j` and `j < n`, which is the gasket row `k`.

    >>> feeds(6, 8).astype(int).tolist()
    [1, 0, 1, 0, 1, 0, 1, 0]
    """
    j = np.arange(int(n))
    return (j & int(k)) == j


def square(arr):
    """The session's `holder`: the rows packed into the upper triangle of an (n-1)×(n-1)
    array, row `k` starting at column `k`. Everything under the diagonal is zero.

    >>> square([2, 6, 8, 11, 12, 5, 4])[:, -1].tolist()
    [1, 8, 6, 12, 15, 2]
    """
    rows = triangle(arr)
    m = len(rows)
    out = np.zeros((m, m), dtype=np.int64)
    for k, r in enumerate(rows):
        out[k, k:] = r
    return out


def left_edge(arr):
    """The first entry of every row, `arr[0]` first: the triangle's left side.

    >>> left_edge([2, 6, 8, 11, 12, 5, 4]).tolist()
    [2, 4, 10, 7, 14, 13, 2]
    """
    a = _arr(arr)
    return np.array([int(a[0])] + [int(r[0]) for r in triangle(a)], dtype=np.int64)


def right_edge(arr):
    """The last entry of every row, `arr[-1]` first: the triangle's right side, which is the
    session's `holder[:, -1]` with the input's own last value in front.

    >>> right_edge([2, 6, 8, 11, 12, 5, 4]).tolist()
    [4, 1, 8, 6, 12, 15, 2]
    """
    a = _arr(arr)
    return np.array([int(a[-1])] + [int(r[-1]) for r in triangle(a)], dtype=np.int64)


def from_left_edge(edge):
    """The input back from its left edge: `row_k[j+1] = row_k[j] ^ row_{k+1}[j]`, built up
    from the apex. Exact inverse of `left_edge`, so the edge alone is the whole triangle.

    >>> from_left_edge([2, 4, 10, 7, 14, 13, 2]).tolist()
    [2, 6, 8, 11, 12, 5, 4]
    """
    e = _arr(edge)
    row = e[-1:].copy()
    for k in range(e.size - 2, -1, -1):
        nxt = np.empty(row.size + 1, dtype=np.int64)
        nxt[0] = e[k]
        nxt[1:] = row  # placeholder, filled by the scan below
        # row_k[j+1] = row_k[j] ^ row_{k+1}[j]: a prefix XOR of the row below, seeded by the edge
        nxt[1:] = e[k] ^ np.bitwise_xor.accumulate(row)
        row = nxt
    return row


def from_right_edge(edge):
    """The input back from its right edge, by mirror: the triangle of a reversed input is the
    reversed triangle, so the right edge is the left edge of `arr[::-1]`.

    >>> from_right_edge([4, 1, 8, 6, 12, 15, 2]).tolist()
    [2, 6, 8, 11, 12, 5, 4]
    """
    return from_left_edge(edge)[::-1]


def taps(k):
    """The submasks of ``k``, in order: the offsets row ``k`` reads from the input. There are
    ``2 ** popcount(k)`` of them, and ``k`` itself is the last.

    >>> taps(5)
    [0, 1, 4, 5]
    >>> taps(4)
    [0, 4]
    """
    k = int(k)
    out, m = [], 0
    while True:
        out.append(m)
        if m == k:
            return out
        m = (m - k) & k   # the next submask: subtract 1 inside k's bits


def explain(arr, k=None, j=0):
    """**Show the XOR.** Cell ``j`` of row ``k`` as the taps it is made of: each submask ``m`` of
    ``k``, its bits, the input position ``j + m`` it reaches, the value there, and the running
    XOR down the column, ending at the cell. ``k`` defaults to the apex, ``len(arr) - 1``.

    >>> explain([1, 11, 9, 11, 3, 15])
       m bits  index  value  running
    0  0  000      0      1        1
    1  1  001      1     11       10
    2  4  100      4      3        9
    3  5  101      5     15        6
    """
    import pandas as pd
    a = _arr(arr)
    n = a.size
    k = n - 1 if k is None else int(k)
    if not 0 <= k < n or not 0 <= j < n - k:
        raise ValueError("row k needs 0 <= k < len(arr) and 0 <= j < len(arr) - k")
    width = max(1, k.bit_length())
    rows, run = [], 0
    for m in taps(k):
        run ^= int(a[j + m])
        rows.append({"m": m, "bits": format(m, f"0{width}b"), "index": j + m, "value": int(a[j + m]), "running": run})
    return pd.DataFrame(rows)
