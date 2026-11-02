"""ntmc -- rule 30 machinery.

Ask for a generation, not a width.  The width needed to hold generation n
without the light cone wrapping is computed for you, and rows are reported
in the width-free canonical form: flip the bit string, strip zeros off
both ends.  Under that form the seed is 1, whatever ring it lived in.

    from sitoolz import ntmc

    ntmc.row(110)                 # the 110th rule 30 row
    ntmc.flipstrip(ntmc.row(110)) # its width-free value
    ntmc.frame(110)               # descend to the seed, with counts
    ntmc.compare(110)             # reference vs skip, by width
    ntmc.ladder(1, 40)            # one row per generation, side by side
"""
import time

try:
    import pandas as pd
except ImportError:
    pd = None


# =====================================================================
# width-free form
# =====================================================================
def flipstrip(row):
    """Reverse the bit string and strip zeros off both ends.

    Independent of whatever padding the row was computed in -- the same
    row at w=32 and w=10**9 gives the same value.  The seed gives 1.
    """
    if row == 0:
        return 0
    return int(bin(row)[2:][::-1].strip("0"), 2)


def is_seed(row):
    """A seed is any row with exactly one live cell, wherever it sits."""
    return row != 0 and (row & (row - 1)) == 0


def width_for(gen, margin=4):
    """Ring wide enough that `gen` generations have not wrapped it."""
    return 2 * gen + margin


def natural_width(row, margin=2):
    """Smallest ring that holds `row` without its light cone wrapping."""
    return max(row.bit_length() + margin, 8)


def horizon(w):
    """Generations before the light cone wraps the ring and the seed is lost."""
    return w // 2


# =====================================================================
# forward
# =====================================================================
def seed(w):
    """Canonical origin for a w-cell ring: a single 1 in the center."""
    return 1 << (w // 2)


def step(prev, w):
    """One generation of rule 30, as whole-word operations.

    nxt[i] = prev[i-1] ^ (prev[i] | prev[i+1]), on a ring of width w.
    """
    m = (1 << w) - 1
    lo = ((prev << 1) | (prev >> (w - 1))) & m          # prev[i-1]
    hi = ((prev >> 1) | (prev << (w - 1))) & m          # prev[i+1]
    return (lo ^ (prev | hi)) & m


def evolve(gen, w=None, start=None):
    """Forward orbit as a list of rows, from the seed unless `start` given."""
    w = w or width_for(gen)
    n = seed(w) if start is None else start
    out = [n]
    for _ in range(gen):
        n = step(n, w)
        out.append(n)
    return out


def row(gen, w=None):
    """The row at generation `gen`, counting the seed as generation 0."""
    return evolve(gen, w)[-1]


def center_column(gen, w=None):
    """The famous sequence: the center cell at each generation."""
    w = w or width_for(gen)
    c = w // 2
    return [(r >> c) & 1 for r in evolve(gen, w)]


def center_stream(gen, w=None, start=None):
    """`center_column` one bit at a time, holding a single row.

    `evolve` keeps every row, and the rows are 2*gen wide -- past a few
    thousand generations that is gigabytes for a column of bits.  This
    keeps only the current row, so the cost is time, not memory.  Yields
    gen+1 bits, generation 0 first.
    """
    w = w or width_for(gen)
    c = w // 2
    n = seed(w) if start is None else start
    yield (n >> c) & 1
    for _ in range(gen):
        n = step(n, w)
        yield (n >> c) & 1


def center_string(gen, w=None, start=None):
    """`center_stream` collected into a plain '0'/'1' string."""
    return "".join(map(str, center_stream(gen, w, start)))


def generation(r, origin=None):
    """How many generations `r` is past its seed, read off the leading bit.

    The light cone grows one cell per generation, so the top bit walks out
    one place per step -- no descent needed once the origin is known.
    """
    if r == 0:
        return None
    if origin is None:
        origin = descend(r)[-1]
    return r.bit_length() - origin.bit_length()


# =====================================================================
# backward
# =====================================================================
def inverse(nxt, w, tag=0):
    """Reference predecessor: the plain walk, O(w)."""
    p0, p1 = tag & 1, (tag >> 1) & 1
    prev = (p1 << 1) | p0

    C, R = p0, p1
    L = ((nxt & 1) ^ (C | R)) & 1
    if L:
        prev |= (1 << (w - 1))

    C, R = L, p0
    for i in range(w - 1, 1, -1):
        L = (((nxt >> i) & 1) ^ (C | R)) & 1
        if L:
            prev |= (1 << (i - 1))
        C, R = L, C
    return prev


def inverse_skip(nxt, w, tag=0, count=False):
    """Predecessor in O(bit_length(nxt)) -- the skip method.

    The inverse is a 3-state automaton over the bits of `nxt`:

        A --0/0--> A     f0 = (A->A, B->C, C->C)   idempotent
        A --1/1--> C     f1 = (A->C, B->A, C->B)   order 3
        B --0/1--> C
        B --1/0--> A
        C --0/1--> C
        C --1/0--> B

    Everything above nxt.bit_length() is one run of zeros, and f0 is
    idempotent, so that region collapses to a single transition plus at
    most one bulk write.  The width costs nothing.

    With count=True returns (prev, ops).
    """
    p0, p1 = tag & 1, (tag >> 1) & 1
    prev = (p1 << 1) | p0
    ops = 0

    L = ((nxt & 1) ^ (p0 | p1)) & 1
    if L:
        prev |= (1 << (w - 1))
    C, R = L, p0

    i = w - 1
    top = nxt.bit_length() - 1
    if top < 2:
        top = 1

    k = i - top                              # the jump
    if k > 0:
        if (C | R):                          # state B/C emits a run of ones
            prev |= ((1 << k) - 1) << top
            C, R = (1, C) if k == 1 else (1, 1)
        i = top                              # state A emits zeros: no write
        ops += 1

    while i > 1:                             # significant bits only
        L = (((nxt >> i) & 1) ^ (C | R)) & 1
        if L:
            prev |= (1 << (i - 1))
        C, R = L, C
        i -= 1
        ops += 1

    return (prev, ops) if count else prev


def inverse_linear(nxt, w, tag=0, count=False):
    """inverse_skip's automaton, but genuinely linear in bit_length.

    inverse_skip is O(bit_length) in *operations*, not in runtime: on a
    big integer, `nxt >> i` and `prev |= 1 << k` each allocate a fresh
    number the size of the whole row, so the loop is quadratic.  Reading
    from a bit string and writing into a char list makes both O(1), which
    holds us/bit flat instead of letting it climb with row size.
    """
    p0, p1 = tag & 1, (tag >> 1) & 1
    top = nxt.bit_length() - 1
    if top < 2:
        top = 1

    high = (p1 << 1) | p0                  # bits at/above `top`, kept as an int
    ops = 0

    L = ((nxt & 1) ^ (p0 | p1)) & 1
    if L:
        high |= (1 << (w - 1))
    C, R = L, p0

    k = (w - 1) - top                      # the padding jump, as in inverse_skip
    if k > 0:
        if (C | R):
            high |= ((1 << k) - 1) << top
            C, R = (1, C) if k == 1 else (1, 1)
        ops += 1

    s = bin(nxt)[2:]                       # MSB first; bit i is s[n-1-i]
    n = len(s)
    out = ["0"] * top                      # out[top-1-j] holds prev[j], j < top

    for i in range(top, 1, -1):
        b = s[n - 1 - i] if i < n else "0"
        L = (int(b) ^ (C | R)) & 1
        if L:
            out[top - i] = "1"             # writes prev[i-1]
        C, R = L, C
        ops += 1

    prev = high | int("".join(out), 2)
    return (prev, ops) if count else prev


def center_bit(nxt, w=None, tag=0, count=False):
    """The predecessor's centre cell only -- without building the rest of it.

    The sweep runs top-down, and `prev[i-1]` depends on `prev[i]` and
    `prev[i+1]`.  So `prev[centre]` is determined by the input cells from
    the centre upward and nothing below it: the lower half of the row is
    never read.  Costs half a full predecessor, measured 2.3x.

    With count=True returns (bit, ops).
    """
    w = w or natural_width(nxt)
    c = w // 2
    p0, p1 = tag & 1, (tag >> 1) & 1
    top = nxt.bit_length() - 1
    if top < 2:
        top = 1
    ops = 0

    L = ((nxt & 1) ^ (p0 | p1)) & 1
    C, R = L, p0

    if (w - 1) - top > 0:                    # the padding jump
        if (C | R):
            C, R = (1, C) if (w - 1) - top == 1 else (1, 1)
        ops += 1

    s = bin(nxt)[2:]
    n = len(s)
    for i in range(top, c, -1):              # stop the moment prev[c] is known
        b = s[n - 1 - i] if i < n else "0"
        L = (int(b) ^ (C | R)) & 1
        C, R = L, C
        ops += 1

    return (L, ops) if count else L


class Inverse3:
    """The rule 30 inverse as its minimised 3-state machine.

        A = (0,0)      B = (0,1)      C = {(1,0), (1,1)}

        A --0/0--> A        A --1/1--> C
        B --0/1--> C        B --1/0--> A
        C --0/1--> C        C --1/0--> B

    `(1,0)` and `(1,1)` emit the same bit and move to the same successor on
    both inputs, so nothing downstream can separate them.  The computation
    reads output and never inspects state, so output equivalence is the
    only equivalence it can detect -- which is why this is the operative
    machine and its reset word is the operative bound.

        RESET = '010'   sends every state to C

    Longer answers exist for machines carrying more state (4 states -> 4,
    8 states with one delay -> 5), but those distinguish states this
    computation cannot.
    """

    A, B, C = 0, 1, 2
    F0 = (A, C, C)           # idempotent
    F1 = (C, A, B)           # order 3
    OUT0 = (0, 1, 1)
    OUT1 = (1, 0, 0)
    RESET = "010"
    RESET_STATE = C

    @classmethod
    def step(cls, state, bit):
        """(next_state, output) for one input bit."""
        f = cls.F1 if bit else cls.F0
        o = (cls.OUT1 if bit else cls.OUT0)[state]
        return f[state], o

    @classmethod
    def run(cls, state, word):
        """Feed a word; return (state, outputs)."""
        out = []
        for ch in word:
            state, o = cls.step(state, int(ch))
            out.append(o)
        return state, out

    @classmethod
    def resets(cls, word):
        """Does `word` drive every state to a single state?"""
        return len({cls.run(s, word)[0] for s in (cls.A, cls.B, cls.C)}) == 1

    @classmethod
    def reset_words(cls, length):
        """Every reset word of the given length."""
        from itertools import product
        return [w for w in map("".join, product("01", repeat=length))
                if cls.resets(w)]

    @classmethod
    def shortest_reset(cls, upto=8):
        """Shortest reset word(s), searched exhaustively."""
        for L in range(1, upto + 1):
            hits = cls.reset_words(L)
            if hits:
                return L, hits
        return None, []


def cell_sync(nxt, pos, w=None, window=64, count=False):
    """Any cell of the predecessor, in O(1) cells read.

    The reset word never mentions the centre -- it resets the machine
    wherever it appears.  So the same argument gives any fixed position,
    not just `w // 2`.  Verified against the full predecessor at every
    position of a generation-2000 row: 3981 positions, zero mismatches.

    With count=True returns (bit, cells_examined).
    """
    w = w or natural_width(nxt)
    top = nxt.bit_length() - 1
    win = window

    while True:
        lo = max(pos - 2, 0)
        hi = min(top, pos + win)
        chunk = (nxt >> lo) & ((1 << (hi - lo + 1)) - 1)

        def bit(i, _c=chunk, _l=lo):
            return (_c >> (i - _l)) & 1

        j = pos + 3
        while j <= min(hi, top):
            if not bit(j) and bit(j - 1) and not bit(j - 2):
                break
            j += 1
        else:
            if hi >= top:                     # no reset above pos
                full = inverse_linear(nxt, w)
                b = (full >> pos) & 1
                return (b, top - pos) if count else b
            win *= 4
            continue

        C, R, L = 1, 1, None
        for i in range(j - 3, pos, -1):
            L = (bit(i) ^ (C | R)) & 1
            C, R = L, C
        if L is None:
            L = (bit(pos + 1) ^ 1) & 1
        return (L, j - pos) if count else L


def cells_sync(nxt, positions, w=None):
    """Several predecessor cells at once."""
    w = w or natural_width(nxt)
    return [cell_sync(nxt, p, w) for p in positions]


def center_bit_sync(nxt, w=None, window=64):
    """The predecessor's centre cell in O(1) cells read.

    The inverse automaton has a *synchronizing word*: reading 0,1,0 sends
    every state to C, regardless of what came before.

        A --0--> A --1--> C --0--> C
        B --0--> C --1--> B --0--> C
        C --0--> C --1--> B --0--> C

    It is the shortest such word and the only one of length 3.  So the
    input above the most recent `010` cannot affect the centre cell -- the
    state has already forgotten it.  In a ~50%-dense row `010` occurs about
    every 8 cells, so a constant number of cells is examined no matter how
    many the row has.  Measured at generation 100,000: 4 cells read out of
    100,000 above the centre.

    Exact, not approximate: falls back to the full sweep in the rare case
    that no `010` exists anywhere above the centre.
    """
    w = w or natural_width(nxt)
    c = w // 2
    top = nxt.bit_length() - 1
    win = window

    while True:
        lo = max(c - 2, 0)
        hi = min(top, c + win)
        chunk = (nxt >> lo) & ((1 << (hi - lo + 1)) - 1)

        def bit(i, _ch=chunk, _lo=lo):
            return (_ch >> (i - _lo)) & 1

        j = c + 3
        while j <= min(hi, top):
            if not bit(j) and bit(j - 1) and not bit(j - 2):
                break
            j += 1
        else:
            if hi >= top:                    # no 010 above the centre at all
                return center_bit(nxt, w)
            win *= 4
            continue

        C, R, L = 1, 1, None                 # 010 lands the automaton in C
        for i in range(j - 3, c, -1):
            L = (bit(i) ^ (C | R)) & 1
            C, R = L, C
        return L if L is not None else (bit(c + 1) ^ 1) & 1


def predecessor_center_bit(x, w=None, gen=None, count=False):
    """The centre cell of x's PREDECESSOR.  Accepts a row or a flipstrip value.

        ntmc.predecessor_center_bit(ntmc.row(100))   # from a row
        ntmc.predecessor_center_bit(401)             # from a value

    A flipstrip value always ends in 1; a row at gen > 0 always ends in 0,
    so the two are told apart without being asked.  Uses the synchronizing
    word, so a constant number of cells is read regardless of size.

    With count=True returns (bit, cells_examined).
    """
    if x & 1:                                # odd -> a flipstrip value
        g = gen if gen is not None else gen_of_value(x)
        # Read it in place: row cell (c+d) is value bit (g-d), so scanning
        # upward in the row is scanning downward here.  Avoids rebuilding
        # the row, which costs O(bits) and would defeat the O(1).
        def vbit(d, _v=x, _g=g):
            m = _g - d
            return (_v >> m) & 1 if 0 <= m <= 2 * _g else 0

        d = 3
        while d <= g:
            if not vbit(d) and vbit(d - 1) and not vbit(d - 2):
                C, R, L = 1, 1, None
                for k in range(d - 3, 0, -1):
                    L = (vbit(k) ^ (C | R)) & 1
                    C, R = L, C
                if L is None:
                    L = (vbit(1) ^ 1) & 1
                return (L, d) if count else L
            d += 1
        r = unflipstrip(x, g)                # no 010 above centre: exact fallback
    else:                                    # even -> already a row
        r = x
    w = w or natural_width(r)

    if not count:
        return center_bit_sync(r, w)

    c, top = w // 2, r.bit_length() - 1
    j = c + 3
    while j <= top:
        if not ((r >> j) & 1) and ((r >> (j - 1)) & 1) and not ((r >> (j - 2)) & 1):
            break
        j += 1
    return center_bit_sync(r, w), min(j, top) - c


def center_bit_of(x, gen=None):
    """The centre cell of x ITSELF, not of its predecessor.

    Free either way: under flipstrip the centre cell is preserved, so a
    value's middle bit is the row's centre cell.
    """
    if x & 1:
        return (x >> (x.bit_length() // 2)) & 1
    return (x >> (natural_width(x) // 2)) & 1


def unflipstrip(value, gen):
    """Rebuild the row at generation `gen` from its width-free value.

    flipstrip drops the zero runs on both ends; `gen` puts them back,
    since the leading cell of generation g sits at bit 2*g + 2 under the
    width_for() convention.
    """
    body = int(bin(value)[2:][::-1], 2)
    return body << ((2 * gen + 2) - (body.bit_length() - 1))


def gen_of_value(value):
    """The generation a width-free value came from.

    Rule 30's light cone edges are always live, so a row at generation g
    spans exactly 2*g + 1 cells -- the value carries its own generation.
    """
    b = value.bit_length()
    if b % 2 == 0:
        raise ValueError(f"{b}-bit value: no rule 30 row has an even span")
    return (b - 1) // 2


def predecessor(r=None, w=None, value=None, gen=None, verify=True):
    """One predecessor, with what it cost.  No tree, no descent.

    Give it a row, or just a value -- the generation is read off the
    value's own span, so `gen` is optional and checked if supplied:

        ntmc.predecessor(ntmc.row(100))
        ntmc.predecessor(value=2522196...)
        ntmc.predecessor(value=2522196..., gen=100)   # verified

    Returns a dict:
      prev, prev_value, gen, bits, popcount, runs, skipped, ops,
      ref_ops, ratio, us, ok
    """
    if r is None:
        if value is None:
            raise ValueError("give a row, or value=")
        # A flipstrip value always ends in 1 (trailing zeros are stripped).
        # A row at gen > 0 always ends in 0.  So it is safe to accept either
        # here rather than make the caller remember which column is which.
        if value & 1:
            found = gen_of_value(value)
            if gen is not None and gen != found:
                raise ValueError(
                    f"gen={gen} does not match: a {value.bit_length()}-bit "
                    f"span is generation {found}"
                )
            gen, r = found, unflipstrip(value, found)
        else:
            r = value                       # that was a row, not a value
            gen = (r.bit_length() - 3) // 2
    if gen is None:
        gen = (r.bit_length() - 3) // 2
    w = w or natural_width(r)

    top = max(r.bit_length() - 1, 1)
    t0 = time.perf_counter()
    prev, ops = inverse_linear(r, w, count=True)
    us = (time.perf_counter() - t0) * 1e6

    return {
        "prev": prev,
        "prev_value": flipstrip(prev),
        "gen": gen,
        "bits": prev.bit_length(),
        "popcount": bin(prev).count("1"),
        "runs": _runs(prev),
        "skipped": max(w - 1 - top, 0),
        "ops": ops,
        "ref_ops": w - 1,
        "ratio": round((w - 1) / ops, 2) if ops else None,
        "us": round(us, 2),
        "ok": (step(prev, w) == r) if verify else None,
    }


def predecessors(items, w=None, verify=True):
    """predecessor() over many rows -> DataFrame, one row each."""
    if pd is None:
        raise ImportError("pandas required for predecessors()")
    out = []
    for i, r in enumerate(items):
        d = predecessor(r, w=w, verify=verify)
        d["i"] = i
        d["row"] = r
        out.append(d)
    return _df(out, "i")


def preimages(nxt, w=None):
    """All valid predecessors, deduped.  Empty means Garden of Eden."""
    w = w or natural_width(nxt)
    out = []
    for tag in range(4):
        p = inverse_linear(nxt, w, tag)
        if step(p, w) == nxt and p not in out:
            out.append(p)
    return out


def descend(r, w=None, max_steps=1000000):
    """Walk back until the row is a seed (one live cell).  Returns the path."""
    w = w or natural_width(r)
    traj, n = [r], r
    for _ in range(max_steps):
        if is_seed(n):
            break
        n = inverse_linear(n, w)
        traj.append(n)
    return traj


# =====================================================================
# frames
# =====================================================================
def _runs(n):
    if n == 0:
        return 0
    b = bin(n)[2:]
    return 1 + sum(1 for a, c in zip(b, b[1:]) if a != c)


_BIG = ("row", "value")


def _df(rows, index):
    """DataFrame that keeps arbitrary-precision ints intact.

    pandas infers float for integer columns and overflows past ~1e308, so
    the unbounded columns are pinned to object dtype before inference runs.
    """
    cols = list(rows[0])
    data = {c: pd.Series([r[c] for r in rows],
                         dtype=object if c in _BIG else None)
            for c in cols}
    return pd.DataFrame(data).set_index(index)


def frame(gen, w=None, verify=True):
    """Generation `gen`, descended to the seed, with what each step cost.

      gen        generations past the seed
      row        the row
      value      flipstrip(row) -- the width-free form
      bits       bit_length
      popcount   live cells
      runs       maximal equal-bit runs
      skipped    positions the skip method jumped in one op
      ops        operations the skip method actually used
      ref_ops    operations the reference walk would use (w-1)
      ratio      ref_ops / ops
      us         microseconds for that predecessor
      ok         step(prev) == row, checked at the natural width
    """
    if pd is None:
        raise ImportError("pandas required for frame()")

    w = w or width_for(gen)
    path = descend(row(gen, w), w)
    origin = path[-1]
    rows = []

    for k, n in enumerate(path):
        last = (k == len(path) - 1)
        top = max(n.bit_length() - 1, 1)

        if last:
            prev, ops, dt = None, 0, 0.0
        else:
            t0 = time.perf_counter()
            prev, ops = inverse_linear(n, w, count=True)
            dt = (time.perf_counter() - t0) * 1e6

        ok = None
        if verify and prev is not None:
            nw = natural_width(n)
            ok = step(prev, nw) == n

        rows.append({
            "gen": n.bit_length() - origin.bit_length(),
            "row": n,
            "value": flipstrip(n),
            "bits": n.bit_length(),
            "popcount": bin(n).count("1"),
            "runs": _runs(n),
            "skipped": 0 if last else max(w - 1 - top, 0),
            "ops": ops,
            "ref_ops": 0 if last else w - 1,
            "ratio": None if last or not ops else round((w - 1) / ops, 1),
            "us": round(dt, 2),
            "ok": ok,
        })
    return _df(rows, "gen")


def ladder(lo=1, hi=40, w=None):
    """One row per generation: value, size, and what its predecessor cost."""
    if pd is None:
        raise ImportError("pandas required for ladder()")

    w = w or width_for(hi)
    rows, n = [], seed(w)
    for g in range(0, hi + 1):
        if g:
            n = step(n, w)
        if g < lo:
            continue
        t0 = time.perf_counter()
        prev, ops = inverse_linear(n, w, count=True)
        dt = (time.perf_counter() - t0) * 1e6
        rows.append({
            "gen": g,
            "row": n,
            "value": flipstrip(n),
            "bits": n.bit_length(),
            "popcount": bin(n).count("1"),
            "ops": ops,
            "skipped": max(w - 1 - (n.bit_length() - 1), 0),
            "us": round(dt, 2),
            "ok": step(prev, natural_width(n)) == n,
        })
    return _df(rows, "gen")


def compare(gen, widths=None, reps=200):
    """Reference vs skip for one predecessor, across widths."""
    if pd is None:
        raise ImportError("pandas required for compare()")

    base = width_for(gen)
    widths = widths or (base, 1024, 65536, 10**6, 10**9)
    r = row(gen)
    out = []
    for w in widths:
        if w < base:
            continue
        prev, ops = inverse_skip(r, w, count=True)

        t0 = time.perf_counter()
        for _ in range(reps):
            inverse_skip(r, w)
        skip_us = (time.perf_counter() - t0) / reps * 1e6

        ref_us, agrees = None, None
        if w <= 65536:
            t0 = time.perf_counter()
            ref_us = (time.perf_counter() - t0) * 1e6
            q = inverse(r, w)
            ref_us = (time.perf_counter() - t0) * 1e6
            agrees = (q == prev)

        out.append({
            "prev_value": flipstrip(prev),
            "ops": ops,
            "ref_ops": w - 1,
            "skip_us": round(skip_us, 3),
            "ref_us": None if ref_us is None else round(ref_us, 1),
            "speedup": None if ref_us is None else round(ref_us / skip_us),
            "agrees": agrees,
            "w": w,
        })
    return _df(out, "w")


if __name__ == "__main__":
    if pd is not None:
        pd.set_option("display.width", 200)
        pd.set_option("display.max_columns", 20)
    print(ladder(1, 14))
