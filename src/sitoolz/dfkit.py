"""dfkit -- dataframes for the rule 30 / driver-space libraries, pandas optional.

Every `*_frame()` function in `driver_lab`, `parity_lab`, `txs` and
`edgelaw_deep` returns whatever this module's `frame()` gives it:

    pandas installed      ->  a real pandas.DataFrame, nothing else changes
    pandas not installed  ->  a Frame: the same table, printable and sliceable

So the same code works either way, and installing pandas later upgrades every
one of them without an edit.  That is the point -- the existing house pattern
(`try: import pandas / except: pd = None`, then `raise ImportError` at the call
site) leaves those functions dead on a machine without pandas, which is this
one.

    from sitoolz.dfkit import frame

    f = frame([{"j": 0, "p": 0}, {"j": 1, "p": 1}])
    f.sort_values("p", ascending=False)
    f[f["p"] > 0]                      # boolean-ish mask, or f.filter(fn)
    f.to_markdown()  /  f.to_csv("out.csv")  /  f.describe()
"""
try:
    import pandas as pd
except ImportError:                                   # pragma: no cover
    pd = None

HAVE_PANDAS = pd is not None


# ------------------------------------------------------------------ the table
class Frame:
    """A minimal tabular container: ordered columns over a list of row dicts.

    Deliberately not a pandas clone.  It covers looking at a result, sorting
    it, slicing a column out and writing it somewhere -- which is what these
    tables are for.
    """

    def __init__(self, rows, columns=None, index=None):
        self.rows = [dict(r) for r in rows]
        if columns is None:
            columns = []
            for r in self.rows:
                for k in r:
                    if k not in columns:
                        columns.append(k)
        self.columns = list(columns)
        self.index = index
        if index and index in self.columns:
            self.columns.remove(index)
            self.columns.insert(0, index)

    # ---------------------------------------------------------- basic protocol
    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    @property
    def shape(self):
        return (len(self.rows), len(self.columns))

    def __getitem__(self, key):
        if isinstance(key, str):
            return Column(self, key)
        if isinstance(key, (list, tuple)):
            return Frame([{k: r.get(k) for k in key} for r in self.rows], list(key))
        if callable(key):
            return self.filter(key)
        if isinstance(key, slice):
            return Frame(self.rows[key], self.columns, self.index)
        if isinstance(key, Mask):
            return Frame([r for r, keep in zip(self.rows, key.bits) if keep],
                         self.columns, self.index)
        if isinstance(key, Column):
            return Frame([r for r, keep in zip(self.rows, key.values) if keep],
                         self.columns, self.index)
        raise TypeError(f"cannot index a Frame with {type(key).__name__}")

    # ------------------------------------------------------------- reshaping
    def head(self, n=10):
        return Frame(self.rows[:n], self.columns, self.index)

    def tail(self, n=10):
        return Frame(self.rows[-n:], self.columns, self.index)

    def filter(self, fn):
        """Rows where fn(row) is true."""
        return Frame([r for r in self.rows if fn(r)], self.columns, self.index)

    def where(self, **kw):
        """Rows matching every column=value given."""
        return self.filter(lambda r: all(r.get(k) == v for k, v in kw.items()))

    def sort_values(self, by, ascending=True):
        keys = [by] if isinstance(by, str) else list(by)

        def k(r):
            return tuple(_sortable(r.get(c)) for c in keys)
        return Frame(sorted(self.rows, key=k, reverse=not ascending),
                     self.columns, self.index)

    def drop(self, *cols):
        keep = [c for c in self.columns if c not in cols]
        return Frame([{k: r.get(k) for k in keep} for r in self.rows], keep,
                     self.index if self.index in keep else None)

    def assign(self, **kw):
        """New columns from callables of the row, or constants."""
        rows = []
        for r in self.rows:
            r = dict(r)
            for name, v in kw.items():
                r[name] = v(r) if callable(v) else v
            rows.append(r)
        return Frame(rows, self.columns + [k for k in kw if k not in self.columns],
                     self.index)

    def set_index(self, col):
        return Frame(self.rows, self.columns, col)

    # --------------------------------------------------------------- summaries
    def value_counts(self, col):
        counts = {}
        for r in self.rows:
            counts[r.get(col)] = counts.get(r.get(col), 0) + 1
        return Frame([{col: k, "count": v} for k, v in
                      sorted(counts.items(), key=lambda kv: -kv[1])])

    def groupby(self, col):
        """dict of value -> Frame."""
        out = {}
        for r in self.rows:
            out.setdefault(r.get(col), []).append(r)
        return {k: Frame(v, self.columns, self.index) for k, v in out.items()}

    def describe(self):
        out = []
        for c in self.columns:
            vals = [r.get(c) for r in self.rows]
            nums = [v for v in vals if isinstance(v, (int, float))
                    and not isinstance(v, bool)]
            if not nums:
                continue
            out.append({"column": c, "count": len(nums), "min": min(nums),
                        "max": max(nums),
                        "mean": round(sum(nums) / len(nums), 4),
                        "nulls": sum(1 for v in vals if v is None)})
        return Frame(out)

    # ------------------------------------------------------------------ output
    def to_dicts(self):
        return [dict(r) for r in self.rows]

    def to_csv(self, path=None, sep=","):
        lines = [sep.join(self.columns)]
        for r in self.rows:
            lines.append(sep.join("" if r.get(c) is None else str(r.get(c))
                                  for c in self.columns))
        text = "\n".join(lines) + "\n"
        if path:
            with open(path, "w") as fh:
                fh.write(text)
            return path
        return text

    def to_markdown(self):
        w = self._widths()
        head = "| " + " | ".join(c.ljust(w[c]) for c in self.columns) + " |"
        sep = "|" + "|".join("-" * (w[c] + 2) for c in self.columns) + "|"
        body = ["| " + " | ".join(_fmt(r.get(c)).ljust(w[c])
                                  for c in self.columns) + " |" for r in self.rows]
        return "\n".join([head, sep] + body)

    def to_pandas(self):
        if pd is None:
            raise ImportError("pandas is not installed")
        df = pd.DataFrame(self.rows, columns=self.columns)
        return df.set_index(self.index) if self.index else df

    # -------------------------------------------------------------- rendering
    def _widths(self):
        w = {}
        for c in self.columns:
            w[c] = max(len(str(c)),
                       max((len(_fmt(r.get(c))) for r in self.rows), default=0))
        return w

    def __repr__(self):
        if not self.rows:
            return f"Frame: 0 rows x {len(self.columns)} columns {self.columns}"
        show, elided = self.rows, 0
        if len(self.rows) > 22:
            show, elided = self.rows[:14] + self.rows[-5:], len(self.rows) - 19
        w = self._widths()
        right = {c: all(isinstance(r.get(c), (int, float)) or r.get(c) is None
                        for r in self.rows) for c in self.columns}
        out = ["  ".join((str(c).rjust(w[c]) if right[c] else str(c).ljust(w[c]))
                         for c in self.columns)]
        out.append("  ".join("-" * w[c] for c in self.columns))
        for i, r in enumerate(show):
            if elided and i == 14:
                out.append(f"  ... {elided} more rows ...")
            out.append("  ".join(
                (_fmt(r.get(c)).rjust(w[c]) if right[c] else _fmt(r.get(c)).ljust(w[c]))
                for c in self.columns))
        out.append(f"\n[{len(self.rows)} rows x {len(self.columns)} columns]")
        return "\n".join(out)


class Column:
    """One column, comparable so that `f[f["p"] > 3]` works."""

    def __init__(self, frame, name):
        self.frame, self.name = frame, name
        self.values = [r.get(name) for r in frame.rows]

    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.values)

    def __getitem__(self, i):
        return self.values[i]

    def __repr__(self):
        return f"Column({self.name!r}, {self.values!r})"

    def tolist(self):
        return list(self.values)

    def _cmp(self, other, op):
        """Element-wise against another column or list, broadcast against a scalar."""
        if isinstance(other, Column):
            other = other.values
        if isinstance(other, (list, tuple)):
            if len(other) != len(self.values):
                raise ValueError(f"length mismatch: {len(self.values)} vs {len(other)}")
            return Mask([op(a, b) for a, b in zip(self.values, other)])
        return Mask([op(v, other) for v in self.values])

    @staticmethod
    def _ord(op):
        """Ordering comparisons are false wherever either side is missing."""
        return lambda a, b: a is not None and b is not None and op(a, b)

    def __gt__(self, o): return self._cmp(o, self._ord(lambda a, b: a > b))
    def __ge__(self, o): return self._cmp(o, self._ord(lambda a, b: a >= b))
    def __lt__(self, o): return self._cmp(o, self._ord(lambda a, b: a < b))
    def __le__(self, o): return self._cmp(o, self._ord(lambda a, b: a <= b))
    def __eq__(self, o): return self._cmp(o, lambda a, b: a == b)
    def __ne__(self, o): return self._cmp(o, lambda a, b: a != b)

    def _mask(self):
        return Mask([bool(v) for v in self.values])

    def __and__(self, o): return self._mask() & _as_mask(o)
    def __or__(self, o): return self._mask() | _as_mask(o)
    def __rand__(self, o): return _as_mask(o) & self._mask()
    def __ror__(self, o): return _as_mask(o) | self._mask()
    def __invert__(self): return ~self._mask()

    def isin(self, options):
        opts = set(options)
        return Mask([v in opts for v in self.values])

    def unique(self):
        seen, out = set(), []
        for v in self.values:
            if v not in seen:
                seen.add(v)
                out.append(v)
        return out


class Mask:
    """A boolean selector, combinable with & and |."""

    def __init__(self, bits):
        self.bits = list(bits)

    def __and__(self, o):
        return Mask([a and b for a, b in zip(self.bits, _as_mask(o).bits)])

    def __or__(self, o):
        return Mask([a or b for a, b in zip(self.bits, _as_mask(o).bits)])
    def __invert__(self): return Mask([not a for a in self.bits])
    def __repr__(self): return f"Mask({sum(self.bits)} of {len(self.bits)} selected)"


def _as_mask(o):
    """A Mask, a boolean Column, or a list of booleans -- all usable as one."""
    if isinstance(o, Mask):
        return o
    if isinstance(o, Column):
        return o._mask()
    if isinstance(o, (list, tuple)):
        return Mask([bool(v) for v in o])
    raise TypeError(f"cannot combine a mask with {type(o).__name__}")


def _fmt(v):
    if v is None:
        return "·"
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


def _sortable(v):
    """Sort keys that never raise on mixed types or None."""
    if v is None:
        return (0, 0, "")
    if isinstance(v, bool):
        return (1, int(v), "")
    if isinstance(v, (int, float)):
        return (1, v, "")
    return (2, 0, str(v))


# ------------------------------------------------------------------- factory
def frame(rows, columns=None, index=None):
    """A pandas.DataFrame if pandas is importable, otherwise a Frame."""
    if pd is not None:
        df = pd.DataFrame(list(rows), columns=columns)
        return df.set_index(index) if index else df
    return Frame(rows, columns, index)


def backend():
    return "pandas" if HAVE_PANDAS else "dfkit.Frame"
