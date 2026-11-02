# Looking at all of it as tables

**Lars Rocha**

Pull every result in `MIXTURE_AND_SPECTRUM.md` and `TRIANGULAR_XOR_SYSTEMS.md` as a
`*_frame()` call. Import, look, sort, filter, export.

```python
from sitoolz import numtools as nt
from sitoolz import driver_lab as dl, parity_lab as pl, txs, edgelaw_deep as ed

dl.driver_frame()          # 128 drivers: structure, prediction, measured verdict
dl.spectrum_frame(3, 1)    # one driver, coordinate by coordinate
dl.shiftset_frame(5)       # every linear shift set, and how each rule scores
pl.parity_frame()          # p, sigma, delta, J, and the two lemmas, per coordinate
txs.clock_frame((1, 3))    # the clock bank: period and family per coordinate
txs.rendezvous_frame()     # which of the 128 drivers re-synchronise
txs.class_frame(5)         # the spectrum equivalence classes
ed.table_frame()           # the edge law, rungs 0..43
```

## pandas optional, and that is the point

```python
from sitoolz.dfkit import backend
backend()        # 'pandas' if it is installed, else 'dfkit.Frame'
```

`frame()` returns a real `pandas.DataFrame` when pandas is importable and a
`Frame` when it is not. Same call, same columns, same code — installing pandas
later upgrades every one of these without an edit.

The pattern in `ntfun` and `ntmc` is `try: import pandas / except: pd = None`,
and then `raise ImportError("pandas required")` at the call site — so on a
machine without pandas, every dataframe function in those two modules is dead
rather than degraded. `dfkit` degrades instead.

`Frame` supports what you actually do to a result table:

```python
df = dl.driver_frame()

df[df["tame"]]                             # boolean column as a mask
df[(df["subcube"]) & (df["ones"] == False)]
df[df["tame"] != df["predicted_tame"]]     # element-wise, column against column
df.sort_values("n_linear", ascending=False).head(10)
df.where(linear=True)
df.value_counts("n_linear")
df.describe()
df.to_markdown()
df.to_csv("drivers.csv")
df.to_pandas()                             # if you install it later
```

## What each table is good for

**`dl.driver_frame()`** — the tame-set result as data. The whole claim is one
line:

```python
len(df[df["tame"] != df["predicted_tame"]])     # 0, over 128 drivers
```

**`dl.shiftset_frame(5)`** — the closed form against the exact rule, per shift
set. The criterion is one line:

```python
all(r["closed_form_ok"] == r["min_eq_gcd"] for r in sf)     # True
```

**`pl.parity_frame()`** — the antipodal lemma is visible as a column: every
row with `doubled=True` has `half_weight=True` and `ratio=0.5000`. Every row
with `fresh=True` has `fresh_ok=True`, which is `δ = J`.

**`txs.rendezvous_frame()`** — the open question as a crosstab:

| G(1,1,1) | drivers | rendezvous | not |
|---|---|---|---|
| False | 64 | 30 | 34 |
| True | 64 | 46 | 18 |

Correlated, nowhere near determinative. The frame also makes plain that
rendezvous and tameness are **independent** — `a+b` is tame and fails it
(8 of 24 coordinates), while plenty of wild drivers hold it.

**`ed.table_frame()`** — the edge law with plateaus computed, and a `source`
column separating the rungs that were already in `edgelaw.TABLE` from the ones
the edge sweep produced:

```python
ed.table_frame()[lambda r: r["plateau"] >= 4]      # the six big plateaus
```

## Files

`sitoolz.dfkit` is standalone and stdlib-only. The `*_frame()` functions import it
lazily, so nothing in the libraries depends on it at import time.
