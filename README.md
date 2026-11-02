# sitoolz

**Build a 19,538-digit prime from M127 in nine rolls. Prove a twin on both sides in under a
second. Hold every prime to 268 million in 15 MB and prove anything to 72 quadrillion from it.
Read rule 30's edge off a table verified to 2⁴⁸.** Lars Rocha's number-theory playground,
every module in one package.

```
pip install sitoolz
```

```python
from sitoolz import showcase, primebase as pb, twin_k2n256, numtools as nt, numtoolz as ntz, xortri

showcase.big_prime()                           # the 19,538-digit prime, built not found
pb.prove.build_up(2**31 - 1, 1024)             # roll M31 to 1024 bits, every rung proven, seconds
twin_k2n256.verify()                           # the 82-digit twin, both sides proven: 0.08 s
pb.is_prime(9_999_999_999_999_917)             # True, by elimination against the net
nt.verdict(1717)                               # ('composite', 'fold-wheel k=16: 17x101', 'free')
ntz.buildsystem(ntz.structured_index(128))     # the XOR Cayley table from its own row 0
xortri.square([2, 6, 8, 11, 12, 5, 4])[:, -1]  # [1, 8, 6, 12, 15, 2]
```

The import name and the distribution are both `sitoolz`. The parts with real compute, the
net, the ladder, the forge, the lanes, the edge law and the triangle, are also a Rust crate,
[`rusttoolz`](https://github.com/oppressionslayer/rusttoolz), held to the same numbers: it rolls M127
to 4,096 bits in four seconds and proves 21 twins of 82 digits in a tenth of one.

## The builds

**[Built, not found](docs/BUILT_NOT_FOUND.md).** Start from a prime you can prove, form
`N = 2·R·q + 1`, and one Pocklington witness proves `N`. Roll nine times from M127 and you hold a
64,903-bit prime whose certificate is the way it came up. The number ships with the package.

**[The twin forge](docs/TWIN_FORGE.md).** `k·2ⁿ ± 1` makes both members of a twin provable,
because `k·2ⁿ` is factored the moment you write it. The forge sieves, screens and proves in
parallel: an 82-digit pair proven both ways, 623-digit and 1,240-digit pairs with the `+1` side
proven, 187 proven Proth primes, and the hunt records to read.

**[The prime net](docs/PRIME_NET.md).** `primebase` keeps every prime to 268,499,993 as half-gaps
and answers from it: the n-th prime, every prime in a window past a trillion, a proof for any
number to 72 quadrillion, its factors, exact counts, and the ladder down to a prime's seed.

**[Rails, lanes, gates and folds](docs/RAILS_AND_LANES.md).** The primes as structure you can
walk. Free lanes settle 88% of every range with no test; the chain-walks and the fold-wheel
factor what slips through; the base-prime bootstrap squares its own reach.

The whole gallery, rule 30's edge law included, is in [docs/README.md](docs/README.md).

## Credit

The ideas are Lars Rocha's, the XOR triangle's closed form among them, solved years before this
package existed. Claude assisted with packaging and the tests; that does not imply ownership.

## The modules

Every module loads on first use, so `import sitoolz` is cheap. The session aliases carry over:
`from sitoolz import numtools as nt, numtoolz as ntz, numtoolsd as ntd, numtoolsm as ntm`.

| group | modules |
|---|---|
| the builds | `showcase` (the numbers), `twin_k2n256` (the proven twin), `twin_forge`, `primebase` (`prove`, `rollball`, `certify`, `counting`, `query`, `neighborhood`) |
| the main library | `numtools` (163 functions: `creator`, `climb`, `rail_walk`, `verdict`, `dashboard`, `roll_prime`, `fold_wheel`, the frames), `numtoolz` (the build system, `structured_index`, the sieves, the Mersenne facts), `numtoolsd` (XOR and Diophantine), `numtoolsm` (perfect machines) |
| rule 30 | `ntfun`, `ntfunextended`, `ntmc`, `ntcomb`, `ntprobe`, `centerrandom`, `fast30`, `edgelaw`, `edgelaw_deep`, `edgeframe`, `edgeproof`, `edgerandom`, `edgereach`, `oeis_crosscheck` |
| drivers and triangular XOR systems | `tri_xor`, `tri_frame`, `txs`, `driver_lab`, `parity_lab`, `dfkit`, `anchors`, `backwalk` |
| the square family | `ntpar`, `ntpar_joint`, `squarewalk`, `square_frame`, `square_rng`, `square_store`, `squareverse`, `antispoof`, `randbattery`, `roguelike`, `shopworld` |
| primes | `qnr2`, `gap_forge`, `prime_check` |
| the city and the car | `city` (the addressed city: streets, lots, cards), `motion` (the one tick law of speed, and the whole flat-ground step), `contact` (first contact along a segment, exact fractions; two moving cars as one sweep), `district` (district drive: the three composed, with checkpoints, XOR-edge replay and branches), `convoy` (the braking reserve between two cars), `zipper` (two queues through one gate); all held to the files the Spindexin site and the Rust engine share |
| the triangle | `xortri` |
| the sequence project | `sequence`: `regen`, `twotables`, `xorsquares`, `rulefinder`, `frame`, `better` (needs sympy) |

`numpy` and `pandas` are required. `gmpy2`, `numba` and `scipy` are used when present
(`pip install "sitoolz[fast]"`) and never required.

## The city, the car and the wall

`sitoolz.city` is Lars's addressed city: a district number lays its streets and raises 64
buildings from their source primes, every building a packed lot card that opens back to its
address. `sitoolz.motion` is the one tick law the site and the engine drive by: accelerate, drag
the new speed, clamp, with the settling speed `throttle·accel·(1/drag − dt)` and the closed form
for a straight run. `sitoolz.contact` is first contact: where a footprint first enters a solid
along a segment, as an exact fraction. `sitoolz.district` composes the three: a district in
metres, the car's clearance square, the stop with a micrometre of skin, a checkpoint every 32
ticks with the block between kept as one packed word of five-bit lanes holding the XOR triangle's
left edge (any command read in at most five shift/XORs), `seek` back to any tick through the
commands it reads, the sweep tested only against the lots the street layout puts under the car,
and branches that share a prefix. `sitoolz.convoy` is the braking reserve between two cars on
one lane: the closure, the tail, the required and fitted gap, and first contact. `sitoolz.zipper`
admits two queues of packed car cards through one gate with exact exit gaps. All six are the same law in JavaScript
on the site, in Rust in the engine and in `rusttoolz`, and the tests hold them to one set of
reference files: the studio's six portable drives replay at zero difference.

## The build system

`xor_prefix(index)` is the prefix-XOR of an index with a leading zero, the anchor row `a`.
`buildsystem(index)` is the table `T[i, j] = a_i ^ a_j`: symmetric, zero on the diagonal, and
the whole of it rebuilds from its own row 0, so the storage is the index alone, itself a
formula (`structured_index`). `recover_index` reads the index back from the row or the table.

## The XOR difference triangle

Start from a row of integers, XOR each neighbour pair, and repeat down to one number:

```
[2, 6, 8, 11, 12, 5, 4]
   [4, 14, 3, 7, 9, 1]
     [10, 13, 4, 14, 8]
       [7, 9, 10, 6]
        [14, 3, 12]
         [13, 15]
           [2]
```

That is Pascal's triangle with XOR for addition, and since a binomial coefficient is odd
exactly when the lower index is a submask of the upper, row `k` has a closed form:

```
row_k[j] = XOR of arr[j + m]  over every m with  m & k == m
```

which is the gasket, `(a & b) == 0`, read as arithmetic. `xortri.row_at(arr, k)` gives any row
straight from the input in `popcount(k)` shift-XORs, the apex is the XOR of the inputs at the
submask positions of `n − 1`, and either edge rebuilds the whole triangle (`from_left_edge`,
`from_right_edge`). `square` is the session's `holder`, and `right_edge` its last column.

## Check

```
pip install -e ".[test]"
pytest -q tests          # every module imports, the doctests, the proofs replay, the docs stay clean
```

## Changes

**0.1.0** (2026-11-01) — the city, motion and contact laws shared with the engine; first packaging of the session modules, primebase and the sequence
project, with the showcase numbers and the docs gallery.
