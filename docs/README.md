# The sitoolz gallery

What Lars Rocha's number-theory playground built, and how to run each piece yourself.

## The builds

- **[Built, not found](BUILT_NOT_FOUND.md).** A 19,538-digit prime rolled up from M127 in nine
  rolls, its proof carried up with it. `showcase.big_prime()` hands you the number;
  `primebase.prove.build_up` rolls your own.
- **[The twin forge](TWIN_FORGE.md).** Twin primes at any size: an 82-digit pair proven on both
  sides in under a second, 623-digit and 1,240-digit pairs with the `+1` side proven, 187
  proven Proth primes along the way. `twin_k2n256.verify()`, `twin_forge.hunt(n)`.
- **[The prime net](PRIME_NET.md).** Every prime to 268,499,993 in 15 MB, proving, factoring
  and generating to 72 quadrillion. `primebase.is_prime`, `primes`, `pi`, `descend`.
- **[Rails, lanes, gates and folds](RAILS_AND_LANES.md).** The primes as structure you can walk:
  88% of every range settled by free logic, the slip band factored by walks, the base-prime
  bootstrap that squares its own reach. `numtools.verdict`, `walk_rail`, `fold_sharp`.

## The prime kit

- **[The prime toolkit](PRIME_TOOLKIT.md).** Every tool in `numtools` and the structures it
  found: the mod-4 fork, the two-square theorem walked, the fold multiplier's octaves.
- **[The cookbook](COOKBOOK.md).** Ten recipes with the output they print.
- **[Mersenne beacons](MERSENNE_BEACONS.md).** The Mersenne primes, perfect numbers and
  Lucas–Lehmer, as rulers and beacons.
- **[Squaring and factors](SQUARING_AND_FACTORS.md).** Squaring keeps a number's factors; XOR
  does not, and that is the whole difference.
- **[The bits of powers](THE_BITS_OF_POWERS.md).** The bits a square cannot have, used as gates.
- **[Data frames](DATAFRAMES.md).** The `*_frame()` views, one table per structure.

## Rule 30 and the drivers

- **[The edge law](EDGE_LAW.md).** Rule 30's edge read off a lookup table, with the period
  spectrum, the proofs, the prior art, and what is still open.
- **[The odometer](THE_ODOMETER.md).** From measurement to law: how the edge was clocked.
- **[The centre law](CENTER_LAW.md)**, **[column 53209](COLUMN_53209.md)**,
  **[the spectrum breakthrough](THE_SPECTRUM_BREAKTHROUGH.md)**,
  **[the v₂ law proof](PROOF_V2_LAW.md)**, **[the proofs](PROOFS.md)**.
- **[An invitation](AN_INVITATION.md).** The doubling, and the open problems around it.
- **[Driver space](DRIVER_SPACE.md).** The two-input and three-input driver families, what was
  measured, and what held: the spectrum is a period ladder.
- **[Triangular XOR systems](TRIANGULAR_XOR_SYSTEMS.md).** The register class and its drivers.

## The sequence project

- **[The reveal](sequence/reveal.md).** How the puzzle sequences work: the product form, the
  parity bridge, and `n = m² ⊕ 5` as a bijection. The tables are beside it as `.dat` files.

## Running any of it

```
pip install sitoolz
```

```python
from sitoolz import numtools as nt, numtoolz as ntz, primebase as pb, showcase, xortri
```

Every module loads on first use. The session aliases carry over:
`nt`, `ntz`, `numtoolsd as ntd`, `numtoolsm as ntm`.
