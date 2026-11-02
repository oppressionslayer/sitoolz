"""oeis_crosscheck -- verify A094603/A094604/A094605/A094606 against each other.

The four rule 30 edge sequences are not independent data.  Any one of them
determines the others, so a term computed at great expense in one can be
checked for free from another -- without re-running a single generation.

The identity, from the light-cone edge law (proof: PROOF_V2_LAW.md):

    run(t) = #{ b : p(b) <= v2(t+1) }

where p(b) is the base-2 log of the period of diagonal b.  Reading it as a
cumulative count of the period spectrum gives:

    A094604[v] = #{ b : p(b) <= v }                    (records = cumulative)
    A094606[b] = p(b)                                  (the spectrum)
    A094605[b] = 2^p(b)                                (the periods)
    A094603[n] = A094604[v2(n)]                        (every row, not just records)

so v occurs exactly A094604[v] - A094604[v-1] times in A094606.

    from sitoolz import oeis_crosscheck as x

    x.spectrum(x.A094604)          # -> A094606, 137 terms from 58
    x.periods(x.A094604)           # -> A094605
    x.rows(x.A094604, 40)          # -> A094603, first 40 rows
    x.records(x.A094606)           # the inverse direction
    x.check()                      # everything against everything

Stdlib only.  No simulation anywhere in this file.
"""

# ---- published data, fetched 2026-08-22 -------------------------------------
A094604 = [1, 3, 4, 6, 7, 9, 15, 16, 24, 25, 27, 29, 34, 36, 37, 39, 41, 43,
           48, 49, 51, 54, 55, 58, 60, 63, 64, 66, 69, 70, 72, 74, 77, 79, 80,
           82, 84, 86, 90, 91, 93, 100, 103, 104, 106, 108, 111, 112, 114, 119,
           120, 122, 124, 127, 128, 133, 135, 137]

A094606 = [0, 1, 1, 2, 3, 3, 4, 5, 5, 6, 6, 6, 6, 6, 6, 7, 8, 8, 8, 8, 8, 8, 8,
           8, 9, 10, 10, 11, 11, 12, 12, 12, 12, 12, 13, 13, 14, 15, 15, 16, 16,
           17, 17, 18, 18, 18, 18, 18, 19, 20, 20, 21, 21, 21, 22, 23, 23, 23,
           24, 24, 25, 25, 25, 26, 27, 27, 28, 28, 28, 29]

A094605 = [1, 2, 2, 4, 8, 8, 16, 32, 32, 64, 64, 64, 64, 64, 64, 128, 256, 256,
           256, 256, 256, 256, 256, 256, 512, 1024, 1024, 2048, 2048, 4096,
           4096, 4096, 4096, 4096, 8192, 8192, 16384, 32768, 32768, 65536,
           65536, 131072, 131072, 262144, 262144, 262144, 262144]

A094603 = [1, 3, 1, 4, 1, 3, 1, 6, 1, 3, 1, 4, 1, 3, 1, 7, 1, 3, 1, 4, 1, 3, 1,
           6, 1, 3, 1, 4, 1, 3, 1, 9, 1, 3, 1, 4, 1, 3, 1, 6, 1, 3, 1, 4, 1, 3,
           1, 7, 1, 3, 1, 4, 1, 3, 1, 6, 1, 3, 1, 4, 1, 3, 1, 15]


def v2(n):
    """2-adic valuation -- OEIS A007814."""
    return (n & -n).bit_length() - 1


# ---- the four directions -----------------------------------------------------
def spectrum(records=None):
    """A094606 from A094604.  v occurs records[v] - records[v-1] times."""
    records = records or A094604
    out, prev = [], 0
    for v, r in enumerate(records):
        out += [v] * (r - prev)
        prev = r
    return out


def periods(records=None):
    """A094605 from A094604 -- the spectrum, exponentiated."""
    return [1 << p for p in spectrum(records)]


def rows(records=None, upto=64):
    """A094603 from A094604: the run in row n, for every n, not just records."""
    records = records or A094604
    return [records[v2(n)] for n in range(1, upto + 1)]


def records(spec=None):
    """A094604 from A094606 -- the inverse direction, a cumulative count."""
    spec = spec or A094606
    out = []
    for v in range(max(spec) + 1):
        out.append(sum(1 for p in spec if p <= v))
    return out


# ---- invariants the proof forces ---------------------------------------------
def invariants(records=None):
    """Structural checks any correct A094604 must satisfy.

    From Lemma 0 of the edge law: p is nondecreasing and steps by 0 or 1, and
    is unbounded.  So every exponent v is attained by at least one diagonal,
    which forces A094604 to be *strictly* increasing.  A repeated value would
    mean some period 2^v is attained by no diagonal at all -- impossible.
    """
    records = records or A094604
    d = [records[i] - records[i - 1] for i in range(1, len(records))]
    return {
        "strictly_increasing": all(x >= 1 for x in d),
        "first_term_is_1": records[0] == 1,
        "min_plateau": min(d), "max_plateau": max(d),
        "violations": [(i, records[i - 1], records[i])
                       for i, x in enumerate(d, 1) if x < 1],
    }


# ---- the cross-check ----------------------------------------------------------
def cross_check(a4=None, a6=None, a5=None, a3=None, verbose=True):
    """Every published sequence against every other, via the identity."""
    a4 = a4 or A094604; a6 = a6 or A094606
    a5 = a5 or A094605; a3 = a3 or A094603
    d6, d5, d3 = spectrum(a4), periods(a4), rows(a4, len(a3))
    r4 = records(a6)
    res = {
        "A094606": {"published": len(a6), "derived": len(d6),
                    "agree_on_overlap": d6[:len(a6)] == a6,
                    "new_terms": max(0, len(d6) - len(a6))},
        "A094605": {"published": len(a5), "derived": len(d5),
                    "agree_on_overlap": d5[:len(a5)] == a5,
                    "new_terms": max(0, len(d5) - len(a5))},
        "A094603": {"published": len(a3), "derived": len(d3),
                    "agree_on_overlap": d3[:len(a3)] == a3,
                    "new_terms": 0},
        "A094604_from_A094606": {"agree_on_overlap": r4 == a4[:len(r4)]},
        "invariants": invariants(a4),
    }
    if verbose:
        print("cross-check of the four rule 30 edge sequences")
        print("(no simulation -- everything below is the identity)\n")
        for k in ("A094606", "A094605", "A094603"):
            r = res[k]
            print(f"  {k}: published {r['published']:>4} terms, derived "
                  f"{r['derived']:>4}, overlap agrees: {r['agree_on_overlap']}"
                  + (f"  -> {r['new_terms']} NEW" if r["new_terms"] else ""))
        print(f"  A094604 rebuilt from A094606: "
              f"{res['A094604_from_A094606']['agree_on_overlap']}")
        inv = res["invariants"]
        print(f"\n  invariants: strictly increasing {inv['strictly_increasing']}, "
              f"plateaus {inv['min_plateau']}..{inv['max_plateau']}, "
              f"violations {len(inv['violations'])}")
    return res


def new_terms(a4=None, a6=None, a5=None):
    """The terms the identity supplies that are not yet published."""
    a4 = a4 or A094604; a6 = a6 or A094606; a5 = a5 or A094605
    return {"A094606": spectrum(a4)[len(a6):], "A094605": periods(a4)[len(a5):]}


def check():
    r = cross_check()
    n = new_terms()
    print(f"\nA094606, terms {len(A094606)+1}-{len(spectrum())}:")
    print("  " + ", ".join(map(str, n["A094606"])))
    print(f"\nA094605, terms {len(A094605)+1}-{len(periods())}:")
    print("  " + ", ".join(map(str, n["A094605"][:24])) + ", ...")
    return r


if __name__ == "__main__":
    check()
