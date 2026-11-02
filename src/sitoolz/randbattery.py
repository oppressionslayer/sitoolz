#!/usr/bin/env python3
"""
randbattery.py -- a small, dependency-free randomness battery for uint word
                  streams.  Point it at any source of integers and it runs the
                  classic tests with real p-values, so you can compare a
                  candidate generator against a known-good PRNG and against
                  physical random data on equal footing.

    from sitoolz import randbattery as rb
    rb.report(words, nbits=32, name="my source")     # one source -> a table
    rb.compare({"mine": w1, "MT19937": w2, ...})      # many sources, side by side

Tests (NIST SP800-22 style where applicable), p-value < 0.01 => FAIL:
    monobit        bit frequency (are 0s and 1s balanced?)
    runs           run count vs expected (too few/many switches?)
    serial_corr    lag-1 autocorrelation of the word stream
    byte_chisq     byte-value uniformity over 256 bins
    bitpos         worst per-bit-position bias across the word width

Pure stdlib (uses math.erfc); scipy is used only if present for the chi-square
tail, else a normal approximation.
"""
import math

try:
    from scipy.stats import chi2 as _chi2
    def _chisq_p(stat, dof): return float(_chi2.sf(stat, dof))
except Exception:
    def _chisq_p(stat, dof):                 # Wilson-Hilferty normal approximation
        z = ((stat / dof) ** (1 / 3) - (1 - 2 / (9 * dof))) / math.sqrt(2 / (9 * dof))
        return 0.5 * math.erfc(z / math.sqrt(2))

def _p_from_z(z):
    return math.erfc(abs(z) / math.sqrt(2))


def _bits(words, nbits):
    for w in words:
        for i in range(nbits - 1, -1, -1):
            yield (w >> i) & 1


# ------------------------------------------------------------------ the tests
def monobit(words, nbits):
    n = len(words) * nbits
    ones = sum(bin(w & ((1 << nbits) - 1)).count("1") for w in words)
    s = abs(2 * ones - n) / math.sqrt(n)
    return dict(test="monobit", stat=round(ones / n, 5), p=round(math.erfc(s / math.sqrt(2)), 4))

def runs(words, nbits):
    b = list(_bits(words, nbits))
    n = len(b)
    pi = sum(b) / n
    if abs(pi - 0.5) >= 2 / math.sqrt(n):      # precondition fails -> not balanced enough
        return dict(test="runs", stat=round(pi, 5), p=0.0)
    vobs = 1 + sum(1 for i in range(1, n) if b[i] != b[i - 1])
    num = abs(vobs - 2 * n * pi * (1 - pi))
    den = 2 * math.sqrt(2 * n) * pi * (1 - pi)
    return dict(test="runs", stat=vobs, p=round(math.erfc(num / den), 4))

def serial_corr(words, nbits):
    n = len(words) - 1
    x = words[:-1]; y = words[1:]
    mx = sum(x) / n; my = sum(y) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    vx = sum((a - mx) ** 2 for a in x); vy = sum((b - my) ** 2 for b in y)
    r = cov / math.sqrt(vx * vy) if vx and vy else 0.0
    z = r * math.sqrt(n)                       # under H0, r ~ N(0, 1/n)
    return dict(test="serial_corr", stat=round(r, 5), p=round(_p_from_z(z), 4))

def byte_chisq(words, nbits):
    nbytes = nbits // 8
    counts = [0] * 256
    for w in words:
        for k in range(nbytes):
            counts[(w >> (8 * k)) & 0xFF] += 1
    total = sum(counts); exp = total / 256
    stat = sum((c - exp) ** 2 / exp for c in counts)
    return dict(test="byte_chisq", stat=round(stat, 1), p=round(_chisq_p(stat, 255), 4))

def bitpos(words, nbits):
    n = len(words)
    ones = [0] * nbits
    for w in words:
        for i in range(nbits):
            ones[i] += (w >> i) & 1
    # aggregate chi-square across all positions: each z_i ~ N(0,1), sum z_i^2 ~ chi2(nbits)
    zs = [(2 * o - n) / math.sqrt(n) for o in ones]
    stat = sum(z * z for z in zs)
    worst = max(range(nbits), key=lambda i: abs(zs[i]))
    return dict(test="bitpos", stat=f"chi2={stat:.1f} worst=bit{worst}:{ones[worst]/n:.4f}",
                p=round(_chisq_p(stat, nbits), 4))


TESTS = [monobit, runs, serial_corr, byte_chisq, bitpos]

def report(words, nbits=32, name="source"):
    words = list(words)
    rows = []
    for t in TESTS:
        r = t(words, nbits)
        r["pass"] = "PASS" if r["p"] >= 0.01 else "FAIL"
        r["source"] = name
        rows.append(r)
    try:
        import pandas as pd
        return pd.DataFrame(rows)[["source", "test", "stat", "p", "pass"]]
    except Exception:
        return rows

def compare(sources, nbits=32):
    """sources: {name: word_list}.  Returns one tidy table, p-value per test."""
    frames = [report(w, nbits, name) for name, w in sources.items()]
    try:
        import pandas as pd
        big = pd.concat(frames, ignore_index=True)
        piv = big.pivot(index="test", columns="source", values="p")
        return piv
    except Exception:
        return frames


if __name__ == "__main__":
    import random
    random.seed(12345)
    N = 65536
    mt = [random.getrandbits(32) for _ in range(N)]
    # a deliberately bad source: a linear counter (should FAIL serial_corr / bitpos)
    lcg = [(i * 2654435761) & 0xFFFFFFFF for i in range(N)]
    print(compare({"MT19937": mt, "counter*phi": lcg}))
