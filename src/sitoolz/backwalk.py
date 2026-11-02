"""backwalk -- the family test run backwards, without knowing the index.

`squarewalk` asks: given i, did the xor keep the root inside index i?  That is
cheap because you already hold i.

Backwards is the hard direction: given a number x and nothing else, is x in the
family -- is x = (i*i ^ m)**2 for some i?  Naively that is two square roots
(root of x, then root of root^m), and a square root is the expensive operation.

The opening is that **m is small, so everything the xor does lives in the low
bits**, and the low bits of a square root are recoverable 2-adically -- by
Hensel lifting -- with no square root anywhere.  So a filter that costs a fixed
number of machine words, no matter how many digits x has, can reject almost
everything before a single expensive root is taken.

    from sitoolz import backwalk as bw

    bw.sqrt_mod(x)             # low bits of a square root, 2-adically
    bw.filter_member(x, 5)     # cheap reject, constant cost in |x|
    bw.is_member(x, 5)         # exact: filter first, confirm only if it passes
    bw.rejection(5)            # what fraction the filter kills
    bw.bench_backwards(5)      # filter+confirm vs the naive two roots
"""
from math import isqrt

K = 64
MASK = (1 << K) - 1


def sqrt_mod(n, k=K):
    """A square root of n modulo 2^k, or None if n is not one there.

    Bit-by-bit Hensel lifting.  Every operation is on k-bit numbers, so the
    cost does not depend on how large the original number was -- only the low
    k bits of n are ever read.
    """
    m = (1 << k) - 1
    n &= m
    if n == 0:
        return 0
    # strip factors of 4; an odd square is 1 mod 8
    s = 0
    while n & 3 == 0:
        n >>= 2
        s += 1
        k -= 2
        m = (1 << k) - 1
    if n & 1 == 0:
        return None                      # 2 * odd is never a square
    if n & 7 != 1:
        return None
    r = 1
    for j in range(3, k):
        if ((r * r - n) >> j) & 1:
            r |= 1 << (j - 1)
        r &= m
    if (r * r - n) & m:
        return None
    return (r << s) & ((1 << (k + 2 * s)) - 1)


def low_square(n, k=K):
    """Is n a square modulo 2^k?  The cheap necessary condition."""
    return sqrt_mod(n, k) is not None


SQ63 = frozenset((i * i) % 63 for i in range(63))
SQ65 = frozenset((i * i) % 65 for i in range(65))
SQ11 = frozenset((i * i) % 11 for i in range(11))


def small_gate(x):
    """Classic residue filter: a square lands in 16/63, 21/65 and 6/11.

    Costs one single-limb division each, so it is O(limbs) rather than the
    2-adic gate's O(1) -- but the constant is tiny next to a square root, and
    it takes the pass rate from 1/24 to about 1/600.
    """
    return (x % 63 in SQ63) and (x % 65 in SQ65) and (x % 11 in SQ11)


def filter_member(x, mask=5, k=K, small=True):
    """Reject x as a family member using only its low bits.

    Two gates, both 2-adic and both constant cost in the size of x:

      1. x must be a square mod 2^k.
      2. its root r must have r ^ mask a square mod 2^k -- and because mask
         is small, r ^ mask is decided by the same low bits we already have.

    With `small=True` the residue gate above runs first, which is the bigger
    filter in practice.  A pass is not a proof; `is_member` confirms.
    """
    if small and not small_gate(x):
        return False
    r = sqrt_mod(x, k)
    if r is None:
        return False
    m2 = (1 << k) - 1
    for cand in (r, (-r) & m2):          # both roots mod 2^k
        if low_square((cand ^ mask) & m2, k):
            return True
    return False


def is_member(x, mask=5, k=K, small=True):
    """Exact membership: cheap filter first, expensive confirm only on a pass."""
    if not filter_member(x, mask, k, small):
        return False
    r = isqrt(x)
    if r * r != x:
        return False
    inner = r ^ mask
    s = isqrt(inner)
    return s * s == inner


def naive_member(x, mask=5):
    """The direct version: two square roots, always paid."""
    r = isqrt(x)
    if r * r != x:
        return False
    inner = r ^ mask
    s = isqrt(inner)
    return s * s == inner


def members(mask=5, lo=0, hi=64):
    """Actual family members, for testing."""
    return [(i * i ^ mask) ** 2 for i in range(lo, hi)]


def rejection(mask=5, n=20000, bits=256, k=K, seed=12345, small=True):
    """What fraction of random non-members the filter kills, and false passes."""
    rnd = seed
    killed = passed = 0
    for _ in range(n):
        rnd = (rnd * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        x = (rnd << (bits - 64)) | rnd
        if filter_member(x, mask, k, small):
            passed += 1
        else:
            killed += 1
    return {"tested": n, "rejected": killed, "passed": passed,
            "rejection_rate": killed / n,
            "survivors_needing_a_root": passed / n}


def sanity(mask=5, hi=200, k=K):
    """Every real member must pass the filter, and is_member must agree."""
    bad_filter = bad_exact = 0
    for x in members(mask, 0, hi):
        if not filter_member(x, mask, k):
            bad_filter += 1
        if is_member(x, mask, k) != naive_member(x, mask):
            bad_exact += 1
    return {"members_checked": hi, "filter_false_negatives": bad_filter,
            "is_member_disagreements": bad_exact}


def bench_backwards(mask=5, bits=4096, n=3000, seed=99, k=K):
    """Filter-then-confirm against the naive two roots, on random inputs."""
    import time
    rnd = seed
    xs = []
    for _ in range(n):
        rnd = (rnd * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        xs.append(((rnd << (bits - 64)) | rnd) | 1)
    t0 = time.perf_counter()
    a = sum(1 for x in xs if naive_member(x, mask))
    t1 = time.perf_counter()
    b = sum(1 for x in xs if is_member(x, mask, k))
    t2 = time.perf_counter()
    return {"inputs": n, "bits": bits, "same_answer": a == b,
            "naive_sec": round(t1 - t0, 4), "filtered_sec": round(t2 - t1, 4),
            "speedup": round((t1 - t0) / (t2 - t1), 1)}


if __name__ == "__main__":
    print("sanity (real members must all pass):", sanity())
    print("2-adic gate only:", rejection(small=False))
    print("both gates:      ", rejection(small=True))
    for bits in (512, 4096, 32768, 262144):
        print(f"bench {bits:>7} bits:", bench_backwards(bits=bits, n=400))
