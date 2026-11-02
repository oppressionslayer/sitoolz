"""Two tables, one offset -- your idea, wired to the chain.

Walk every number, keep a prime list, read off two indexes:

    coun  = len(primes[:-1])      drop the last, take the length
    counn = coun + 1              just for fun, a one offset

Then: the cursor sits on the prime  q = p_(2*futurex - 1),  and the offset
decides whether A' takes q's own index or the one before it.  That is the
whole "-1 or -2" branch -- it was never a pattern, it was which table you read.
"""

from sympy import isprime, prime

LIMIT = 400


def build(limit=LIMIT):
    """One pass over all the numbers, exactly as described."""
    primes, coun, counn = [], [], []
    for m in range(limit + 1):
        if isprime(m):
            primes.append(m)
        c = len(primes[:-1])          # drop the last, take the length
        coun.append(c)
        counn.append(c + 1)           # the one offset
    return primes, coun, counn


if __name__ == "__main__":
    primes, coun, counn = build()

    print("the two tables step at every prime, counn always one ahead:")
    print("    m :", " ".join(f"{m:3}" for m in range(20, 32)))
    print("  coun:", " ".join(f"{coun[m]:3}" for m in range(20, 32)))
    print(" counn:", " ".join(f"{counn[m]:3}" for m in range(20, 32)))

    src = open("/home/wiliamrocha/sequence/regen.py").read().split("if __name__")[0]
    exec(src)

    print("\nthe cursor's reference prime is q = p_(2*futurex - 1).")
    print("A' reads one of the two tables at q:\n")
    print("   k  futurex    q   coun[q] counn[q]   A'   reads    lands on"
          "    odd(W)")
    bits = []
    for k in range(1, 11):
        q = prime(2 * futurex(k) - 1)
        which = "counn" if A(k) == counn[q] else "coun "
        bits.append(1 if which == "counn" else 0)
        N = odd_part(W("B", k + 1))
        tag = "prime" if isprime(N) else ("unit" if N == 1 else "COMPOSITE")
        print(f"  {k:2}  {futurex(k):5} {q:6}  {coun[q]:6} {counn[q]:7} "
              f"{A(k):5}   {which}   {prime(A(k)):6}    {N:5} {tag}")

    print(f"\n  read counn -> take q.  read coun -> take the prime before q.")
    print(f"  bit pattern: {bits}")
    lo = [odd_part(W("B", k + 1)) for k in range(1, 11) if not bits[k - 1]]
    print(f"  every coun read has a prime or unit odd part: {lo}  "
          f"-> {all(isprime(v) or v == 1 for v in lo)}")
    hi = [odd_part(W("B", k + 1)) for k in range(1, 11) if bits[k - 1]]
    print(f"  counn reads are mixed: {hi}  <- both composites live here")
