"""twin_k2n256 -- the 82-digit twin pair k*2^256 +/- 1, k = 19902, proven on both sides.

Pure Python, standard library only, a few seconds. The +1 side is proved by generalized
Pocklington (N-1 = k*2^256 is factored by construction, Proth 1878 / Brillhart-Lehmer-Selfridge
1975), the -1 side by the Lucas N+1 test (N+1 = k*2^256, the same factorization). The Lucas
prover self-tests against every prime and composite in [3, 4000) before it is trusted.

    from sitoolz import twin_k2n256
    twin_k2n256.verify()        # {'upper_proven': True, 'lower_proven': True, ...}

or ``python -m sitoolz.twin_k2n256`` for the printed report.
"""
from math import gcd, isqrt

k, n = 19902, 256
B = 1 << n
lo, hi = k*B - 1, k*B + 1
FACTORS = {2: 5, 3: 2, 31: 2, 107: 2}          # primes of k*2^n, with Pocklington witnesses
FPRIMES = [2, 3, 31, 107]                       # distinct primes of N-1 = N+1 = k*2^n

# ---------- Miller-Rabin (a quick screen, not the proof) ----------
def miller_rabin(num, bases=(2,3,5,7,11,13,17,19,23,29,31,37)):
    if num < 2: return False
    d, r = num-1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in bases:
        if a % num == 0: continue
        x = pow(a, d, num)
        if x in (1, num-1): continue
        for _ in range(r-1):
            x = x*x % num
            if x == num-1: break
        else: return False
    return True

# ---------- N-1 : generalized Pocklington (proves the +1 side) ----------
def pocklington_ok(N, factors):
    if not all(pow(a, N-1, N) == 1 and gcd(pow(a, (N-1)//p, N) - 1, N) == 1
               for p, a in factors.items()):
        return False
    return B*B > N                              # factored part 2^n alone exceeds sqrt(N)

# ---------- N+1 : Lucas / BLS (proves the -1 side) ----------
def jacobi(a, nn):
    a %= nn; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if nn % 8 in (3, 5): r = -r
        a, nn = nn, a
        if a % 4 == 3 and nn % 4 == 3: r = -r
        a %= nn
    return r if nn == 1 else 0

def lucas_uv(m, P, Q, mod):
    D = P*P - 4*Q; inv2 = pow(2, -1, mod)
    U, V, Qk = 0, 2, 1
    for b in bin(m)[2:]:
        U, V = (U*V) % mod, (V*V - 2*Qk) % mod; Qk = (Qk*Qk) % mod
        if b == '1':
            U, V = ((P*U + V)*inv2) % mod, ((D*U + P*V)*inv2) % mod; Qk = (Qk*Q) % mod
    return U, V

def lucas_np1_prime(N, fprimes, tries=60):
    if N < 3: return N == 2
    if N % 2 == 0: return False
    if isqrt(N)**2 == N: return False
    Ds, D = [], 5
    while len(Ds) < tries:
        j = jacobi(D, N)
        if j == 0 and abs(D) % N != 0: return False
        if j == -1: Ds.append(D)
        D = -D + 2 if D < 0 else -D - 2
    for D in Ds:
        for P in range(1, 40, 2):
            Q = (P*P - D) // 4
            if (P*P - D) % 4 or gcd(N, 2*Q*D) != 1: continue
            if lucas_uv(N+1, P, Q, N)[0] != 0: continue
            if all(gcd(lucas_uv((N+1)//q, P, Q, N)[0], N) == 1 for q in fprimes):
                return True
    return False

# ---------- self-test the Lucas prover before trusting it ----------
def _sieve(n):
    s = [True]*n; s[0]=s[1]=False
    for i in range(2, isqrt(n)+1):
        if s[i]:
            for j in range(i*i, n, i): s[j] = False
    return s
def _factorset(x):
    f=set(); d=2
    while d*d<=x:
        while x%d==0: f.add(d); x//=d
        d+=1
    if x>1: f.add(x)
    return f
def selftest():
    prime = _sieve(4000)
    for N in range(3, 4000):
        if lucas_np1_prime(N, list(_factorset(N+1))) != prime[N]:
            return False, N
    return True, None



def verify():
    """Prove both sides and return the facts as a dict."""
    import time
    t = time.time()
    ok, bad = selftest()
    out = {
        "k": k, "n": n, "digits": len(str(hi)),
        "lower": lo, "upper": hi,
        "difference": hi - lo,
        "midpoint_divisible_by_3": (lo + 1) % 3 == 0,
        "lucas_selftest": ok,
        "upper_proven": pocklington_ok(hi, FACTORS),
        "lower_proven": lucas_np1_prime(lo, FPRIMES),
        "witnesses": dict(FACTORS),
        "seconds": round(time.time() - t, 2),
    }
    return out


if __name__ == "__main__":
    ok, bad = selftest()
    print(f"Lucas N+1 prover self-test on [3,4000): {'PASS' if ok else f'FAIL at {bad}'}")
    print()
    print("pair straddles k*2^256, difference == 2 :", hi - lo == 2)
    print("midpoint divisible by 3                 :", (lo+1) % 3 == 0)
    print("UPPER  k*2^256 + 1  Pocklington-PROVEN   :", pocklington_ok(hi, FACTORS))
    print("LOWER  k*2^256 - 1  Lucas/BLS-PROVEN     :", lucas_np1_prime(lo, FPRIMES))
    print("(Miller-Rabin screen, both):",
          miller_rabin(lo) and miller_rabin(hi))
