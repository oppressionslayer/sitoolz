"""Map the factor anatomy around a number -- the smooth-head x giant-prime split.

Near any N, ~ln 2 (69%) of the offsets N+k are ``(small smooth head) x (one giant
prime > sqrt)`` -- the largest-prime-factor law (Dickman). This walks a window and
shows, per offset: its factorization, the largest prime factor and its exponent
``log_N(largest)``, whether that factor is a "giant" (> sqrt), and the smooth head
left when you divide the giant out. It also flags the nearest giant-prime offset --
a handy proof anchor when N-1 itself is too smooth for a one-step Pocklington/Lucas
certificate.
"""
from __future__ import annotations

import math

from .certify import factorize, is_prime
from ._base import REACH


def neighborhood(N: int, width: int = 10) -> dict:
    """Factor map of ``[N-width, N+width]``.

    Returns ``{center, prime, offsets, anchor}`` where ``offsets`` is a list of
    per-``k`` records::

        {k, value, factors, largest_prime, exponent, giant, smooth_head}

    ``exponent`` is ``log_value(largest_prime)`` (``> 0.5`` means a giant factor).
    ``anchor`` is the nearest offset (smallest ``|k|``) whose largest factor is a
    giant -- the best nearby handle for an N+-1 style proof -- or ``None``.
    """
    N = int(N)
    if N + width > REACH:
        raise ValueError(f"N+width exceeds provable reach {REACH:,}")
    offsets = []
    for k in range(-width, width + 1):
        m = N + k
        if m < 2:
            continue
        f = factorize(m)
        big = max(f)
        offsets.append({
            "k": k,
            "value": m,
            "factors": f,
            "largest_prime": big,
            "exponent": math.log(big) / math.log(m),
            "giant": big * big > m,
            "smooth_head": m // big,
        })
    # a proof anchor is a NEARBY offset (k != 0) carrying a giant prime factor
    giants = [r for r in offsets if r["giant"] and r["k"] != 0]
    anchor = min(giants, key=lambda r: (abs(r["k"]), r["k"])) if giants else None
    return {
        "center": N,
        "prime": is_prime(N) if N <= REACH else None,
        "offsets": offsets,
        "anchor": anchor,
    }


def show(N: int, width: int = 10) -> None:
    """Print :func:`neighborhood` as an aligned table."""
    info = neighborhood(N, width)
    print(f"N = {info['center']:,}  prime={info['prime']}")
    for r in info["offsets"]:
        tag = "  <== GIANT" if r["giant"] else ("  smooth" if r["exponent"] < 0.4 else "")
        print(f"{r['k']:+3d}  {r['value']:>18,}  head {r['smooth_head']:>16,} x "
              f"{r['largest_prime']:>18,}  N^{r['exponent']:.3f}{tag}")
    a = info["anchor"]
    if a:
        print(f"proof anchor: k={a['k']:+d}  giant prime {a['largest_prime']:,}")
