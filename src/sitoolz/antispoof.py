#!/usr/bin/env python3
"""
antispoof.py -- a covert spoof-detection lane library, powered by square_rng.

The idea in one breath: run a SEPARATE, deterministic sequence that LOOKS RANDOM
(so an attacker can't forge or predict it) but is PREDETERMINED for whoever holds
the seed (so they can verify it).  It never rejects or blocks anything -- it only
INFORMS: "that item was spoofed."  Nothing breaks, so the attacker gets no
feedback and never learns they were caught.

    from sitoolz import antispoof as af
    toks  = af.chained_tokens("seed", messages)      # the lane, sent alongside
    bad   = af.covert_check("seed", messages, toks)  # -> indices that were spoofed

Stream patterns: chained_tokens / covert_check (order+content+history),
inband_tokens, session_digest (whole-transcript receipt), timeboxed_tokens
(anti-replay), majority_spoofed (survives one leaked seed).

Image watermark: embed_watermark / tamper_map -- an invisible looks-random lane
spread over the pixels that localizes tampering without altering the picture.

Honest scope: a DETECTION/INFORMANT rail, not confidentiality and not a hardened
MAC.  "Looks statistically random" stops a blind spoofer cold; a determined
*predicting* adversary needs the seed step crypto-backed.  See SQUARE_RNG.md.
"""
from sitoolz import square_rng


# ---------------------------------------------------------------- core primitive
def token(seed, *parts, bits=32):
    """A looks-random-but-predetermined token from (seed, parts...)."""
    key = "|".join([str(seed)] + [str(p) for p in parts])
    return square_rng.SquareRandom(key).getrandbits(bits)


# ---------------------------------------------------------------- stream patterns
def chained_tokens(seed, payloads, bits=32):
    """One token per item, chained so it binds content + position + history."""
    out, prev = [], "INIT"
    for i, p in enumerate(payloads):
        t = token(seed, prev, i, p, bits=bits)
        out.append(t)
        prev = t
    return out


def covert_check(seed, payloads, tokens, bits=32):
    """INFORM which indices were spoofed.  Never rejects, never raises -- returns
    the list of bad indices (empty = all clear).  The stream is untouched."""
    expected = chained_tokens(seed, payloads, bits=bits)
    bad = [i for i, (g, e) in enumerate(zip(tokens, expected)) if g != e]
    if len(payloads) != len(tokens):
        bad.append(min(len(payloads), len(tokens)))
    return bad


def inband_tokens(seed, payloads, bits=32):
    """Unchained per-item tags (position + content). Ride inside each frame."""
    return [token(seed, i, p, bits=bits) for i, p in enumerate(payloads)]


def session_digest(seed, payloads, bits=64):
    """One rolling value proving the WHOLE transcript -- any edit changes it."""
    prev = "INIT"
    for i, p in enumerate(payloads):
        prev = token(seed, prev, i, p, bits=bits)
    return prev


def timeboxed_tokens(seed, epoch, payloads, bits=32):
    """Tokens bound to a time window -- a replay in another epoch goes stale."""
    return [token(seed, epoch, i, p, bits=bits) for i, p in enumerate(payloads)]


def majority_spoofed(seeds, payloads, lanes, bits=32):
    """N independent lanes, majority vote -- survives ONE leaked seed.
    lanes: list of token-lists, one per seed.  Returns indices the majority flags."""
    n = len(payloads)
    flagged = []
    for i in range(n):
        votes = sum(1 for s, L in zip(seeds, lanes)
                    if i >= len(L) or L[i] != token(s, i, payloads[i], bits=bits))
        if votes > len(seeds) // 2:
            flagged.append(i)
    return flagged


# ---------------------------------------------------------------- image watermark
def watermark_bits(seed, h, w):
    """A looks-random-but-predetermined bit for every pixel, from the seed."""
    import numpy as np
    rng = square_rng.SquareRandom(f"{seed}|watermark")
    raw = rng.bytes((h * w + 7) // 8)
    return np.unpackbits(np.frombuffer(raw, dtype=np.uint8))[:h * w].reshape(h, w)


def embed_watermark(img, seed, channel=2):
    """Set each pixel's least-significant bit to the predetermined watermark bit.
    Invisible to the eye, and it changes nothing about how the image displays --
    the lane rides *inside* the picture without breaking it."""
    import numpy as np
    a = np.array(img, dtype=np.uint8).copy()
    h, w = a.shape[:2]
    wm = watermark_bits(seed, h, w).astype(np.uint8)
    if a.ndim == 3:
        a[..., channel] = (a[..., channel] & 0xFE) | wm
    else:
        a = (a & 0xFE) | wm
    return a


def tamper_map(img, seed, block=16, channel=2):
    """INFORM where the watermark stopped matching -> which blocks were tampered.
    Returns a (H/block, W/block) array of match rates (1.0 = pristine, ~0.5 =
    tampered, since edited pixels' LSBs no longer track the lane).  Never rejects."""
    import numpy as np
    a = np.array(img, dtype=np.uint8)
    h, w = a.shape[:2]
    wm = watermark_bits(seed, h, w)
    lsb = (a[..., channel] if a.ndim == 3 else a) & 1
    match = (lsb == wm)
    H, W = h // block, w // block
    grid = np.zeros((H, W))
    for by in range(H):
        for bx in range(W):
            grid[by, bx] = match[by*block:(by+1)*block, bx*block:(bx+1)*block].mean()
    return grid


def is_authentic(img, seed, block=16, channel=2, thresh=0.9):
    """Overall verdict + the list of tampered block coordinates (informational)."""
    import numpy as np
    grid = tamper_map(img, seed, block, channel)
    bad = [(int(y), int(x)) for y, x in zip(*np.where(grid < thresh))]
    return (len(bad) == 0), bad, grid
