#!/usr/bin/env python3
"""
square_rng.py -- a drop-in random generator powered by the square-index whitener.

Seed it with *anything* (a word, a dict, a file's bytes) and it turns that seed
into an endless PRNG-quality stream by whitening it through Lars's square
index (n = r**2 + rem, recursed to leaves).  Because it subclasses
random.Random, the WHOLE stdlib API comes along for free:

    from sitoolz import square_rng
    rng = square_rng.SquareRandom("lars")
    rng.random()            # float in [0, 1)
    rng.randint(1, 100)     # ints
    rng.choice("abcdef")    # picks
    rng.shuffle(deck)       # in-place shuffle
    rng.sample(pool, 5)     # sampling
    rng.gauss(0, 1)         # normal variates
    rng.bytes(16)           # raw random bytes

Same seed -> same sequence (reproducible).  Different seed -> different stream.
Self-contained except for `square_frame` (the whitener); pure stdlib otherwise.
"""
import os
import pickle
import random

from sitoolz import square_frame as sf


class SquareRandom(random.Random):
    """random.Random whose entropy is the square-index whitening of the seed."""
    BLOCK_BYTES = 65536         # seed material whitened per refill (~8192 words)
    WORD = 32
    SKIP = 32                   # drop each block's leading isqrt-chain leaves (low-entropy top)

    # -- random.Random hooks -------------------------------------------------
    def seed(self, a=None, version=2):
        """Seed from any object.  bytes/str/int used directly; else pickled."""
        if a is None:
            a = os.urandom(16)
        if isinstance(a, (bytes, bytearray)):
            b = bytes(a)
        elif isinstance(a, str):
            b = a.encode()
        elif isinstance(a, int):
            a = abs(a)
            b = a.to_bytes((a.bit_length() + 7) // 8 or 1, "big")
        else:
            b = pickle.dumps(a, protocol=pickle.HIGHEST_PROTOCOL)
        self._seedbytes = b or b"\x00"
        self._counter = 0
        self._material = None      # feedback material from the previous block
        self._words = []
        self._i = 0
        # prime: whiten one block only to seed the feedback, then discard it, so the
        # first visible word already comes from an XOR-fed (avalanched) block -- else
        # near-identical seeds would share their leading words.
        self._refill()
        self._i = len(self._words)

    _M64 = 0xFFFFFFFFFFFFFFFF

    def _mix(self, salt):
        """A 64-bit fold of the WHOLE seed + counter (FNV-1a variant).  Used to
        build the value's header so the leftmost leaves (which read the top bits)
        depend on every seed byte -> tiny seed changes avalanche from word 0."""
        h = (0x811C9DC5 ^ salt ^ (self._counter * 0x9E3779B97F4A7C15)) & self._M64
        for byte in self._seedbytes:
            h = ((h ^ byte) * 0x0100000001B3) & self._M64
        return (h ^ (h >> 29)) & self._M64

    def _refill(self):
        """Fresh words by whitening a block.  The FIRST block is the seed tiled to
        size; every later block whitens the *previous block's random output*
        (output feedback), so each input is genuinely different -> the blocks are
        independent, not near-copies of the same seed."""
        base = self._seedbytes
        tile = (base * (self.BLOCK_BYTES // len(base) + 1))[:self.BLOCK_BYTES]
        if self._material is None:                       # first block: the seed tile
            buf = tile
        else:                                            # later: tile XOR prior output
            m = self._material                           # (any length) cycled to full size
            mrep = (m * (self.BLOCK_BYTES // len(m) + 1))[:self.BLOCK_BYTES]
            buf = bytes(a ^ b for a, b in zip(tile, mrep))
        # domain-separate with a counter + seed head so distinct seeds never collide
        # header = 32 bytes folded from the ENTIRE seed + counter, so the top of
        # val (hence the leftmost leaves) avalanches on any seed change.
        head = b"".join(self._mix(s).to_bytes(8, "big") for s in range(4))
        val = int.from_bytes(head + buf, "big")
        self._counter += 1
        # the leftmost leaves are the isqrt-chain (low-entropy low bits); the high
        # entropy lives in the rem leaves, later in pre-order -- so drop the head.
        self._words = sf.leaf_words(val, self.WORD, auto=True)[self.SKIP:]
        self._material = b"".join(w.to_bytes(4, "big") for w in self._words)
        self._i = 0

    def _next_word(self):
        if self._i >= len(self._words):
            self._refill()
        w = self._words[self._i]
        self._i += 1
        return w

    def getrandbits(self, k):
        if k <= 0:
            raise ValueError("number of bits must be positive")
        out = 0
        for _ in range((k + self.WORD - 1) // self.WORD):
            out = (out << self.WORD) | self._next_word()
        return out & ((1 << k) - 1)

    def random(self):
        return self.getrandbits(53) / (1 << 53)      # 53-bit float in [0, 1)

    # -- state (reproducible snapshots) --------------------------------------
    def getstate(self):
        return (self._seedbytes, self._counter, self._material, self._i, tuple(self._words))

    def setstate(self, state):
        self._seedbytes, self._counter, self._material, self._i, words = state
        self._words = list(words)

    # -- convenience ---------------------------------------------------------
    def bytes(self, n):
        return self.getrandbits(8 * n).to_bytes(n, "big") if n else b""

    def __repr__(self):
        s = self._seedbytes[:12]
        return f"SquareRandom(seed={s!r}{'...' if len(self._seedbytes) > 12 else ''})"


def SquareRandom_new(seed=None):
    return SquareRandom(seed)


if __name__ == "__main__":
    rng = SquareRandom("lars rocha")
    print("repr        :", rng)
    print("random()    :", round(rng.random(), 6))
    print("randint     :", [rng.randint(1, 100) for _ in range(8)])
    print("choice      :", [rng.choice("ABCDEFG") for _ in range(8)])
    deck = list(range(10)); rng.shuffle(deck)
    print("shuffle     :", deck)
    print("sample      :", rng.sample(range(1000), 6))
    print("gauss       :", [round(rng.gauss(0, 1), 2) for _ in range(6)])
    print("bytes(8)    :", rng.bytes(8).hex())

    # reproducibility: same seed -> same first draws
    a = SquareRandom(42); b = SquareRandom(42)
    same = [a.randint(0, 10**9) for _ in range(5)] == [b.randint(0, 10**9) for _ in range(5)]
    print("reproducible:", same)
    c = SquareRandom(43)
    diff = [a.random() for _ in range(3)] != [c.random() for _ in range(3)]
    print("seed-sensitive:", diff)
