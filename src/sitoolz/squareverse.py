#!/usr/bin/env python3
"""
squareverse.py -- a whole reproducible universe from any seed.

Seed it with a word, and out comes a deterministic world: terrain, loot, dice,
names, whatever -- all powered by the square-index generator (square_rng), all
recreatable.  Same seed -> same universe, on any machine, forever.  "Saving" a
world is just remembering the word that made it.

    from sitoolz import squareverse as sv
    w = sv.World("dragon")                 # the universe of 'dragon'
    w.stream("combat").roll("3d6+2")       # dice, from an independent sub-stream
    w.stream("loot").pick({"gold": 70, "gem": 25, "artifact": 5})   # weighted
    w.stream("npc").name()                 # a procedural name
    w.fbm(x, y)                            # smooth terrain noise at a point

Why it's fun (and useful for sims/gaming):
  * seed with ANYTHING (word, number, dict) -> a shareable world seed;
  * independent SUB-STREAMS per subsystem, so terrain and loot never desync;
  * fully reproducible -> ship the seed, regenerate the world, no save files;
  * good-quality randomness (passes a full randomness battery -- see randbattery).

Non-cryptographic by design -- deterministic and recreatable is the whole point.
"""
import math
from sitoolz import square_rng


# ---------------------------------------------------------------- a richer stream
class Stream(square_rng.SquareRandom):
    """A SquareRandom with game/sim conveniences bolted on."""

    def roll(self, spec="1d6"):
        """Dice notation: '3d6', '1d20+5', '2d8-1'."""
        spec = spec.replace(" ", "").lower()
        n, _, rest = spec.partition("d")
        n = int(n or 1)
        mod = 0
        for sign in "+-":
            if sign in rest:
                rest, _, m = rest.partition(sign)
                mod = int(m) * (1 if sign == "+" else -1)
                break
        faces = int(rest)
        return sum(self.randint(1, faces) for _ in range(n)) + mod

    def pick(self, table):
        """Weighted choice.  table: {value: weight} or a list of (value, weight)."""
        items = list(table.items()) if isinstance(table, dict) else list(table)
        total = sum(w for _, w in items)
        r = self.random() * total
        upto = 0
        for value, weight in items:
            upto += weight
            if r < upto:
                return value
        return items[-1][0]

    def chance(self, p):
        """True with probability p."""
        return self.random() < p

    def name(self, syllables=3):
        """A pronounceable procedural name."""
        cons, vow = "bcdfghjklmnprstvw", "aeiou"
        s = "".join(self.choice(cons) + self.choice(vow) for _ in range(syllables))
        return s.capitalize()


# ---------------------------------------------------------------- the world
class World:
    """A deterministic universe.  Everything derives from one seed."""

    def __init__(self, seed):
        self.seed = seed
        # a permutation table (shuffled by the seed) drives fast coherent noise
        r = square_rng.SquareRandom(f"{seed}::squareverse::perm")
        p = list(range(256))
        r.shuffle(p)
        self._perm = p + p

    def stream(self, name):
        """An INDEPENDENT reproducible sub-stream, addressed by name.
        world.stream('terrain') and world.stream('loot') never interfere."""
        return Stream(f"{self.seed}::{name}")

    # -- coherent value noise (terrain, clouds, textures) --------------------
    def _lattice(self, ix, iy):
        h = self._perm[(self._perm[ix & 255] + (iy & 255)) & 255]
        return h / 255.0

    def noise2d(self, x, y):
        """Smooth value noise in [0,1) at a continuous point."""
        x0, y0 = math.floor(x), math.floor(y)
        fx, fy = x - x0, y - y0
        sx = fx * fx * (3 - 2 * fx)
        sy = fy * fy * (3 - 2 * fy)
        v00 = self._lattice(x0, y0);     v10 = self._lattice(x0 + 1, y0)
        v01 = self._lattice(x0, y0 + 1); v11 = self._lattice(x0 + 1, y0 + 1)
        top = v00 * (1 - sx) + v10 * sx
        bot = v01 * (1 - sx) + v11 * sx
        return top * (1 - sy) + bot * sy

    def fbm(self, x, y, octaves=5, lac=2.0, gain=0.5):
        """Fractal (multi-octave) noise -- natural-looking terrain."""
        amp, freq, total, norm = 1.0, 1.0, 0.0, 0.0
        for _ in range(octaves):
            total += amp * self.noise2d(x * freq, y * freq)
            norm += amp
            amp *= gain
            freq *= lac
        return total / norm

    def heightmap(self, h, w, scale=8.0, octaves=5):
        """A (h, w) numpy array of fbm terrain in [0,1]."""
        import numpy as np
        out = np.empty((h, w), dtype=float)
        for j in range(h):
            for i in range(w):
                out[j, i] = self.fbm(i / w * scale, j / h * scale, octaves)
        lo, hi = out.min(), out.max()
        return (out - lo) / (hi - lo + 1e-9)

    def __repr__(self):
        return f"World({self.seed!r})"


# ---------------------------------------------------------------- self-test / demo
if __name__ == "__main__":
    # IDIOM: hold a stream, then draw from it many times.  Creating a stream
    # costs ~0.08s (it primes); drawing is cheap.  Never make a stream in a loop.
    w = World("dragon")
    print("world:", w)
    combat = w.stream("combat")
    print("combat rolls (3d6+2):", [combat.roll("3d6+2") for _ in range(6)])
    loot_rng = w.stream("loot")
    loot = {"gold": 70, "gem": 25, "artifact": 5}
    print("loot drops        :", [loot_rng.pick(loot) for _ in range(8)])
    npc = w.stream("npc")
    print("npc names         :", [npc.name() for _ in range(5)])

    # reproducibility flex: recreate the world from the seed, get the same draws
    sA = World("atlantis").stream("x"); a = [sA.roll("1d100") for _ in range(5)]
    sB = World("atlantis").stream("x"); b = [sB.roll("1d100") for _ in range(5)]
    print("same seed -> same world:", a == b, a)
    sC = World("mu").stream("x")
    print("different seed differs :", a != [sC.roll("1d100") for _ in range(5)])

    # independent sub-streams don't interfere (hold each stream once)
    w2 = World("iso")
    terr = w2.stream("terrain"); t1 = [terr.random() for _ in range(3)]
    loot2 = w2.stream("loot"); [loot2.random() for _ in range(5000)]       # hammer loot
    terr2 = World("iso").stream("terrain"); t2 = [terr2.random() for _ in range(3)]
    print("sub-streams independent:", t1 == t2)

    # a tiny ascii terrain map
    print("\nterrain from the word 'dragon':")
    hm = World("dragon").heightmap(14, 44, scale=5)
    ramp = " .:-=+*#%@"
    for row in hm:
        print("  " + "".join(ramp[min(len(ramp)-1, int(v*len(ramp)))] for v in row))
