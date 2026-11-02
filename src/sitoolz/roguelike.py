#!/usr/bin/env python3
"""
roguelike.py -- a seeded roguelike dungeon, generated entirely from a word.

Same word -> same dungeon, same loot, same monsters, forever.  The game CORE is
pure Python (renderer-agnostic): it just knows tiles, entities, and turns.  The
RENDERERS are swappable -- terminal to play now, PNG to look at, and pygame / web
(pygbag) / Taichi / OpenGL later drop in over the exact same core, no logic change.

    from sitoolz import roguelike as rl
    d = rl.Dungeon("dragon")     # a dungeon from the word 'dragon'
    print(d.render_ascii())
    d.render_png("dungeon.png")
    d.play()                     # interactive: w/a/s/d to move, q to quit

Powered by squareverse -> square_rng (deterministic, non-crypto, recreatable).
"""
from sitoolz import squareverse as sv


class Dungeon:
    WALL, FLOOR = "#", "."

    def __init__(self, seed, w=52, h=26, max_rooms=10):
        self.seed = seed
        self.w, self.h = w, h
        r = sv.World(seed).stream("dungeon")     # ONE held stream drives generation
        self.grid = [[self.WALL] * w for _ in range(h)]
        self.rooms = []
        for _ in range(max_rooms):
            rw, rh = r.randint(4, 9), r.randint(3, 6)
            rx, ry = r.randint(1, w - rw - 1), r.randint(1, h - rh - 1)
            if any(self._overlap(rx, ry, rw, rh, o) for o in self.rooms):
                continue
            self._carve_room(rx, ry, rw, rh)
            if self.rooms:                        # corridor to the previous room
                self._corridor(self.rooms[-1], (rx + rw // 2, ry + rh // 2), r)
            self.rooms.append((rx + rw // 2, ry + rh // 2))
        # entities
        self.px, self.py = self.rooms[0]
        self.sx, self.sy = self.rooms[-1]         # stairs (the exit)
        self.gold = set()
        self.monsters = {}
        loot = sv.World(seed).stream("loot")
        for cx, cy in self.rooms[1:]:
            for _ in range(loot.randint(0, 3)):
                gx, gy = cx + loot.randint(-2, 2), cy + loot.randint(-1, 1)
                if self._is_floor(gx, gy) and (gx, gy) != (self.sx, self.sy):
                    self.gold.add((gx, gy))
            if loot.chance(0.6):
                mx, my = cx + loot.randint(-2, 2), cy + loot.randint(-1, 1)
                if self._is_floor(mx, my) and (mx, my) != (self.px, self.py):
                    self.monsters[(mx, my)] = loot.randint(3, 8)   # monster HP
        self.hp, self.score, self.won, self.dead = 20, 0, False, False
        self.log = []
        self.fight = sv.World(seed).stream("combat")

    # -- generation helpers --------------------------------------------------
    def _overlap(self, x, y, w, h, c):
        return abs(c[0] - (x + w // 2)) < w and abs(c[1] - (y + h // 2)) < h
    def _carve_room(self, x, y, w, h):
        for j in range(y, y + h):
            for i in range(x, x + w):
                self.grid[j][i] = self.FLOOR
    def _corridor(self, a, b, r):
        (x0, y0), (x1, y1) = a, b
        if r.chance(0.5):
            for x in range(min(x0, x1), max(x0, x1) + 1): self.grid[y0][x] = self.FLOOR
            for y in range(min(y0, y1), max(y0, y1) + 1): self.grid[y][x1] = self.FLOOR
        else:
            for y in range(min(y0, y1), max(y0, y1) + 1): self.grid[y][x0] = self.FLOOR
            for x in range(min(x0, x1), max(x0, x1) + 1): self.grid[y1][x] = self.FLOOR
    def _is_floor(self, x, y):
        return 0 <= x < self.w and 0 <= y < self.h and self.grid[y][x] == self.FLOOR

    # -- game logic (renderer-agnostic) --------------------------------------
    def step(self, dx, dy):
        if self.won or self.dead:
            return
        nx, ny = self.px + dx, self.py + dy
        if not self._is_floor(nx, ny):
            return
        if (nx, ny) in self.monsters:                      # bump = fight
            you, foe = self.fight.roll("1d6+3"), self.fight.roll("1d6")
            if you >= foe:
                self.monsters[(nx, ny)] -= you
                if self.monsters[(nx, ny)] <= 0:
                    del self.monsters[(nx, ny)]
                    self.log.append("you slay a monster!")
                else:
                    self.log.append(f"you hit it ({you} dmg)")
            else:
                self.hp -= foe
                self.log.append(f"it hits you for {foe} (hp {self.hp})")
                if self.hp <= 0:
                    self.dead = True; self.log.append("you died.")
            return
        self.px, self.py = nx, ny
        if (nx, ny) in self.gold:
            self.gold.discard((nx, ny)); self.score += 10
            self.log.append(f"+10 gold (score {self.score})")
        if (nx, ny) == (self.sx, self.sy):
            self.won = True; self.log.append("you reached the stairs! you win!")

    # -- renderers (swap these; the core above never changes) ----------------
    def _glyph(self, x, y):
        if (x, y) == (self.px, self.py): return "@"
        if (x, y) == (self.sx, self.sy): return ">"
        if (x, y) in self.monsters: return "M"
        if (x, y) in self.gold: return "$"
        return self.grid[y][x]

    def render_ascii(self):
        rows = ["".join(self._glyph(x, y) for x in range(self.w)) for y in range(self.h)]
        status = f"  seed={self.seed!r}  hp={self.hp}  score={self.score}  " \
                 f"{'WON' if self.won else 'DEAD' if self.dead else 'exploring'}"
        return "\n".join(rows) + "\n" + status

    def render_png(self, path, tile=18):
        from PIL import Image, ImageDraw
        cols = {"#": (38, 42, 54), ".": (176, 170, 150), "@": (70, 200, 235),
                ">": (90, 220, 120), "M": (220, 70, 70), "$": (240, 205, 70)}
        img = Image.new("RGB", (self.w * tile, self.h * tile + tile), cols["#"])
        dr = ImageDraw.Draw(img)
        for y in range(self.h):
            for x in range(self.w):
                g = self._glyph(x, y)
                base = cols["."] if g in "@>M$" else cols[g]
                dr.rectangle([x*tile, y*tile, x*tile+tile-1, y*tile+tile-1], fill=base)
                if g in "@>M$":
                    dr.ellipse([x*tile+3, y*tile+3, x*tile+tile-4, y*tile+tile-4], fill=cols[g])
        dr.text((6, self.h*tile+3),
                f"seed='{self.seed}'  hp={self.hp}  score={self.score}", fill=(230, 230, 230))
        img.save(path)
        return path

    # -- interactive terminal renderer (for YOU to run) ----------------------
    def play(self):
        moves = {"w": (0, -1), "s": (0, 1), "a": (-1, 0), "d": (1, 0)}
        print("seeded roguelike -- w/a/s/d move, q quit.  @ you  $ gold  M monster  > exit")
        while not (self.won or self.dead):
            print("\n" + self.render_ascii())
            if self.log: print("  " + " | ".join(self.log[-3:]))
            c = input("> ").strip().lower()[:1]
            if c == "q": print("bye!"); return
            if c in moves: self.step(*moves[c])
        print("\n" + self.render_ascii())
        print("  " + (self.log[-1] if self.log else ""))


# ---------------------------------------------------------------- autopilot + demo
def _bfs(d, start, goal):
    from collections import deque
    q = deque([start]); prev = {start: None}
    while q:
        cur = q.popleft()
        if cur == goal: break
        for dx, dy in ((0,1),(0,-1),(1,0),(-1,0)):
            nx, ny = cur[0]+dx, cur[1]+dy
            if d._is_floor(nx, ny) and (nx, ny) not in prev:
                prev[(nx, ny)] = cur; q.append((nx, ny))
    if goal not in prev: return []
    path, cur = [], goal
    while cur: path.append(cur); cur = prev[cur]
    return path[::-1]

def autoplay(seed, verbose=False):
    """A deterministic bot walks to the exit -- proves the core plays & replays."""
    d = Dungeon(seed)
    path = _bfs(d, (d.px, d.py), (d.sx, d.sy))
    for (nx, ny) in path[1:]:
        d.step(nx - d.px, ny - d.py)
        if d.dead: break
    return d


if __name__ == "__main__":
    import sys
    seed = sys.argv[1] if len(sys.argv) > 1 else "dragon"
    d = Dungeon(seed)
    print(d.render_ascii())
    print(f"\n  rooms={len(d.rooms)}  gold={len(d.gold)}  monsters={len(d.monsters)}")
    # reproducibility: regenerate and compare
    same = Dungeon(seed).render_ascii() == d.render_ascii()
    print(f"  same seed -> same dungeon: {same}")
    # autoplay + replay determinism
    a, b = autoplay(seed), autoplay(seed)
    print(f"  autopilot result: score={a.score} hp={a.hp} won={a.won}; replay identical: "
          f"{(a.score,a.hp,a.won)==(b.score,b.hp,b.won)}")
