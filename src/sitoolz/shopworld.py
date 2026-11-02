#!/usr/bin/env python3
"""
shopworld.py -- the "clerk at a store" game, as an INDEX WORLD.

A shop is a mess of evolving state: inventory, a queue of random customers, cash,
reputation, events.  Normally you'd SAVE all of that -- a big, fragile save file.

Here we save NOTHING but the seed + the log of your actions.  The entire world,
and its exact status at any moment, is *recomputed* by replaying that log through
a deterministic index.  That is the whole trick of a recreatable open world:

    STORE THE HISTORY, NOT THE STATE.
    world_state = replay(seed, actions)      # a pure function

Consequences you get for free:
  * tiny saves (a word + a few letters) that recreate a rich world exactly;
  * "what if?" -- change one past action, replay, get a different world (a fork);
  * adaptivity -- the world reacts to how you've played, yet stays recreatable,
    because the reaction is keyed on your history (which is part of the address).

The "index" here is a cheap stable hash -- the seed IS the randomness.  Swap in
edge_at / square_rng / anything deterministic; the principle doesn't change.
"""
import hashlib

# ------------------------------------------------------------ the index lookup
def idx(*parts):
    """A stable deterministic lookup: any key -> a big integer.  This is the
    'index' -- content is ADDRESSED by (seed, ...), never stored."""
    key = "|".join(str(p) for p in parts).encode()
    return int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "big")

def pick(options, *key):
    return options[idx(*key) % len(options)]

def rng(lo, hi, *key):
    return lo + idx(*key) % (hi - lo + 1)


# ------------------------------------------------------------ the world content
ITEMS = ["apple", "bread", "milk", "egg", "cheese"]
NAMES = ["Mira", "Jax", "Odell", "Suki", "Bram", "Wren", "Cato", "Lux", "Ves", "Pia"]
RESTOCK_COST, RESTOCK_QTY, START_CASH = 2, 6, 25

def customer(seed, n, reputation):
    """Customer n is ADDRESSED, not stored.  Note: their generosity is keyed on
    your REPUTATION so far -- so the world ADAPTS to how you've played, while
    staying perfectly recreatable (reputation is itself derived from the log)."""
    tier = 0 if reputation < 3 else (1 if reputation < 8 else 2)   # rep opens up richer trade
    return {
        "n": n,
        "name": pick(NAMES, seed, "name", n),
        "item": pick(ITEMS, seed, "item", n, tier),
        "qty":  rng(1, 3, seed, "qty", n),
        "offer": rng(2, 4, seed, "offer", n) + tier,             # better rep -> better offers
    }


# ------------------------------------------------------------ the pure world function
def replay(seed, actions):
    """Recompute the ENTIRE world + status from just (seed, actions).  No state
    was ever stored; this rebuilds it deterministically."""
    inv = {it: rng(2, 5, seed, "start", it) for it in ITEMS}     # starting stock, addressed
    cash, rep, events = START_CASH, 0, []
    for n, act in enumerate(actions):
        c = customer(seed, n, rep)
        if act == "sell":
            if inv[c["item"]] >= c["qty"]:
                inv[c["item"]] -= c["qty"]; cash += c["offer"] * c["qty"]; rep += 1
                events.append(f"day {n}: sold {c['qty']} {c['item']} to {c['name']} (+${c['offer']*c['qty']})")
            else:
                rep -= 1
                events.append(f"day {n}: {c['name']} wanted {c['item']} -- OUT OF STOCK (rep-)")
        elif act.startswith("restock:"):
            it = act.split(":", 1)[1]
            inv[it] += RESTOCK_QTY; cash -= RESTOCK_COST * RESTOCK_QTY
            events.append(f"day {n}: restocked {it} (+{RESTOCK_QTY}, -${RESTOCK_COST*RESTOCK_QTY})")
        elif act == "skip":
            rep -= 1
            events.append(f"day {n}: skipped {c['name']} (rep-)")
    return {"day": len(actions), "cash": cash, "reputation": rep,
            "inventory": inv, "events": events}


# ------------------------------------------------------------ tiny save / load
def save(seed, actions):
    """The ENTIRE save file: a word and a few letters."""
    return seed + "|" + ",".join(actions)

def load(blob):
    seed, _, acts = blob.partition("|")
    return seed, (acts.split(",") if acts else [])


# ------------------------------------------------------------ interactive play
def play(seed="dragon"):
    actions = []
    print(f"SHOPWORLD -- seed '{seed}'.  serve customers, restock, survive.")
    print("commands: [enter]=sell to customer  r <item>=restock  x=skip  q=quit")
    while True:
        st = replay(seed, actions)
        c = customer(seed, len(actions), st["reputation"])
        print(f"\n--- day {st['day']}  cash=${st['cash']}  rep={st['reputation']} ---")
        print("  stock: " + "  ".join(f"{k}:{v}" for k, v in st["inventory"].items()))
        print(f"  >> {c['name']} wants {c['qty']} {c['item']} @ ${c['offer']} each")
        cmd = input("  action> ").strip().lower()
        if cmd == "q":
            print("\nSAVE FILE (recreates everything):", save(seed, actions)); return
        elif cmd.startswith("r "):
            it = cmd[2:].strip()
            actions.append(f"restock:{it}" if it in ITEMS else "skip")
        elif cmd == "x":
            actions.append("skip")
        else:
            actions.append("sell")


# ------------------------------------------------------------ the demonstration
if __name__ == "__main__":
    seed = "dragon"
    acts = ["sell", "restock:apple", "sell", "sell", "skip",
            "restock:milk", "sell", "sell", "sell", "restock:bread", "sell"]

    st = replay(seed, acts)
    print(f"played {st['day']} days of shop '{seed}':")
    for e in st["events"]:
        print("   " + e)
    print(f"\n  FINAL: cash=${st['cash']}  rep={st['reputation']}  stock={st['inventory']}")

    blob = save(seed, acts)
    print(f"\n  the ENTIRE save file ({len(blob)} bytes): {blob!r}")

    # 1. recreatable: reload the tiny blob, replay, identical world
    s2, a2 = load(blob)
    print(f"  reload + replay identical? {replay(s2, a2) == st}")

    # 2. what-if fork: change ONE early action -> the whole timeline diverges,
    #    because customers ADAPT to reputation, so an early rep change cascades.
    forked = acts.copy(); forked[0] = "skip"          # day 0: don't serve the first customer
    fk = replay(seed, forked)
    first_diff = next((i for i, (a, b) in enumerate(zip(st["events"], fk["events"])) if a != b), 0)
    print(f"  fork (day 0 sell->skip): cash ${st['cash']}->${fk['cash']}, "
          f"rep {st['reputation']}->{fk['reputation']}")
    print(f"    timelines diverge from day {first_diff}:")
    print(f"      original @day{first_diff}: {st['events'][first_diff]}")
    print(f"      forked   @day{first_diff}: {fk['events'][first_diff]}")
    print(f"    -> one changed choice = a whole different world, still 88 bytes to store")

    # 3. adaptivity: the SAME customer slot changes with reputation (world adapts)
    lowrep = customer(seed, 6, 0)["offer"]
    hirep  = customer(seed, 6, 10)["offer"]
    print(f"  adaptive: customer #6 offers ${lowrep} at low rep vs ${hirep} at high rep "
          f"-- world reacts to your history, still recreatable")
