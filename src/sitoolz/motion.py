"""motion -- the one tick law of speed, shared by the site, the engine and here.

    a      = v + throttle*accel*dt
    b      = a - a*drag*dt
    v_next = clamp(b, -reverse_cap, forward_cap)

Accelerate, then drag the new speed, then clamp, as three statements, at 120 Hz. For constant
inputs below the cap the speed approaches ``throttle*accel*(1/drag - dt)`` as q^n with
``q = 1 - drag*dt``; the cap is a ceiling, the approach is acceleration against drag.

    from sitoolz import motion
    motion.settles_at(1, 12, 0.7)        # 17.042857142857...
    motion.cap_tick(20, 16)              # 363: the tick V=20, A=16 reaches its cap
"""
import math

HZ = 120
DT = 1.0 / HZ


def longitudinal_step(speed, throttle, accel, drag, dt=DT, reverse_cap=None, forward_cap=None):
    nxt = speed
    nxt += throttle * accel * dt
    nxt -= nxt * drag * dt
    if reverse_cap is not None and nxt < -reverse_cap:
        nxt = -reverse_cap
    if forward_cap is not None and nxt > forward_cap:
        nxt = forward_cap
    return nxt


def drag_of(brake, wet=0.0):
    """The studio's drag: 3 - 1.3*wet braking, 0.7 coasting."""
    return (3.0 - 1.3 * wet) if brake else 0.7


def settles_at(throttle, accel, drag, dt=DT):
    """The speed a constant throttle approaches below the cap."""
    return throttle * accel * (1.0 / drag - dt)


def run(v0, throttle, accel, drag, ticks, dt=DT, reverse_cap=None, forward_cap=None):
    """The speed after ``ticks`` ticks, by the recurrence."""
    v = v0
    for _ in range(ticks):
        v = longitudinal_step(v, throttle, accel, drag, dt, reverse_cap, forward_cap)
    return v


def closed_form(v0, throttle, accel, drag, n, dt=DT):
    """The same speed in one step, below the cap: v_n = v_inf + (v0 - v_inf) q^n."""
    q = 1.0 - drag * dt
    v_inf = settles_at(throttle, accel, drag, dt)
    return v_inf + (v0 - v_inf) * q ** n


def distance(v0, throttle, accel, drag, n, dt=DT):
    """Signed distance after n ticks along a fixed heading, position taking the new speed each
    step: x_n = h*[n*v_inf + (v0 - v_inf) q (1 - q^n) / (1 - q)]."""
    q = 1.0 - drag * dt
    v_inf = settles_at(throttle, accel, drag, dt)
    return dt * (n * v_inf + (v0 - v_inf) * q * (1 - q ** n) / (1 - q))


def cap_tick(speed_cap, accel, drag=0.7, dt=DT, limit=100000):
    """The first tick at which full throttle from rest reaches the cap, or None."""
    v = 0.0
    for n in range(1, limit + 1):
        v = longitudinal_step(v, 1.0, accel, drag, dt, 0.45 * speed_cap, speed_cap)
        if v >= speed_cap:
            return n
    return None


def decode_input(code):
    """The studio's input code: (throttle+1) + 3*(steer+1) + 9*brake -> (throttle, steer, brake)."""
    return code % 3 - 1, (code // 3) % 3 - 1, code // 9 == 1


def initial(x=0.0, z=0.0, yaw=0.0):
    """The ten fields of a flat-ground pose, at rest."""
    return {"x": x, "z": z, "yaw": yaw, "speed": 0.0, "vx": 0.0, "vz": 0.0,
            "steer_vis": 0.0, "lean": 0.0, "squat": 0.0, "braking": False}


def wrap(a):
    """An angle into (-pi, pi], as the engine's ``wrap``."""
    a = math.fmod(a + math.pi, 2 * math.pi)
    if a < 0:
        a += 2 * math.pi
    return a - math.pi


def step(pose, code, handling, dt=DT):
    """One tick of the whole flat-ground law: accelerate, drag, clamp, turn, the cosmetics, blend
    the velocity, move. ``handling`` has ``speed, accel, turn, grip, wet``; ``code`` is the studio's
    input code. Returns the next pose; the engine's ``step_drive`` on flat ground, bit for bit.

    >>> h = {"speed": 30, "accel": 12, "turn": 1.2, "grip": 1, "wet": 0}
    >>> p = initial(0, 0, 0)
    >>> for _ in range(120): p = step(p, 5, h)
    >>> round(p["speed"], 6), round(p["x"], 6)
    (8.596954, 4.833136)
    """
    throttle, steer, brake = decode_input(code)
    n = dict(pose)
    grip = handling["grip"] if handling["grip"] >= 1 else handling["grip"] - 0.055 * handling["wet"]
    top = handling["speed"] * (0.45 if throttle < 0 else 1)
    n["speed"] = longitudinal_step(n["speed"], throttle, handling["accel"], drag_of(brake, handling["wet"]), dt,
                                   handling["speed"] * 0.45, top)
    turn = handling["turn"] * steer * min(abs(n["speed"]) / 5, 1) * (-1 if n["speed"] < 0 else 1) * (1.5 if brake else 1)
    n["yaw"] = wrap(n["yaw"] + turn * dt)
    k, sp = min(dt * 8, 1), min(abs(n["speed"]) / 22, 1)
    n["steer_vis"] += (steer * 0.5 - n["steer_vis"]) * k
    n["lean"] += (steer * sp * 0.11 * (1 if handling["turn"] > 2.2 else 0.6) - n["lean"]) * min(dt * 5, 1)
    n["squat"] += ((-0.05 if brake else throttle * 0.035 * (1 - sp)) - n["squat"]) * min(dt * 6, 1)
    n["braking"] = bool(brake or throttle < 0)
    fx, fz = math.cos(n["yaw"]), math.sin(n["yaw"])
    g = math.pow(grip, dt * 60)
    n["vx"] = n["vx"] * g + fx * n["speed"] * (1 - g)
    n["vz"] = n["vz"] * g + fz * n["speed"] * (1 - g)
    if grip >= 1:
        n["vx"], n["vz"] = fx * n["speed"], fz * n["speed"]
    n["x"] += n["vx"] * dt
    n["z"] += n["vz"] * dt
    return n
