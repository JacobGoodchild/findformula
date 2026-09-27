"""Search polynomial continued fractions  b(0) + a(1)/(b(1) + a(2)/(b(2) + ...))
and try to identify their limits as Mobius transforms of famous constants:
    x = (p + q*c) / (r + s*c)   with small integers p,q,r,s.
Pure computation, no network, no file writes except the results file.
"""
import sys, itertools, json
from multiprocessing import Pool
import mpmath as mp

mp.mp.dps = 60
CONSTS = {
    "pi": mp.pi, "pi^2": mp.pi ** 2, "zeta3": mp.zeta(3), "log2": mp.log(2),
    "catalan": mp.catalan, "e": mp.e, "log3": mp.log(3), "sqrt3*pi": mp.sqrt(3) * mp.pi,
    "zeta5": mp.zeta(5), "pi^4": mp.pi ** 4, "log2^2": mp.log(2) ** 2, "pi*log2": mp.pi * mp.log(2),
    "sqrt2": mp.sqrt(2), "euler_gamma": mp.euler, "1/pi": 1 / mp.pi, "pi^3": mp.pi ** 3,
    "e^2": mp.e ** 2, "sqrt(e)": mp.sqrt(mp.e), "log(1+sqrt2)": mp.log(1 + mp.sqrt(2)),
    "sqrt2*pi": mp.sqrt(2) * mp.pi,
}


def cf_value(a, b, depth):
    p0, p1 = 1, b(0)
    q0, q1 = 0, 1
    for n in range(1, depth + 1):
        an, bn = a(n), b(n)
        p0, p1 = p1, bn * p1 + an * p0
        q0, q1 = q1, bn * q1 + an * q0
    return p1, q1


def evaluate(a, b):
    try:
        P1, Q1 = cf_value(a, b, 400)
        P2, Q2 = cf_value(a, b, 600)
    except Exception:
        return None
    if Q1 == 0 or Q2 == 0:
        return None
    x1 = mp.mpf(P1) / Q1
    x2 = mp.mpf(P2) / Q2
    if not mp.isfinite(x2) or abs(x2) > 1e6:
        return None
    if abs(x1 - x2) > mp.mpf(10) ** -40 * max(1, abs(x2)):
        return None  # not converged fast enough to identify
    return x2


def identify(x):
    mp.mp.dps = 40
    rel = mp.pslq([x, 1], maxcoeff=10 ** 4, maxsteps=2000)
    if rel is not None:
        return None  # rational, boring
    for name, c in CONSTS.items():
        rel = mp.pslq([x * c, x, c, 1], maxcoeff=2000, maxsteps=4000)
        if rel is not None:
            s, r, q, p = rel  # s*x*c + r*x + q*c + p = 0 -> x = -(p+q c)/(r+s c)
            if s == 0 and r == 0:
                continue
            # verify with 58 digits
            mp.mp.dps = 60
            chk = abs(s * x * CONSTS[name] + r * x + q * CONSTS[name] + p)
            mp.mp.dps = 40
            if chk < mp.mpf(10) ** -50:
                return name, [int(t) for t in rel]
    return None


def work(params):
    fam, bc, ac = params
    b = lambda n: sum(c * n ** i for i, c in enumerate(bc))
    if fam == "lin":
        a = lambda n: sum(c * n ** i for i, c in enumerate(ac))
    else:
        s, us = ac[0], ac[1:]
        def a(n, s=s, us=us):
            v = s
            for u in us:
                v *= (n + u)
            return v
    mp.mp.dps = 60
    x = evaluate(a, b)
    if x is None:
        return None
    ident = identify(x)
    if ident is None:
        return None
    return {"fam": fam, "b": bc, "a": ac, "x": mp.nstr(x, 30), "const": ident[0], "rel": ident[1]}


def tasks(which):
    if which == "A":  # b linear, a quadratic
        for b1 in range(1, 6):
            for b0 in range(-3, 6):
                for a0, a1, a2 in itertools.product(range(-5, 6), repeat=3):
                    if a2 == 0 and a1 == 0:
                        continue
                    yield ("lin", (b0, b1), (a0, a1, a2))
    if which == "B":  # b quadratic, a = s * prod (n+u_i), 4 factors
        for b2 in range(1, 7):
            for b1 in range(-6, 7):
                for b0 in range(-3, 7):
                    for us in itertools.combinations_with_replacement(range(-1, 3), 4):
                        for s in (-4, -2, -1, 1, 2, 4):
                            yield ("prod", (b0, b1, b2), (s,) + us)
    if which == "C":  # b cubic, a = s * prod (n+u_i), 6 factors (Apery-like)
        for b3 in range(1, 41):
            for b2 in range(0, 3 * b3 // 2 + 1):
                for b1 in range(0, b2 + 1):
                    for b0 in range(0, b1 + 1):
                        for us in [(0,) * 6, (0, 0, 0, 1, 1, 1), (0,) * 3 + (-1,) * 3, (0, 0, 1, 1, 1, 1)]:
                            for s in (-1, -2, -4, -8, -16):
                                yield ("prod", (b0, b1, b2, b3), (s,) + us)


if __name__ == "__main__":
    which = sys.argv[1]
    out = open(f"cf_hits_{which}.jsonl", "w")
    n = 0
    with Pool(3) as pool:
        for res in pool.imap_unordered(work, tasks(which), chunksize=200):
            n += 1
            if res:
                out.write(json.dumps(res) + "\n")
                out.flush()
            if n % 20000 == 0:
                print(which, n, flush=True)
    print("done", which, n, flush=True)
