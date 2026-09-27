"""Flip-flop Grover trapping on the LINE GRAPH of the FCC lattice via FCC Green functions.
Bond vectors (cubic coords): (1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1); sum_t cos(k.e_t) = 2s,
s = c1c2 + c2c3 + c3c1.  P = 6, degree 22, flat band mu = -1/11.
  L-kind: E[cos(k.x)/(2 sum(1-cos))]  = E[cos/(4(3 - s))]   (singular at k=0; computed as G(0) - Delta(x))
  S-kind: E[cos(k.x)/(2 sum(1+cos))]  = E[cos/(4(3 + s))]
  3-kind: E[cos(k.x)/(2(16 + sum cos))] = E[cos/(4(8 + s))]"""
import mpmath as mp, itertools, sys, pickle
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
P = 6
E = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
Wf = 3 * mp.gamma(mp.mpf(1) / 3)**6 / (2**(mp.mpf(14) / 3) * mp.pi**4)    # E[1/(1 - s/3)]


def base(x, t, sgn):
    """E[cos(k.x)/(t - sgn*s)] with k3 integrated exactly (t outside the band)."""
    m1, m2, m3 = [abs(v) for v in x]
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2)
        A = t - sgn * c1 * c2; B = sgn * (c1 + c2)
        root = mp.sqrt(A * A - B * B)
        if m3 == 0: e3 = 1 / root
        elif B == 0: e3 = mp.mpf(0)
        else: e3 = ((A - root) / B)**m3 / root
        return mp.cos(m1 * k1) * mp.cos(m2 * k2) * e3
    return mp.quad(f, [0, mp.pi / 2, mp.pi], [0, mp.pi / 2, mp.pi]) / mp.pi**2


def delta(x):
    """E[(1 - cos k.x)/(3 - s)] (regular)."""
    m1, m2, m3 = [abs(v) for v in x]
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2)
        A = 3 - c1 * c2; B = c1 + c2
        root = mp.sqrt(A * A - B * B)
        if root == 0: return mp.mpf(0)
        if m3 == 0: r = 1
        elif B == 0: r = mp.mpf(0)
        else: r = ((A - root) / B)**m3
        return (1 - mp.cos(m1 * k1) * mp.cos(m2 * k2) * r) / root
    return mp.quad(f, [0, mp.pi / 2, mp.pi], [0, mp.pi / 2, mp.pi]) / mp.pi**2


cache = {}
def G(kind, x):
    a = sorted(abs(v) for v in x[:2]) + [abs(x[2])]
    key = (kind, tuple(sorted(abs(v) for v in x)))     # full cubic symmetry of s
    if key in cache: return cache[key]
    xx = key[1]
    if kind == "L":
        v = ((Wf / 3) - (delta(xx) if any(xx) else 0)) / 4
    elif kind == "S":
        v = base(xx, 3, -1) / 4
    else:
        v = base(xx, 8, -1) / 4
    cache[key] = v
    return v


add = lambda a, b: tuple(p + q for p, q in zip(a, b)); neg = lambda a: tuple(-p for p in a)
Z = (0, 0, 0)
def BdB(kind, i, Ri, j, Rj):
    x = add(Rj, neg(Ri))
    return mp.fsum(G(kind, add(add(x, a), b)) for a, b in [(Z, Z), (neg(E[i]), Z), (Z, E[j]), (neg(E[i]), E[j])])
def Linv(i, Ri, j, Rj): return ((1 if (i, Ri) == (j, Rj) else 0) + BdB("L", i, Ri, j, Rj)) / (4 * P)
def Qinv(i, Ri, j, Rj): return ((1 if (i, Ri) == (j, Rj) else 0) - BdB("3", i, Ri, j, Rj)) / (4 * P - 4)
def Fm(i, Ri, j, Rj): return (1 if (i, Ri) == (j, Rj) else 0) - BdB("S", i, Ri, j, Rj)
o = (0, Z)
nbrs = []
for t in range(P):
    for R in itertools.product(range(-2, 3), repeat=3):
        s = (t, R)
        if s == o: continue
        if {R, add(R, E[t])} & {Z, E[0]}: nbrs.append(s)
assert len(nbrs) == 4 * P - 2, len(nbrs)
def Yent(Minv, a, b, signless):
    sg = 1 if signless else -1
    return Minv(*o, *o) + sg * Minv(*o, *b) + sg * Minv(*a, *o) + Minv(*a, *b)
deg = 4 * P - 2; mu0 = -mp.mpf(2) / deg; th2 = 1 - mu0**2
classes = {}
for a in nbrs:
    Yr = [Yent(Linv, a, b, False) for b in nbrs]
    key0 = mp.nstr(Yr[nbrs.index(a)], 15) + "|" + mp.nstr(mp.fsum(y * y for y in Yr), 15)
    if key0 in classes: continue
    Qr = [Yent(Qinv, a, b, True) for b in nbrs]
    pp = mp.fsum(((1 if b == a else 0) - y)**2 for b, y in zip(nbrs, Yr)) / 4
    pm = mp.fsum(((1 if b == a else 0) - y)**2 for b, y in zip(nbrs, Qr)) / 4
    fl = 0
    for sg in (1, -1):
        lam = mp.mpc(mu0, sg * mp.sqrt(th2))
        fl += mp.fsum(abs((Fm(*o, *o) - lam * Fm(*b, *o) - mp.conj(lam) * Fm(*o, *a) + Fm(*b, *a)) / (2 * th2 * deg))**2 for b in nbrs)
    classes[key0] = (a, pp, pm, fl, Yr, Qr)
    print(f"L(FCC) start {a}: TOTAL {mp.nstr(pp + pm + fl, 25)} | +1 {mp.nstr(pp, 20)} -1 {mp.nstr(pm, 20)} flat {mp.nstr(fl, 20)} | Yrowsum {mp.nstr(mp.fsum(Yr), 12)}", flush=True)
pickle.dump((nbrs, classes, dict(cache)), open("lgFCC.pkl", "wb"))
