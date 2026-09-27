"""Flip-flop Grover trapping on the LINE GRAPH of the hypercubic lattice Z^d (d >= 3), reduced
exactly to base-lattice Green functions (fast Bessel integrals):
  B = [1 + z^{e_t}]_t (1 x d),  A_LG + 2 = B^dag B,  degree 4d - 2
  L^{-1} = (1/4d)[I + B^dag B/(4d - BB^dag)],        4d - BB^dag     = 2 sum(1 - cos k)
  Q^{-1} = (1/(4d-4))[I - B^dag B/(4d-4 + BB^dag)],  4d - 4 + BB^dag = 2(3d - 2 + sum cos k)
  F      = I - B^dag B/(BB^dag),                     BB^dag          = 2 sum(1 + cos k)
Sites of L(Z^d): (t, R) = the edge from R to R + e_t."""
import mpmath as mp, itertools, sys, pickle

d = int(sys.argv[1]); mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
P = d


def bessel_G(x, a, sgn):
    """E[cos(k.x) / (a + sgn * sum cos k_i)] = int_0^inf e^{-a t} prod_i s_i I_{x_i}(t) dt,
    s_i = 1 for sgn = -1 and (-1)^{x_i} for sgn = +1 (valid when the integral converges)."""
    def f(t):
        return mp.exp(-a * t) * mp.fprod((1 if sgn == -1 else (-1)**abs(xi)) * mp.besseli(abs(xi), t) for xi in x)
    return mp.quad(f, [0, 1, 5, 20, 80, 320, 1280, mp.inf])


cache = {}


def G(kind, x):
    x = tuple(sorted(abs(v) for v in x))
    key = (kind, x)
    if key not in cache:
        if kind == "L":
            cache[key] = bessel_G(x, d, -1) / 2            # E[cos/(2 sum(1 - cos))]
        elif kind == "S":
            cache[key] = bessel_G(x, d, +1) / 2            # E[cos/(2 sum(1 + cos))]
        else:
            cache[key] = bessel_G(x, 3 * P - 2, +1) / 2    # E[cos/(2(3P - 2 + sum cos))]
    return cache[key]


E = [tuple(1 if i == t else 0 for i in range(d)) for t in range(d)]
add = lambda a, b: tuple(x + y for x, y in zip(a, b))
neg = lambda a: tuple(-x for x in a)
Z = tuple([0] * d)


def BdB(kind, i, Ri, j, Rj):
    """E[(B^dag B)_{ij} e^{ik.(Ri - Rj)} / den],  (B^dag B)_{ij} = (1 + z^{-e_i})(1 + z^{e_j})."""
    x = add(Rj, neg(Ri))
    return mp.fsum(G(kind, add(add(x, a), b)) for a, b in [(Z, Z), (neg(E[i]), Z), (Z, E[j]), (neg(E[i]), E[j])])


def Linv(i, Ri, j, Rj): return ((1 if (i, Ri) == (j, Rj) else 0) + BdB("L", i, Ri, j, Rj)) / (4 * P)
def Qinv(i, Ri, j, Rj): return ((1 if (i, Ri) == (j, Rj) else 0) - BdB("3", i, Ri, j, Rj)) / (4 * P - 4)
def Fm(i, Ri, j, Rj): return (1 if (i, Ri) == (j, Rj) else 0) - BdB("S", i, Ri, j, Rj)


o = (0, Z)
nbrs = []
for t in range(d):
    for R in itertools.product(range(-1, 2), repeat=d):
        s = (t, R)
        if s == o:
            continue
        if {R, add(R, E[t])} & {Z, E[0]}:
            nbrs.append(s)
assert len(nbrs) == 4 * P - 2, len(nbrs)


def Yent(Minv, a, b, signless):
    sg = 1 if signless else -1
    return Minv(*o, *o) + sg * Minv(*o, *b) + sg * Minv(*a, *o) + Minv(*a, *b)


deg = 4 * P - 2
mu0 = -mp.mpf(2) / deg
th2 = 1 - mu0**2
classes = {}
for a in nbrs:
    Yr = [Yent(Linv, a, b, False) for b in nbrs]
    Qr = [Yent(Qinv, a, b, True) for b in nbrs]
    pp = mp.fsum(((1 if b == a else 0) - y)**2 for b, y in zip(nbrs, Yr)) / 4
    pm = mp.fsum(((1 if b == a else 0) - y)**2 for b, y in zip(nbrs, Qr)) / 4
    fl = 0
    for sg in (1, -1):
        lam = mp.mpc(mu0, sg * mp.sqrt(th2))
        fl += mp.fsum(abs((Fm(*o, *o) - lam * Fm(*b, *o) - mp.conj(lam) * Fm(*o, *a) + Fm(*b, *a)) / (2 * th2 * deg))**2 for b in nbrs)
    key = mp.nstr(pp + pm + fl, 20)
    if key not in classes:
        classes[key] = (a, pp, pm, fl, Yr, Qr)
        print(f"d={d} start {a}: TOTAL {mp.nstr(pp + pm + fl, 25)} | +1 {mp.nstr(pp, 20)} -1 {mp.nstr(pm, 20)} "
              f"flat {mp.nstr(fl, 20)} | Y row sum {mp.nstr(mp.fsum(Yr), 12)} | Y[a][a] {mp.nstr(Yr[nbrs.index(a)], 20)}", flush=True)
pickle.dump((nbrs, classes, dict(cache)), open(f"lgZ{d}.pkl", "wb"))
