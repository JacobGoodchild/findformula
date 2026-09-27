"""Line graph of the simple cubic lattice: flip-flop Grover trapping (Y row, flat channel)."""
import sympy as sp, mpmath as mp, sys, pickle, itertools
from gen3d import *
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 15
dirs = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
def ends(t, x): d = dirs[t]; return {x, tuple(a + b for a, b in zip(x, d))}
edges = set()
for t0 in range(3):
    for t in range(3):
        for x in itertools.product(range(-2, 3), repeat=3):
            if (t, x) == (t0, (0, 0, 0)): continue
            if ends(t, x) & ends(t0, (0, 0, 0)):
                key = (t0, t, x); rk = (t, t0, tuple(-a for a in x))
                if rk not in edges: edges.add(key)
E = sorted(edges); nv = 3
A, D, deg = bloch(nv, E); print("deg", deg, flush=True)
B = sp.Matrix([[1 + z1**-d[0] * z2**-d[1] * z3**-d[2] for d in dirs]]); Bd = sp.Matrix([[1 + z1**d[0] * z2**d[1] * z3**d[2]] for d in dirs])
if sp.simplify(sp.expand(Bd * B - (A + 2 * sp.eye(3)))) != sp.zeros(3, 3):
    B, Bd = sp.Matrix([[1 + z1**d[0] * z2**d[1] * z3**d[2] for d in dirs]]), sp.Matrix([[1 + z1**-d[0] * z2**-d[1] * z3**-d[2]] for d in dirs])
print("check", sp.simplify(sp.expand(Bd * B - (A + 2 * sp.eye(3)))) == sp.zeros(3, 3), flush=True)
arcs = out_arcs(E, 0); print(arcs, flush=True)
Y, _, _ = transfer_row(nv, E, 0, 0, False)
det = sp.expand((B * Bd)[0]); Fnum = sp.expand(det * sp.eye(3) - Bd * B); av = Avg3(det); O = (0, (0, 0, 0))
F = lambda p, q: av.avg(sp.expand(Fnum[p[0], q[0]] * ph(p[1], q[1])))
Foo = F(O, O); Foa = F(O, arcs[0]); Fbo = [F(hb, O) for hb in arcs]; Fba = [F(hb, arcs[0]) for hb in arcs]
print("F", Foo, Foa, Fbo, Fba, flush=True)
YQ, _, _ = transfer_row(nv, E, 0, 0, True)
pickle.dump((arcs, Y, YQ, Foo, Foa, Fbo, Fba), open("lgsc.pkl", "wb"))
mu0 = -mp.mpf(1)/5; th2 = 1 - mu0**2; n = len(arcs)
pp = mp.fsum(((1 if j == 0 else 0) - Y[j])**2 for j in range(n)) / 4
pm = mp.fsum(((1 if j == 0 else 0) - YQ[j])**2 for j in range(n)) / 4
fl = 0
for sg in (1, -1):
    lam = mp.mpc(mu0, sg * mp.sqrt(th2))
    fl += mp.fsum(abs((Foo - lam * Fbo[j] - mp.conj(lam) * Foa + Fba[j]) / (2 * th2 * 10))**2 for j in range(n))
print("+1", pp, "-1", pm, "flat", fl, "TOTAL", pp + pm + fl)
