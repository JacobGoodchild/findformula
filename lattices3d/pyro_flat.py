import sympy as sp, mpmath as mp, sys, pickle
from gen3d import *
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 15
fcc = [(0, 2, 2), (2, 0, 2), (2, 2, 0)]
nv, E = geometric3(fcc, [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)], 2)
A, D, deg = bloch(nv, E)
cj = {z1: 1/z1, z2: 1/z2, z3: 1/z3}
for Brow in ([1, 1/z1, 1/z2, 1/z3], [1, z1, z2, z3]):
    B = sp.Matrix([[1, 1, 1, 1], Brow])
    Bd = B.T.subs(cj, simultaneous=True)          # conjugate transpose on |z|=1
    diff = sp.simplify(sp.expand(Bd * B - (A + 2 * sp.eye(4))))
    if diff == sp.zeros(4, 4): break
    Bd = B.T; Bc = B.subs(cj, simultaneous=True)
    diff2 = sp.simplify(sp.expand(Bc.T * B - (A + 2 * sp.eye(4))))
print("B^dag B == A + 2 :", diff == sp.zeros(4, 4))
BBd = sp.expand(B * Bd)
det = sp.expand(BBd.det()); adj = BBd.adjugate()
print("det(B B^dag) =", sp.factor(det))
Fnum = sp.expand(det * sp.eye(4) - Bd * adj * B)       # F = Fnum / det
av = Avg3(det)
arcs = out_arcs(E, 0); O = (0, (0, 0, 0))
def F(p, q): return av.avg(sp.expand(Fnum[p[0], q[0]] * ph(p[1], q[1])))
ha = arcs[0]
vals = {"oo": F(O, O), "oa": F(O, ha)}
print("F(o,o) =", vals["oo"], " F(o,h_a) =", vals["oa"], flush=True)
for hb in arcs:
    vals[("bo", hb)] = F(hb, O); vals[("ba", hb)] = F(hb, ha)
    print("  ", hb, vals[("bo", hb)], vals[("ba", hb)], flush=True)
pickle.dump((arcs, vals), open("pyro_flat.pkl", "wb"))
mu0 = -mp.mpf(1)/3; th2 = 1 - mu0**2
tot = 0
for sgn in (1, -1):
    lam = mp.mpc(mu0, sgn * mp.sqrt(th2)); s = 0
    for hb in arcs:
        v = (vals["oo"] - lam * vals[("bo", hb)] - mp.conj(lam) * vals["oa"] + vals[("ba", hb)]) / (2 * th2 * 6)
        s += abs(v)**2
    tot += s
print("flat channels total:", tot)
