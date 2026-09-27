import sympy as sp, mpmath as mp, pickle
from gen import Averager, transfer_matrix, out_arcs, phase, z1, z2
from builders import LATTICES
mp.mp.dps = 25
dirs = [(1, 0), (0, 1), (1, -1)]
nv, E = LATTICES["linegraph_triangular"]()
# flat band: A + 2 = B^dag B, B = [1 + z^delta_t]
B = sp.Matrix([[1 + z1**d[0] * z2**d[1] for d in dirs]])
Bd = sp.Matrix([[1 + z1**-d[0] * z2**-d[1]] for d in dirs])
from gen import bloch_mats
A, D, deg = bloch_mats(nv, E)
ok1 = sp.simplify(sp.expand(Bd * B - (A + 2 * sp.eye(3)))) == sp.zeros(3, 3)
Bd2 = sp.Matrix([[1 + z1**d[0] * z2**d[1]] for d in dirs]); B2 = sp.Matrix([[1 + z1**-d[0] * z2**-d[1] for d in dirs]])
ok2 = sp.simplify(sp.expand(Bd2 * B2 - (A + 2 * sp.eye(3)))) == sp.zeros(3, 3)
print("B^dag B == A+2:", ok1, ok2)
if not ok1: B, Bd = B2, Bd2
det = sp.expand((B * Bd)[0])
Fnum = sp.expand(det * sp.eye(3) - Bd * B)
av = Averager(det)
arcs = out_arcs(E, 0); O = (0, (0, 0))
def F(p, q): return av.avg(sp.expand(Fnum[p[0], q[0]] * phase(p[1], q[1])))
Y, _, _ = transfer_matrix(nv, E, 0, False); YQ, _, _ = transfer_matrix(nv, E, 0, True)
res = {}
mu0 = -mp.mpf(1) / 5; th2 = 1 - mu0**2; dg = 10
for s, ha in [(i, arcs[i]) for i in (0, 2, 3)]:
    pp = mp.fsum(((1 if j == s else 0) - Y[s][j])**2 for j in range(len(arcs))) / 4
    pm = mp.fsum(((1 if j == s else 0) - YQ[s][j])**2 for j in range(len(arcs))) / 4
    Foo = F(O, O); Foa = F(O, ha); fl = 0
    Fb = [(F(hb, O), F(hb, ha)) for hb in arcs]
    for sg in (1, -1):
        lam = mp.mpc(mu0, sg * mp.sqrt(th2))
        fl += mp.fsum(abs((Foo - lam * fb0 - mp.conj(lam) * Foa + fba) / (2 * th2 * dg))**2 for fb0, fba in Fb)
    res[s] = (pp, pm, fl, pp + pm + fl, [x[1] for x in Fb], Foo, Foa)
    print("start", s, arcs[s], "TOTAL", mp.nstr(pp + pm + fl, 20), "| +1", mp.nstr(pp, 16), "-1", mp.nstr(pm, 16), "flat", mp.nstr(fl, 16), flush=True)
    pass
pickle.dump((arcs, Y, YQ, res), open("lgtri.pkl", "wb"))
