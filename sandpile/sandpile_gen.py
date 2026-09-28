"""Height-1 probability of the abelian sandpile at vertex `origin` of a general 2D periodic lattice.
Resistances between the origin and its neighbours (and among neighbours) come from the exact Bloch
engine (lattices/gen.py); P1 = det(I - M), M_jk = (R(0,j) + R(0,k) - R(j,k))/2 over cut neighbours."""
import sys, sympy as sp, mpmath as mp
sys.path.insert(0, '/home/user/findformula/lattices')
from gen import bloch_mats, Averager, out_arcs, z1, z2, phase
from builders import LATTICES
def resistances(nv, E, origin, dps=25):
    mp.mp.dps = dps
    A, D, deg = bloch_mats(nv, E)
    L = D - A
    det = sp.expand(L.det()); adj = L.adjugate(); av = Averager(det)
    O = (origin, (0, 0))
    nbrs = [O] + [(v, o) for (v, o) in out_arcs(E, origin)]
    def G(p, q): return adj[p[0], q[0]] * phase(p[1], q[1])
    R = {}
    for i, p in enumerate(nbrs):
        for j, q in enumerate(nbrs):
            if j <= i: continue
            N = G(p, p) + G(q, q) - G(p, q) - G(q, p)
            R[(i, j)] = R[(j, i)] = av.avg(N)
    for i in range(len(nbrs)): R[(i, i)] = mp.mpf(0)
    return nbrs, R
def P1_numeric(nbrs, R, keep=1):
    cut = [i for i in range(1, len(nbrs)) if i != keep]
    n = len(cut); M = mp.matrix(n, n)
    for a, j in enumerate(cut):
        for b, k in enumerate(cut):
            M[a, b] = (R[(0, j)] + R[(0, k)] - R[(j, k)]) / 2
    return mp.det(mp.eye(n) - M)
if __name__ == "__main__":
    name, origin = sys.argv[1], int(sys.argv[2])
    nv, E = LATTICES[name]()
    nbrs, R = resistances(nv, E, origin)
    print(name, origin, nbrs)
    print("P1 (keep 1st bond) =", P1_numeric(nbrs, R, 1), "  (keep last) =", P1_numeric(nbrs, R, len(nbrs) - 1))
    import pickle; pickle.dump((nbrs, R), open(f"R_{name.split('(')[0]}_{origin}.pkl", "wb"))
