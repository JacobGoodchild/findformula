"""Exact P1 from pickled resistances: PSLQ each R in a lattice-specific basis, then symbolic determinant."""
import sys, pickle, sympy as sp, mpmath as mp
BASES = {
 "dice": (["1", "s3/pi"], lambda: [1, mp.sqrt(3) / mp.pi]),
 "lieb": (["1", "1/pi"], lambda: [1, 1 / mp.pi]),
 "star": (["1", "s3/pi"], lambda: [1, mp.sqrt(3) / mp.pi]),
 "truncated_square": (["1", "1/pi", "c"], lambda: [1, 1 / mp.pi, mp.sqrt(2) * mp.acos(mp.mpf(1) / 3) / mp.pi]),
 "union_jack": (["1", "1/pi", "c"], lambda: [1, 1 / mp.pi, mp.sqrt(2) * mp.acos(mp.mpf(1) / 3) / mp.pi]),
 "shastry_sutherland": (["1", "s3", "1/pi", "s3/pi"], lambda: [1, mp.sqrt(3), 1 / mp.pi, mp.sqrt(3) / mp.pi]),
}
if __name__ == "__main__":
    pkl, lat = sys.argv[1], sys.argv[2]
    nbrs, R = pickle.load(open(pkl, "rb"))
    mp.mp.dps = 25
    names, bf = BASES[lat]; B = bf()
    syms = [sp.Integer(1)] + [sp.Symbol(n) for n in names[1:]]
    RS = {}
    for (i, j), v in R.items():
        if i >= j: continue
        rel = mp.pslq([v] + B, maxcoeff=10**5, maxsteps=10**6)
        if not rel or rel[0] == 0: print("no relation for", i, j, v); sys.exit(1)
        RS[(i, j)] = RS[(j, i)] = sum(sp.Rational(-c, rel[0]) * s for c, s in zip(rel[1:], syms))
    for i in range(len(nbrs)): RS[(i, i)] = 0
    for i in range(1, len(nbrs)): print("R(0,%d) =" % i, RS[(0, i)])
    for keep in (1, len(nbrs) - 1):
        cut = [i for i in range(1, len(nbrs)) if i != keep]
        M = sp.Matrix(len(cut), len(cut), lambda a, b: (RS[(0, cut[a])] + RS[(0, cut[b])] - RS[(cut[a], cut[b])]) / 2)
        P = sp.factor(sp.expand((sp.eye(len(cut)) - M).det()))
        subs = {sp.Symbol(n): v for n, v in zip(names[1:], B[1:])}
        print("keep", keep, ": P1 =", P, " =", sp.N(P.subs({k: sp.Float(str(v), 30) for k, v in subs.items()}), 25))
