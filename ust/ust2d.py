import sys, pickle, sympy as sp, mpmath as mp
sys.path.insert(0, '/home/user/findformula/sandpile')
from sandpile_exact import BASES
from ustdeg import degree_distribution
mp.mp.dps = 25
for pkl, lat in [("R_triangular_0.pkl", "triangular"), ("R_kagome_0.pkl", "kagome"), ("R_truncated_square_0.pkl", "truncated_square"), ("R_union_jack_0.pkl", "union_jack"), ("R_union_jack_1.pkl", "union_jack"),
                 ("R_star_0.pkl", "star"), ("R_dice_0.pkl", "dice"), ("R_dice_1.pkl", "dice"), ("R_lieb_0.pkl", "lieb"),
                 ("R_checkerboard_0.pkl", "checkerboard"), ("R_triakis_triangular_0.pkl", "triakis_triangular"), ("R_triakis_triangular_1.pkl", "triakis_triangular")]:
    nbrs, R = pickle.load(open('/home/user/findformula/sandpile/' + pkl, 'rb'))
    names, bf = BASES[lat]; B = bf()
    syms = [sp.Integer(1)] + [sp.Symbol(n.replace('/', '_').replace('s3', 'r')) for n in names[1:]]
    RS = {}
    for (i, j), v in R.items():
        if i >= j: continue
        rel = mp.pslq([v] + B, maxcoeff=10**5, maxsteps=10**6)
        RS[(i, j)] = RS[(j, i)] = sum(sp.Rational(-c, rel[0]) * s for c, s in zip(rel[1:], syms))
    for i in range(len(nbrs)): RS[(i, i)] = 0
    z = len(nbrs) - 1
    Y = sp.Matrix(z, z, lambda a, b: (RS[(0, a + 1)] + RS[(0, b + 1)] - RS[(a + 1, b + 1)]) / 2)
    P = degree_distribution(Y)
    sub = {s: sp.Float(str(v), 25) for s, v in zip(syms[1:], B[1:])}
    vals = [sp.N(p.subs(sub), 15) for p in P]
    print(pkl[2:-4], "z=%d" % z, "sum", sp.N(sum(vals), 10), "mean", sp.N(sum(k * v for k, v in enumerate(vals)), 10))
    for k, p in enumerate(P):
        if p != 0: print("   deg", k, ":", p, " =", vals[k])
