import sympy as sp, mpmath as mp, pickle
from gen import Averager, out_arcs, phase, z1, z2
from builders import LATTICES
from ident import lin
mp.mp.dps = 25; pi = mp.pi; s3 = mp.sqrt(3)
dirs = [(1, 0), (0, 1), (1, -1)]
nv, E = LATTICES["linegraph_triangular"]()
Bd = sp.Matrix([[1 + z1**d[0] * z2**d[1]] for d in dirs]); B = sp.Matrix([[1 + z1**-d[0] * z2**-d[1] for d in dirs]])
det = sp.expand((B * Bd)[0]); Fnum = sp.expand(det * sp.eye(3) - Bd * B); av = Averager(det)
arcs = out_arcs(E, 0); O = (0, (0, 0))
X = mp.cbrt(4)*mp.gamma(mp.mpf(1)/3)**3/(8*pi**2)
vals = []
for hb in arcs:
    v = av.avg(sp.expand(Fnum[hb[0], O[0]] * phase(hb[1], O[1])))
    vals.append(v); print(hb, v, lin(v, {'1': 1, 'X': X, 's3/(piX)': s3/(pi*X)}), flush=True)
pickle.dump(vals, open('lgtri_Fbo.pkl', 'wb'))
