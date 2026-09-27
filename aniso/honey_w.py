"""Anisotropic honeycomb Szegedy walk: bond probabilities p1, p2, p3 (sum 1).
Weighted transfer currents Y(a,b) = sqrt(p_a p_b) (delta_o - delta_ha)^T L_P^{-1} (delta_o - delta_hb)."""
import sys, sympy as sp, mpmath as mp
sys.path.insert(0, '../lattices')
from gen import Averager, z1, z2
mp.mp.dps = 25
def honey(p):
    p = [sp.nsimplify(x) for x in p]
    offs = [(0, 0), (-1, 0), (0, -1)]     # B neighbours of A(0)
    f = sum(pi * z1**o[0] * z2**o[1] for pi, o in zip(p, offs))
    fb = sum(pi * z1**-o[0] * z2**-o[1] for pi, o in zip(p, offs))
    L = sp.Matrix([[1, -f], [-fb, 1]])
    det = sp.expand(L.det()); adj = L.adjugate()
    av = Averager(det)
    O = (0, (0, 0)); heads = [(1, o) for o in offs]
    def G(pq, qq):
        (pv, po), (qv, qo) = pq, qq
        return adj[pv, qv] * z1**(po[0] - qo[0]) * z2**(po[1] - qo[1])
    Y = [[None] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            N = G(O, O) - G(O, heads[j]) - G(heads[i], O) + G(heads[i], heads[j])
            Y[i][j] = mp.sqrt(mp.mpf(sp.N(p[i] * p[j], 40))) * av.avg(sp.expand(N))
    return Y
if __name__ == "__main__":
    for p in [(sp.Rational(1, 3),) * 3, (sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 4)), (sp.Rational(2, 5), sp.Rational(1, 3), sp.Rational(4, 15)), (sp.Rational(3, 5), sp.Rational(1, 5), sp.Rational(1, 5))]:
        Y = honey(p)
        pf = [float(x) for x in p]
        # triangle angles opposite sides p_i
        import math
        ang = []
        for i in range(3):
            a, b, c = pf[i], pf[(i + 1) % 3], pf[(i + 2) % 3]
            cosv = (b * b + c * c - a * a) / (2 * b * c)
            ang.append(math.acos(max(-1, min(1, cosv))))
        print("p =", p, " Y diag:", [mp.nstr(Y[i][i], 15) for i in range(3)], " 2theta/pi:", [round(2 * t / math.pi, 15) for t in ang])
        print("    Y row0:", [mp.nstr(y, 15) for y in Y[0]])
        pbar = 2 * sum(((1 if j == 0 else 0) - Y[0][j])**2 for j in range(3)) / 4
        print("    pbar(start type 1) =", pbar)
