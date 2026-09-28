"""Checkerboard (planar pyrochlore) lattice: square lattice + crossed diagonals (conductance w) on
alternate plaquettes. Spanning-tree entropy per vertex via Jensen's formula on the 2x2 supercell."""
import sympy as sp, mpmath as mp, sys
from gen import z1, z2
from genus_scan import supercell, lap_det
W = sp.symbols('W', positive=True)
rule = lambda x, y: [(1, 0, 1), (0, 1, 1)] + ([(1, 1, W)] if (x + y) % 2 == 0 else [(-1, 1, W)])
nv, E = supercell(rule)
DET = lap_det(nv, E)
def Elogdet(w, dps=40):
    mp.mp.dps = dps
    d = sp.expand(DET.subs(W, sp.Rational(str(w)) if not isinstance(w, sp.Basic) else w))
    pw = [t.as_powers_dict().get(z2, 0) for t in sp.Add.make_args(d)]
    P = sp.Poly(sp.expand(d * z2**(-min(pw))), z2)
    co = [sp.lambdify(z1, c, "mpmath") for c in P.all_coeffs()]
    def f(k):
        with mp.extradps(dps):
            pc = [mp.mpc(c(mp.expj(k))) for c in co]
            big = max(abs(c) for c in pc)
            while abs(pc[0]) <= big * mp.mpf(10)**(-mp.mp.dps): pc = pc[1:]
            r = mp.polyroots(pc, maxsteps=400, extraprec=4 * mp.mp.prec)
            return mp.log(abs(pc[0])) + mp.fsum(mp.log(abs(x)) for x in r if abs(x) > 1)
    return mp.quad(f, [j * mp.pi / 12 for j in range(-12, 13)]) / (2 * mp.pi)
if __name__ == "__main__":
    for w in sys.argv[1:]:
        v = Elogdet(w); print(w, v, v / 4)
