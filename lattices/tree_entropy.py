"""Spanning-tree entropy (log Mahler measure of the Laplacian determinant) of a 2D periodic graph:
   E_k[ln det L(k)] = (1/2pi) INT dk1 [ ln|lead(z1)| + SUM_{roots r(z1) of det in z2, |r|>1} ln|r| ]   (Jensen)
Usage: python3 tree_entropy.py <lattice name> [dps]"""
import sys, sympy as sp, mpmath as mp
from gen import bloch_mats, z1, z2
from builders import LATTICES
def entropy(nv, E, dps=40, signless=False):
    mp.mp.dps = dps
    A, D, deg = bloch_mats(nv, E)
    d = sp.expand((D + A if signless else D - A).det())
    pw = [t.as_powers_dict().get(z2, 0) for t in sp.Add.make_args(d)]
    P = sp.Poly(sp.expand(d * z2**(-min(pw))), z2)
    co = [sp.lambdify(z1, c, "mpmath") for c in P.all_coeffs()]
    def f(k):
        with mp.extradps(dps):
            zz = mp.expj(k)
            pc = [mp.mpc(c(zz)) for c in co]
            big = max(abs(c) for c in pc)
            while abs(pc[0]) <= big * mp.mpf(10)**(-mp.mp.dps): pc = pc[1:]
            roots = mp.polyroots(pc, maxsteps=400, extraprec=4 * mp.mp.prec)
            return mp.log(abs(pc[0])) + mp.fsum(mp.log(abs(r)) for r in roots if abs(r) > 1)
    pts = [j * mp.pi / 12 for j in range(-12, 13)]
    return mp.quad(f, pts) / (2 * mp.pi), nv
if __name__ == "__main__":
    name = sys.argv[1]; dps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    nv, E = LATTICES[name]()
    v, nv = entropy(nv, E, dps)
    print(name, "E ln det =", mp.nstr(v, dps - 5), " per vertex z =", mp.nstr(v / nv, dps - 5))
