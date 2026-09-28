"""Tree entropies of the genus-0 plaquette decorations found by plaquette_search.py, with PSLQ identification."""
import json, mpmath as mp, sympy as sp, sys
from gen import z1, z2
from plaquette_search import build
from genus_scan2 import lapdet
def entropy_from_det(d, dps):
    mp.mp.dps = dps
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
    res = json.load(open("plaquette_search.json"))
    dps = 45
    for c, deg, odd, flag in res:
        if not flag: continue
        nv, E = build(c)
        v = entropy_from_det(lapdet(nv, E), dps)
        mp.mp.dps = dps
        pi = mp.pi; s2, s3, s5 = mp.sqrt(2), mp.sqrt(3), mp.sqrt(5)
        names = ["G/pi", "ln2", "ln3", "ln(1+s2)", "ln(2+s3)", "ln(s2+s3)", "ln5", "ln phi", "s3Cl(pi/3)/pi", "Cl(pi/4)/pi"]
        vals = [mp.catalan / pi, mp.log(2), mp.log(3), mp.log(1 + s2), mp.log(2 + s3), mp.log(s2 + s3), mp.log(5), mp.log((1 + s5) / 2), s3 * mp.clsin(2, pi / 3) / pi, mp.clsin(2, pi / 4) / pi]
        r = mp.pslq([v] + vals, maxcoeff=10**4, maxsteps=10**7)
        rel = " ".join(f"{-k}/{r[0]}*{n}" for k, n in zip(r[1:], names) if k) if r and r[0] else None
        print(c.replace("0", "."), "E ln det(4 sites) =", mp.nstr(v, 30), "|", rel, flush=True)
