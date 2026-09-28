"""Compute Y, YQ once; then channels for every start arc (incl. flat bands). Save to pickle."""
import mpmath as mp, sys, pickle, time
from gen import transfer_matrix, flat_bands, flat_projector_channel, out_arcs
from builders import LATTICES
import sympy as sp
mp.mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 25
name, origin = sys.argv[1], int(sys.argv[2])
bip = name in ("lieb", "dice", "truncated_square(4.8.8)", "truncated_trihexagonal(4.6.12)")
nv, E = LATTICES[name]()
t = time.time()
Y, arcs, deg = transfer_matrix(nv, E, origin, False)
YQ = Y if bip else transfer_matrix(nv, E, origin, True)[0]
fl, H, _ = flat_bands(nv, E)
out = {"arcs": arcs, "Y": Y, "YQ": YQ, "flat": {}, "bip": bip}
n = len(arcs)
for s in range(n):
    pp = mp.fsum(((1 if j == s else 0) - Y[s][j])**2 for j in range(n)) / 4
    pm = mp.fsum(((1 if j == s else 0) - YQ[s][j])**2 for j in range(n)) / 4
    ex = {}
    for (m0, mult) in fl:
        if abs(sp.N(m0)) < 1:
            ex[str(m0)] = flat_projector_channel(nv, E, origin, m0, H, deg, s)
    tot = pp + pm + mp.fsum(ex.values())
    out["flat"][s] = (pp, pm, ex, tot)
    print(name, "origin", origin, "start", s, arcs[s], "TOTAL", mp.nstr(tot, 22), "| +1", mp.nstr(pp, 18), "-1", mp.nstr(pm, 18), {k: mp.nstr(v, 18) for k, v in ex.items()}, flush=True)
pickle.dump(out, open(f"full_{name.split('(')[0]}_{origin}.pkl", "wb"))
print("time", round(time.time() - t))
