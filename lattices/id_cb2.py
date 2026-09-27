import mpmath as mp, sympy as sp
from gen import trapping, out_arcs
from builders import LATTICES
from ident import lin
mp.mp.dps = 30; pi = mp.pi
nv, E = LATTICES["checkerboard"]()
arcs = out_arcs(E, 0); print(arcs)
st = [i for i, a in enumerate(arcs) if a[0] == 1][0]
tot, pp, pm, extra, Y, YQ = trapping(nv, E, origin=0, start=st, bipartite=False, verbose=False)
fl = list(extra.values())[0]
Bs = {"1": 1, "1/pi": 1/pi, "1/pi^2": 1/pi**2}
m = mp.mpf(1)/4; K = mp.ellipk(m); Ee = mp.ellipe(m)
print("start", st, arcs[st], "total", tot)
print("+1", lin(pp, Bs)); print("flat", lin(fl, Bs))
print("Y row", [lin(y, Bs) for y in Y[st]])
print("YQ row", [lin(y, {"1": 1, "K/pi": K/pi, "E/pi": Ee/pi}) for y in YQ[st]])
