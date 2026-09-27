import mpmath as mp
from gen import trapping
from builders import LATTICES
from ident import lin
mp.mp.dps = 30; pi = mp.pi
nv, E = LATTICES["checkerboard"]()
tot, pp, pm, extra, Y, YQ = trapping(nv, E, origin=0, start=0, bipartite=False, verbose=False)
fl = list(extra.values())[0]
Bs = {"1": 1, "1/pi": 1/pi, "1/pi^2": 1/pi**2}
print("+1", pp, lin(pp, Bs)); print("flat", fl, lin(fl, Bs))
m = mp.mpf(1)/4; K = mp.ellipk(m); Ee = mp.ellipe(m)
BK = {"1": 1, "K/pi": K/pi, "E/pi": Ee/pi, "K^2/pi^2": (K/pi)**2, "E^2/pi^2": (Ee/pi)**2, "KE/pi^2": K*Ee/pi**2}
print("-1", pm, lin(pm, BK))
print("YQ row", [lin(y, {"1": 1, "K/pi": K/pi, "E/pi": Ee/pi}) for y in YQ[0]])
print("Y row", [lin(y, Bs) for y in Y[0]])
print("TOTAL", tot)
