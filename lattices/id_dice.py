import mpmath as mp
from gen import trapping
from builders import LATTICES
from ident import lin
mp.mp.dps = 30; pi = mp.pi; s3 = mp.sqrt(3)
nv, E = LATTICES["dice"]()
tot, pp, pm, extra, Y, YQ = trapping(nv, E, origin=0, start=0, bipartite=True, verbose=False)
fl = list(extra.values())[0]
B = {"1": 1, "s3/pi": s3/pi, "1/pi": 1/pi, "1/pi^2": 1/pi**2, "s3/pi^2": s3/pi**2}
print("flat", fl, lin(fl, B))
print("total", tot, lin(tot, B))
