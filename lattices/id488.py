import mpmath as mp
from gen import trapping
from builders import LATTICES
from ident import lin
mp.mp.dps = 30; pi = mp.pi; s2 = mp.sqrt(2)
ac = mp.acos(mp.mpf(1)/3)
B = {"1": 1, "1/pi": 1/pi, "s2/pi": s2/pi, "ac/pi": ac/pi, "s2ac/pi": s2*ac/pi}
nv, E = LATTICES["truncated_square(4.8.8)"]()
for start in range(3):
    tot, pp, pm, extra, Y, YQ = trapping(nv, E, origin=0, start=start, bipartite=True, verbose=False)
    print("start", start, "total", tot)
    print("  Y row:", [lin(y, B) for y in Y[start]])
B2 = dict(B); B2.update({"1/pi^2": 1/pi**2, "s2/pi^2": s2/pi**2, "ac/pi^2": ac/pi**2, "ac^2/pi^2": ac**2/pi**2, "s2ac/pi^2": s2*ac/pi**2})
print("total (start 0):", lin(trapping(nv, E, 0, 0, True, False)[0], B2))
