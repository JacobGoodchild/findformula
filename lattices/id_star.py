import mpmath as mp, pickle
from ident import lin
from fractions import Fraction
mp.mp.dps = 25; pi = mp.pi
d = pickle.load(open("full_star_0.pkl", "rb"))
for m in [mp.mpf(3)/128, 1 - mp.mpf(3)/128, mp.mpf(3)/32, mp.mpf(1)/9, mp.mpf(1)/16]:
    K = mp.ellipk(m); E = mp.ellipe(m)
    B = {"1": 1, "K/pi": K/pi, "E/pi": E/pi}
    rows = [[lin(y, B, maxc=5000) for y in d["YQ"][s]] for s in (0, 2)]
    ok = all(r is not None for row in rows for r in row)
    print("m =", m, "all identified:", ok)
    if ok:
        for s, row in zip((0, 2), rows): print("  start", s, row)
        break
for s in (0, 2):
    pp, pm, ex, tot = d["flat"][s]
    print("start", s, "+1 =", mp.identify(pp), " flats:", {k: mp.identify(v) for k, v in ex.items()})
