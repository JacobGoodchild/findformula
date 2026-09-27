import mpmath as mp, pickle
from ident import lin
mp.mp.dps = 25; pi = mp.pi
d = pickle.load(open("full_star_0.pkl", "rb"))
m = mp.mpf(3)/128; K = mp.ellipk(m); E = mp.ellipe(m); s6 = mp.sqrt(6)
B = {"1": 1, "s6K/pi": s6*K/pi, "s6E/pi": s6*E/pi}
for s in (0, 1, 2):
    print("start", s, [lin(y, B, maxc=20000) for y in d["YQ"][s]])
