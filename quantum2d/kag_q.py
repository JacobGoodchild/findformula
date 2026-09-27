import mpmath as mp, pickle
from identify import find_lin
mp.mp.dps = 35; pi = mp.pi; s2 = mp.sqrt(2)
_, _, _, _, YQk = pickle.load(open("hk.pkl", "rb"))
m = mp.mpf(5)/32; K = mp.ellipk(m); E = mp.ellipe(m)
names = ["1", "K/pi", "E/pi", "s2 K/pi", "s2 E/pi", "pi/K", "s2 pi/K"]
vals = [1, K/pi, E/pi, s2*K/pi, s2*E/pi, pi/K, s2*pi/K]
for j in range(4):
    print("YQ0%d" % j, YQk[0][j], find_lin(YQk[0][j], names, vals, maxc=10**4))
