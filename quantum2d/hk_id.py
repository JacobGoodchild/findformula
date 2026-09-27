import mpmath as mp, pickle
from lgf2d import b
from identify import find_gamma, find_lin
mp.mp.dps = 35; pi = mp.pi; s3 = mp.sqrt(3)
(hp, _, _), Yh, (kp, kpp, kpm), Yk, YQk = pickle.load(open("hk.pkl", "rb"))
print("honeycomb Y01:", find_lin(Yh[0][1], ["1", "s3/pi"], [1, s3/pi]))
for j in range(4):
    print("kagome Y0%d:" % j, find_lin(Yk[0][j], ["1", "s3/pi"], [1, s3/pi]))
g11 = b((0, 0), lambda k: 11 - mp.cos(k), lambda k: 1 + mp.exp(-1j*k))
print("g(11) = E[1/(11 - c1 - c2 - c12)] =", g11)
print(" gamma search:", find_gamma(g11))
for m in [0]:
    pass
# elliptic: triangular LGF G(E) = ... K(k) with k^2 from E; try K at candidate moduli
for j in range(4):
    print("kagome YQ0%d:" % j, YQk[0][j], find_lin(YQk[0][j], ["1", "g", "1/(pi^2 g)", "s3/pi", "1/pi"], [1, g11, 1/(pi**2*g11), s3/pi, 1/pi]))
