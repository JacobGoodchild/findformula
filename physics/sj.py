import mpmath as mp
mp.mp.dps = 60
# Stimson-Jeffery drag factor for two equal touching spheres (limit alpha->0 of the bispherical series)
f = lambda x: 1 - 2 * (mp.sinh(x) ** 2 - x ** 2) / (mp.sinh(2 * x) + 2 * x)
lam = mp.quad(f, [0, 1, 4, 16, 64, 200]) / 3
print(lam)
# cross-check against the finite-alpha series at small alpha
def series(a, N):
    s = 0
    for n in range(1, N):
        num = 4 * mp.sinh((n + 0.5) * a) ** 2 - (2 * n + 1) ** 2 * mp.sinh(a) ** 2
        den = 2 * mp.sinh((2 * n + 1) * a) + (2 * n + 1) * mp.sinh(2 * a)
        s += mp.mpf(n * (n + 1)) / ((2 * n - 1) * (2 * n + 3)) * (1 - num / den)
    return 4 * mp.sinh(a) / 3 * s
mp.mp.dps = 20
for a in [1.0, 0.1, 0.02]:
    print(a, series(mp.mpf(a), int(60 / a)))
mp.mp.dps = 60
open("sj.txt", "w").write(str(lam))
