import mpmath as mp
mp.mp.dps = 25
def pbar(p1, p2, p3):
    p = [mp.mpf(x) for x in (p1, p2, p3)]
    lam = mp.sqrt(1 / (p[0] * p[1] * p[2]))
    t = [mp.atan(lam * x) / mp.pi for x in p]
    T12 = (mp.mpf(1)/2 - p[2] - p[0]*t[0] - p[1]*t[1] + p[2]*t[2]) / p[0]
    T13 = (mp.mpf(1)/2 - p[1] - p[0]*t[0] - p[2]*t[2] + p[1]*t[1]) / p[0]
    return ((1 - 2*t[0])**2 + p[0]/p[1]*T12**2 + p[0]/p[2]*T13**2) / 2, sum(t)
for p, ref in [(("1/2", "1/4", "1/4"), "0.04680597938234766441775526"), (("2/5", "1/3", "4/15"), "0.06566305101224223533800595"), (("3/5", "1/5", "1/5"), "0.03234542578607244487622299"), (("1/3", "1/3", "1/3"), "0.0833333333333333333333")]:
    v, s = pbar(*[mp.mpf(sp) if "/" not in sp else mp.mpf(sp.split("/")[0]) / mp.mpf(sp.split("/")[1]) for sp in p])
    print(p, v, "engine", ref, "sum t =", s)
