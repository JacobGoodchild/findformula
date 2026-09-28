import mpmath as mp
def closed(eps):
    eps = mp.mpf(eps); P = (1 - eps) / 3; al = (1 + eps) / 2; be = P / 2
    r = 3 * (1 + eps) / (1 - eps); tau = (r * r - 3) / 2; e = r
    m = 16 * e / ((e - 1)**3 * (e + 3)); g = 4 * mp.ellipk(m) / (mp.pi * mp.sqrt((e - 1)**3 * (e + 3)))
    G00 = al / (2 * be * be) * g
    Gn = -((3 + 2 * tau) * g - 2) / (6 * be)
    Gt = al / (2 * be * be) * (tau * g - 1) / 3
    As = (1 - P * (G00 + Gn)) / 2
    Ao = -(P / 4) * (G00 + 2 * Gn + Gt)
    Al = -(mp.sqrt(eps) * mp.sqrt(P) / 2) * (G00 + Gn)
    return mp.mpf(1) / 24 + As**2 + 2 * Ao**2 + Al**2, (As, Ao, Al, g, m)
if __name__ == "__main__":
    mp.mp.dps = 20
    for e in ['1e-12', '0.3', '0.6']: print(e, closed(e)[0])
