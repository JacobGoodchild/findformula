"""BCC sandpile P1. Sites: integer points with all coordinates same parity; neighbours (+-1,+-1,+-1).
L(k) = 8(1 - c1 c2 c3), c_i = cos k_i. R(x) = E[2(1 - cos k.x)/L]."""
import mpmath as mp, sympy as sp
mp.mp.dps = 25
def G(x):   # E[(1 - cos k.x)/(1 - c1c2c3)] with k3 integrated exactly
    m1, m2, m3 = map(abs, x)
    def f(k1, k2):
        a = mp.cos(k1) * mp.cos(k2); r = mp.sqrt(1 - a * a)
        if r == 0: return mp.mpf(0)
        t = ((1 - r) / a)**m3 if (m3 and a != 0) else (mp.mpf(1) if m3 == 0 else mp.mpf(0))
        return (1 - mp.cos(m1 * k1) * mp.cos(m2 * k2) * t) / r
    return mp.quad(f, [0, mp.pi / 2, mp.pi], [0, mp.pi / 2, mp.pi]) / mp.pi**2
Wb = mp.gamma(mp.mpf(1) / 4)**4 / (4 * mp.pi**3)
R = {}
for x in [(1,1,1), (2,0,0), (2,2,0), (2,2,2)]:
    R[x] = G(x) / 4       # 2(1-cos)/(8(1-ccc))
    print(x, R[x], mp.pslq([R[x], 1, Wb, 1 / (mp.pi**2 * Wb)], maxcoeff=10**5, maxsteps=10**6))
import pickle; pickle.dump({k: str(v) for k, v in R.items()}, open('bcc_R.pkl', 'wb'))
