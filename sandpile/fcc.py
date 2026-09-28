"""FCC sandpile P1. Neighbours: permutations of (+-1,+-1,0). L(k) = 12 - 4s, s = c1c2 + c2c3 + c3c1."""
import mpmath as mp, sympy as sp, itertools
mp.mp.dps = 25
def R(x):
    m1, m2, m3 = map(abs, x)
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2); A = 3 - c1 * c2; B = c1 + c2
        r = mp.sqrt(A * A - B * B)
        if r == 0: return mp.mpf(0)
        t = mp.mpf(1) if m3 == 0 else (((A - r) / B)**m3 if B != 0 else mp.mpf(0))
        return (1 - mp.cos(m1 * k1) * mp.cos(m2 * k2) * t) / r
    return mp.quad(f, [0, mp.pi / 2, mp.pi], [0, mp.pi / 2, mp.pi]) / mp.pi**2 / 2
Wf = 3 * mp.gamma(mp.mpf(1) / 3)**6 / (2**(mp.mpf(14) / 3) * mp.pi**4)
vals = {}
for x in [(1,1,0), (2,0,0), (2,1,1), (2,2,0)]:
    vals[x] = R(x); print(x, vals[x], mp.pslq([vals[x], 1, Wf, 1 / (mp.pi**2 * Wf)], maxcoeff=10**5, maxsteps=10**6))
