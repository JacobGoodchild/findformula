"""Sandpile height-1 probability on the simple cubic lattice via Majumdar-Dhar with exact resistances.
R(x) = 2 INT_0^inf e^{-6t} [I0^3 - prod I_{x_i}] dt (unit resistors, L = 6 - 2 sum cos)."""
import mpmath as mp, itertools, sympy as sp
mp.mp.dps = 30
def R(x):
    f = lambda t: mp.exp(-6 * t) * (mp.besseli(0, 2 * t)**3 - mp.fprod(mp.besseli(abs(a), 2 * t) for a in x))
    return 2 * mp.quad(f, [0, 1, 5, 20, 80, 320, mp.inf])
W = mp.sqrt(6) / (32 * mp.pi**3) * mp.gamma(mp.mpf(1)/24) * mp.gamma(mp.mpf(5)/24) * mp.gamma(mp.mpf(7)/24) * mp.gamma(mp.mpf(11)/24)
vals = {}
for x in [(1,0,0), (1,1,0), (2,0,0)]:
    vals[x] = R(x); print(x, vals[x], mp.pslq([vals[x], 1, W, 1/(mp.pi**2 * W)], maxcoeff=10**5, maxsteps=10**6))
