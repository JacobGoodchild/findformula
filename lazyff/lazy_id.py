import mpmath as mp, sympy as sp
from lazy_sq import g
mp.mp.dps = 30
T = sp.symbols('t')
rows = {}
for tn in [sp.Rational(5, 2), 3, sp.Rational(7, 2), 4, 5, 6]:
    t = mp.mpf(sp.N(tn, 40)); k = 2 / t
    K = mp.ellipk(k**2) / mp.pi; E = mp.ellipe(k**2) / mp.pi
    print("t =", tn, " g0 check:", mp.nstr(g((0, 0), t) - 2 * K / t, 5))
    for x in [(1, 1), (2, 0)]:
        r = mp.pslq([g(x, t), 1, K, E], maxcoeff=10**6, maxsteps=10**6)
        rows.setdefault(x, []).append((tn, r))
        print("   ", x, r)
