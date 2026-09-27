import mpmath as mp, sympy as sp
from honey_w import honey
mp.mp.dps = 30
for p in [(sp.Rational(2, 5), sp.Rational(7, 20), sp.Rational(1, 4)), (sp.Rational(11, 25), sp.Rational(8, 25), sp.Rational(6, 25)), (sp.Rational(9, 20), sp.Rational(7, 20), sp.Rational(1, 5))]:
    Y = honey(p)
    pm = [mp.mpf(sp.N(x, 40)) for x in p]
    lam = mp.findroot(lambda l: sum(mp.atan(l * x) for x in pm) - mp.pi, 3)
    t = [mp.atan(lam * x) / mp.pi for x in pm]
    T12 = mp.sqrt(pm[1] / pm[0]) * Y[0][1]
    r = mp.pslq([T12, 1, t[0], t[1]], maxcoeff=10**5, maxsteps=10**6)
    print(p, r, " -> T12 = (", -r[1], "+", -r[2], "t1 +", -r[3], "t2 ) /", r[0])
