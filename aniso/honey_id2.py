import mpmath as mp, sympy as sp, itertools
from honey_w import honey
mp.mp.dps = 30
p = (sp.Rational(3, 7), sp.Rational(2, 7), sp.Rational(2, 7) + 0)   # placeholder, replaced below
for p in [(sp.Rational(9, 20), sp.Rational(7, 20), sp.Rational(1, 5)), (sp.Rational(1, 2), sp.Rational(3, 10), sp.Rational(1, 5))]:
    Y = honey(p)
    pm = [mp.mpf(sp.N(x, 40)) for x in p]
    lam = mp.findroot(lambda l: sum(mp.atan(l * x) for x in pm) - mp.pi, 3)
    th = [mp.atan(lam * x) for x in pm]; ct = [1 / (lam * x) for x in pm]
    print("lambda", lam, mp.identify(lam**2))
    for j in (1, 2):
        T = mp.sqrt(pm[j] / pm[0]) * Y[0][j]
        names = ["1", "t1", "t2", "c1", "c1t1", "c1t2", "c1^2", "c1^2t1", "c1^2t2"]
        vals = [1, th[0]/mp.pi, th[1]/mp.pi, ct[0]/mp.pi, ct[0]*th[0]/mp.pi, ct[0]*th[1]/mp.pi, ct[0]**2/mp.pi, ct[0]**2*th[0]/mp.pi, ct[0]**2*th[1]/mp.pi]
        r = mp.pslq([T] + vals, maxcoeff=10**4, maxsteps=10**7)
        print(" T(1,%d)" % (j + 1), mp.nstr(T, 25), [(n, c) for n, c in zip(["T"] + names, r) if c] if r else None)
