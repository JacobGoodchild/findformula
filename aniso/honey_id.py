import mpmath as mp, sympy as sp
from honey_w import honey
mp.mp.dps = 25
def thetas(p):
    lam = mp.findroot(lambda l: sum(mp.atan(l * x) for x in p) - mp.pi, 3)
    return [mp.atan(lam * x) for x in p], lam
for p in [(sp.Rational(2, 5), sp.Rational(1, 3), sp.Rational(4, 15)), (sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 6))]:
    Y = honey(p)
    pm = [mp.mpf(sp.N(x, 30)) for x in p]
    th, lam = thetas(pm)
    print("p", p, "lambda", lam, "2theta/pi", [mp.nstr(2 * t / mp.pi, 15) for t in th], "Ydiag", [mp.nstr(Y[i][i], 15) for i in range(3)])
    for j in (1, 2):
        T = mp.sqrt(pm[j] / pm[0]) * Y[0][j]          # physical current through bond j
        B = [1, th[0] / mp.pi, th[1] / mp.pi, mp.cot(th[0]) / mp.pi, th[0] * mp.cot(th[0]) / mp.pi, th[1] * mp.cot(th[0]) / mp.pi]
        print("   T(1,%d) =" % (j + 1), mp.nstr(T, 20), mp.pslq([T] + B, maxcoeff=1000, maxsteps=10**6))
