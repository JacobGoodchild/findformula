import mpmath as mp, sympy as sp, pickle, sys
sys.argv = ["x", "32"]
exec(open("lg_bcc.py").read().split("deg = 4 * P - 2")[0])   # reuse definitions (G, Linv, Fm, nbrs ...)
pi = mp.pi
Wn = mp.gamma(mp.mpf(1)/4)**4/(4*pi**3)
Xn = 1/(pi**2*Wn)
W, X = sp.symbols('W X', positive=True)
def ex(v):
    if abs(v) < mp.mpf(10)**-25: return sp.Integer(0)
    r = mp.pslq([v, 1, Wn, Xn], maxcoeff=10**6, maxsteps=10**6)
    assert r and r[0] != 0, v
    return -(r[1] + r[2]*W + r[3]*X) / sp.Integer(r[0])
deg = 14; mu0 = -sp.Rational(1, 7); th2 = 1 - mu0**2
import pickle
_, classes, _ = pickle.load(open('lgBCC.pkl','rb'))
reps = [v[0] for v in classes.values()]
for a in reps:
    Ye = [ex(Yent(Linv, a, b, False)) for b in nbrs]
    pp = sp.expand(sum(((1 if b == a else 0) - y)**2 for b, y in zip(nbrs, Ye)) / 4)
    Foo = ex(Fm(*o, *o)); Foa = ex(Fm(*o, *a))
    fl = 0
    for sg in (1, -1):
        lam = mu0 + sg*sp.I*sp.sqrt(th2)
        for b in nbrs:
            amp = (Foo - lam*ex(Fm(*b, *o)) - sp.conjugate(lam)*Foa + ex(Fm(*b, *a))) / (2*th2*deg)
            fl += sp.expand(amp*sp.conjugate(amp))
    fl = sp.expand(sp.simplify(fl))
    num = lambda e: sp.N(e.subs({W: sp.Float(str(Wn), 35), X: sp.Float(str(Xn), 35)}), 25)
    print("start", a); print("  +1 =", pp, "=", num(pp)); print("  flat =", fl, "=", num(fl))
