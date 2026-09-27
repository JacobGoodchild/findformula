import mpmath as mp, sympy as sp, pickle
from ident import lin
mp.mp.dps = 25; pi = mp.pi; s15 = mp.sqrt(15); s3 = mp.sqrt(3)
arcs, Y, YQ, res = pickle.load(open('lgtri_first.pkl', 'rb'))
ms = mp.mpf(1)/2 - 23*s15/180; Kn = mp.ellipk(ms)/pi; En = mp.ellipe(ms)/pi; rn = mp.root(135, 4)
Xn = mp.cbrt(4)*mp.gamma(mp.mpf(1)/3)**3/(8*pi**2)
K, E, P, X = sp.symbols('K E pi X')     # K = K(m*)/pi, E = E(m*)/pi
r = sp.Integer(135)**sp.Rational(1, 4)
S = {'1': 1, 'K/r': K/r, 's15K/r': sp.sqrt(15)*K/r, 'E/r': E/r, 's15E/r': sp.sqrt(15)*E/r,
     's3/pi': sp.sqrt(3)/P, 'X': X, 's3/(piX)': sp.sqrt(3)/(P*X)}
num = {'1': 1, 'K/r': Kn/rn, 's15K/r': s15*Kn/rn, 'E/r': En/rn, 's15E/r': s15*En/rn, 's3/pi': s3/pi, 'X': Xn, 's3/(piX)': s3/(pi*Xn)}
def exact(v, names):
    rel = lin(v, {n: num[n] for n in names}, maxc=20000)
    d = dict(rel); c = d.pop('v')
    return -sum(sp.Integer(k) * S[n] for n, k in d.items()) / c
BQ = ['1', 'K/r', 's15K/r', 'E/r', 's15E/r']; BL = ['1', 's3/pi']; BF = ['1', 'X', 's3/(piX)']
out = {}
for s in (0,):
    pp, pm, fl, tot, Fba, Foo, Foa = res[s]
    Ye = [exact(y, BL) for y in Y[s]]; Qe = [exact(y, BQ) for y in YQ[s]]
    Fe = [exact(f, BF) for f in Fba]; Fooe = exact(Foo, BF); Foae = exact(Foa, BF)
    Fbo = [exact(f, BF) for f in pickle.load(open('lgtri_Fbo.pkl', 'rb'))]
    n = len(arcs)
    ppe = sum(((1 if j == s else 0) - Ye[j])**2 for j in range(n)) / 4
    pme = sum(((1 if j == s else 0) - Qe[j])**2 for j in range(n)) / 4
    mu0 = -sp.Rational(1, 5); th2 = 1 - mu0**2
    fle = 0
    for sg in (1, -1):
        lam = mu0 + sg*sp.I*sp.sqrt(th2)
        for j in range(n):
            a = (Fooe - lam*Fbo[j] - sp.conjugate(lam)*Foae + Fe[j]) / (2*th2*10)
            fle += sp.expand(a*sp.conjugate(a))
    tote = sp.simplify(sp.expand(ppe + pme + fle))
    subs = {K: Kn, E: En, P: pi, X: Xn}
    f = lambda e: mp.mpf(str(sp.N(e.subs({K: sp.Float(str(Kn), 30), E: sp.Float(str(En), 30), P: sp.pi, X: sp.Float(str(Xn), 30)}), 25)))
    print("+1 exact:", sp.simplify(ppe), " =", f(ppe), " engine", pp)
    print("-1 exact:", sp.collect(sp.expand(pme), [K, E]), " =", f(pme), " engine", pm)
    print("flat exact:", sp.simplify(fle), "=", f(sp.expand(fle)), " engine", fl)
    print("TOTAL =", f(tote), " engine", tot)
    out[s] = (ppe, pme, sp.expand(fle), tote)
pickle.dump(out, open('lgtri_exact.pkl', 'wb'))
