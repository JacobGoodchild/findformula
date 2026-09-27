import mpmath as mp, sympy as sp, pickle, sys
sys.argv = ["x", "34"]
src = open("lg_fcc.py").read().split("deg = 4 * P - 2")[0]
exec(src)
nbrs_, classes, cache_ = pickle.load(open("lgFCC.pkl", "rb"))
cache.update(cache_)
pi = mp.pi; s3 = mp.sqrt(3)
Wsc = mp.sqrt(6)/(32*pi**3)*mp.gamma(mp.mpf(1)/24)*mp.gamma(mp.mpf(5)/24)*mp.gamma(mp.mpf(7)/24)*mp.gamma(mp.mpf(11)/24)
Yn = Wsc * mp.sqrt(2)
def hd(n, t=mp.mpf(8)):
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2)
        A = t + c1 * c2; B = c1 + c2; D = A * A - B * B
        return [1 / mp.sqrt(D), -A / D**1.5, (2 * A * A + B * B) / D**2.5][n]
    return mp.quad(f, [0, mp.pi / 2, mp.pi], [0, mp.pi / 2, mp.pi]) / mp.pi**2
hs = [hd(0), hd(1), hd(2)]
Wf_, Xf_, Y_, XY_, R3_, h0, h1, h2 = sp.symbols('W_f X_f Y X_Y r3 h0 h1 h2', positive=True)   # X_f = 1/(pi^2 W_f), X_Y = 1/(pi^2 Y), r3 = sqrt3/pi
numv = {Wf_: Wf, Xf_: 1/(pi**2*Wf), Y_: Yn, XY_: 1/(pi**2*Yn), R3_: s3/pi, h0: hs[0], h1: hs[1], h2: hs[2]}
sym = {}
for (kind, x), v in cache.items():
    if kind == "L":
        bas = [Wf_, Xf_]; vals = [4 * v]                      # 4v = W_f/3 - Delta
        target = 4 * v
    elif kind == "S":
        bas = [Y_, XY_, R3_]; target = 4 * v                  # 4v = E[cos/(3+s)]
    else:
        bas = [h0, h1, h2]; target = 4 * v                    # 4v = E[cos/(8+s)]
    if kind == "3":
        r = pickle.load(open("fcc_t8_rels.pkl", "rb"))[x]          # verified at 60 digits
    else:
        r = mp.pslq([target, 1] + [numv[b] for b in bas], maxcoeff=10**8, maxsteps=10**7)
    assert r and r[0] != 0, (kind, x, target)
    sym[(kind, x)] = -(r[1] + sum(c * b for c, b in zip(r[2:], bas))) / sp.Integer(r[0]) / 4
    print(kind, x, r, flush=True)
def Gs(kind, x): return sym[(kind, tuple(sorted(abs(t) for t in x)))]
def BdBs(kind, i, Ri, j, Rj):
    x = add(Rj, neg(Ri))
    return sum(Gs(kind, add(add(x, a_), b_)) for a_, b_ in [(Z, Z), (neg(E[i]), Z), (Z, E[j]), (neg(E[i]), E[j])])
def Ls(i, Ri, j, Rj): return ((1 if (i, Ri) == (j, Rj) else 0) + BdBs("L", i, Ri, j, Rj)) / sp.Integer(4 * P)
def Qs(i, Ri, j, Rj): return ((1 if (i, Ri) == (j, Rj) else 0) - BdBs("3", i, Ri, j, Rj)) / sp.Integer(4 * P - 4)
def Fs(i, Ri, j, Rj): return (1 if (i, Ri) == (j, Rj) else 0) - BdBs("S", i, Ri, j, Rj)
deg = 22; mu0 = -sp.Rational(1, 11); th2 = 1 - mu0**2
subs = {k: sp.Float(str(v), 40) for k, v in numv.items()}
out = {}
for key, (a, pp, pm, fl, Yr, Qr) in classes.items():
    Ye = [Ls(*o, *o) - Ls(*o, *b) - Ls(*a, *o) + Ls(*a, *b) for b in nbrs]
    Qe = [Qs(*o, *o) + Qs(*o, *b) + Qs(*a, *o) + Qs(*a, *b) for b in nbrs]
    ppe = sp.expand(sum(((1 if b == a else 0) - y)**2 for b, y in zip(nbrs, Ye)) / 4)
    pme = sp.expand(sum(((1 if b == a else 0) - y)**2 for b, y in zip(nbrs, Qe)) / 4)
    fle = 0
    for sg in (1, -1):
        lam = mu0 + sg * sp.I * sp.sqrt(th2)
        for b in nbrs:
            amp = (Fs(*o, *o) - lam * Fs(*b, *o) - sp.conjugate(lam) * Fs(*o, *a) + Fs(*b, *a)) / (2 * th2 * deg)
            fle += sp.expand(amp * sp.conjugate(amp))
    fle = sp.expand(sp.simplify(fle))
    print("start", a)
    for nm, e, ref in (("+1", ppe, pp), ("flat", fle, fl), ("-1", pme, pm)):
        print("  ", nm, "=", e); print("      =", sp.N(e.subs(subs), 30), " engine", mp.nstr(ref, 30))
    out[a] = (ppe, fle, pme)
pickle.dump((out, {str(k): str(v) for k, v in numv.items()}), open("lgfcc_sym.pkl", "wb"))
