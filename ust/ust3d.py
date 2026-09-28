import sympy as sp, itertools, mpmath as mp
from ustdeg import degree_distribution
def Ymat(nb, R):
    def Rv(x):
        k = tuple(sorted(map(abs, x), reverse=True))
        return 0 if not any(k) else R[k]
    return sp.Matrix(len(nb), len(nb), lambda i, j: (Rv(nb[i]) + Rv(nb[j]) - Rv(tuple(a - b for a, b in zip(nb[i], nb[j])))) / 2)
al, W, q = sp.symbols('alpha W q')
mp.mp.dps = 20
Wsc = mp.sqrt(6)/(32*mp.pi**3)*mp.gamma(mp.mpf(1)/24)*mp.gamma(mp.mpf(5)/24)*mp.gamma(mp.mpf(7)/24)*mp.gamma(mp.mpf(11)/24)
Wb = mp.gamma(mp.mpf(1)/4)**4/(4*mp.pi**3); Wf = 3*mp.gamma(mp.mpf(1)/3)**6/(2**(mp.mpf(14)/3)*mp.pi**4)
cases = {
 "SC": ([(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)],
        {(1,0,0): sp.Rational(1,3), (1,1,0): al, (2,0,0): 2 - 4*al}, {al: (7*Wsc + 54/(mp.pi**2*Wsc))/36}),
 "BCC": (list(itertools.product([1,-1], repeat=3)),
        {(1,1,1): sp.Rational(1,4), (2,0,0): (W - 4*q)/4, (2,2,0): 4*q, (2,2,2): (8 - 3*W - 36*q)/4}, {W: Wb, q: 1/(mp.pi**2*Wb)}),
}
fnb = []
for pp in set(itertools.permutations([1,1,0])):
    for s in itertools.product([1,-1], repeat=3):
        v = tuple(a*b for a, b in zip(pp, s))
        if v not in fnb: fnb.append(v)
cases["FCC"] = (fnb, {(1,1,0): sp.Rational(1,6), (2,0,0): (4*W - 3*q)/6, (2,1,1): (2*W + 3*q - 1)/3, (2,2,0): (8 - 12*W - 9*q)/3}, {W: Wf, q: 1/(mp.pi**2*Wf)})
import sys
for name in sys.argv[1:]:
    nb, R, sub = cases[name]
    Y = Ymat(nb, R)
    p = degree_distribution(Y)
    subs = {k: sp.Float(str(v), 25) for k, v in sub.items()}
    print(name, "sum =", sp.N(sum(p).subs(subs), 15))
    for k, v in enumerate(p):
        if v != 0: print("  deg", k, ":", v, " =", sp.N(v.subs(subs), 18))
