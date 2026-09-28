"""Identify honeycomb+NNN resistances in {1, sqrt3/pi, g, g'} (g = triangular Green function at t = 13/2)."""
import pickle, mpmath as mp, sympy as sp, sys
sys.path.insert(0, '/home/user/findformula/ust')
mp.mp.dps = 25
def gfun(t):
    return mp.quad(lambda k: 1 / mp.sqrt((t - mp.cos(k))**2 - 4 * mp.cos(k / 2)**2), [0, mp.pi]) / mp.pi
with mp.workdps(60):
    g0 = gfun(mp.mpf(13) / 2); g1 = mp.diff(gfun, mp.mpf(13) / 2)
g0, g1 = +g0, +g1
print("g check vs 4K(64/189)/(3 sqrt21 pi):", g0 - 4 * mp.ellipk(mp.mpf(64) / 189) / (3 * mp.sqrt(21) * mp.pi))
nbrs, R = pickle.load(open('R_honeycomb_nnn_0.pkl', 'rb'))
B = [1, mp.sqrt(3) / mp.pi, g0, g1]
names = [sp.Integer(1), sp.Symbol('r'), sp.Symbol('g'), sp.Symbol('gp')]
RS = {}
for (i, j), v in R.items():
    if i >= j: continue
    rel = mp.pslq([v] + B, maxcoeff=10**6, maxsteps=10**6)
    print(i, j, mp.nstr(v, 20), rel)
    if not rel or rel[0] == 0: sys.exit("no relation")
    RS[(i, j)] = RS[(j, i)] = sum(sp.Rational(-c, rel[0]) * s for c, s in zip(rel[1:], names))
for i in range(len(nbrs)): RS[(i, i)] = 0
z = len(nbrs) - 1
Y = sp.Matrix(z, z, lambda a, b: (RS[(0, a + 1)] + RS[(0, b + 1)] - RS[(a + 1, b + 1)]) / 2)
P1 = sp.factor((sp.eye(z - 1) - Y[1:, 1:]).det())
sub = {sp.Symbol('r'): sp.Float(str(B[1]), 25), sp.Symbol('g'): sp.Float(str(g0), 25), sp.Symbol('gp'): sp.Float(str(g1), 25)}
print("P1 =", P1, "=", sp.N(P1.subs(sub), 20))
from ustdeg import degree_distribution
P = degree_distribution(Y)
for k, p in enumerate(P):
    if p != 0: print("deg", k, sp.N(p.subs(sub), 18))
