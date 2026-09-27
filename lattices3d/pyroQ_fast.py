"""Pyrochlore signless transfer currents at high precision via FCC Green functions at t = 15.
Map primitive FCC phases z1 = w2 w3, z2 = w1 w3, z3 = w1 w2 (w_i = e^{i u_i}, u uniform on the cube)."""
import sympy as sp, mpmath as mp, pickle
from gen3d import geometric3, bloch, out_arcs, ph, z1, z2, z3
mp.mp.dps = 45
w1, w2, w3 = sp.symbols('w1 w2 w3')
fcc = [(0, 2, 2), (2, 0, 2), (2, 2, 0)]
nv, E = geometric3(fcc, [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)], 2)
A, D, deg = bloch(nv, E)
Q = D + A
det = sp.expand(Q.det()); adj = Q.adjugate()
sub = {z1: w2 * w3, z2: w1 * w3, z3: w1 * w2}
detw = sp.expand(det.subs(sub, simultaneous=True))
print("det(Q) in cube coords:", sp.factor(detw))
t = 15
cache = {}
def Gx(x):
    """E_u[cos(u.x)/(t - s)], s = c1c2 + c2c3 + c3c1."""
    x = tuple(abs(v) for v in x)
    key = tuple(sorted(x[:2])) + (x[2],)
    if key in cache: return cache[key]
    m1, m2, m3 = x
    def f(u1, u2):
        c1, c2 = mp.cos(u1), mp.cos(u2)
        Aa = t - c1 * c2; Bb = c1 + c2; root = mp.sqrt(Aa * Aa - Bb * Bb)
        e3 = (((Aa - root) / Bb)**m3 / root if Bb != 0 else mp.mpf(0)) if m3 else 1 / root
        return mp.cos(m1 * u1) * mp.cos(m2 * u2) * e3
    v = mp.quad(f, [0, mp.pi], [0, mp.pi]) / mp.pi**2
    cache[key] = v; return v
def expect(N):
    """E[N/det]; det = kappa (t - s) in cube coords."""
    Nw = sp.expand(N.subs(sub, simultaneous=True))
    return Nw
# find kappa: det_w = kappa*(t - s) as Laurent polynomial -> compare constant term
const = sp.Poly(sp.expand(detw * (w1 * w2 * w3)**4), w1, w2, w3).as_dict().get((4, 4, 4), 0)
kappa = sp.Rational(const, t); print("kappa =", kappa)
def avg(N):
    Nw = sp.expand(N.subs(sub, simultaneous=True))
    P_ = sp.Poly(sp.expand(Nw * (w1 * w2 * w3)**8), w1, w2, w3)
    tot = mp.mpf(0)
    for mon, c in P_.as_dict().items():
        x = tuple(e - 8 for e in mon)
        tot += mp.mpf(sp.N(c, 60)) * Gx(x)
    return tot / mp.mpf(sp.N(kappa, 60))
arcs = out_arcs(E, 0); O = (0, (0, 0, 0))
G = lambda p, q: adj[p[0], q[0]] * ph(p[1], q[1])
ha = arcs[0]
row = []
for hb in arcs:
    v = avg(G(O, O) + G(O, hb) + G(ha, O) + G(ha, hb)); row.append(v)
    print(hb, mp.nstr(v, 40), flush=True)
t_ = mp.mpf(15)
def gd(n):
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2)
        Aa = t_ - c1 * c2; Bb = c1 + c2; Dd = Aa * Aa - Bb * Bb
        return [1 / mp.sqrt(Dd), -Aa / Dd**1.5, (2 * Aa * Aa + Bb * Bb) / Dd**2.5][n]
    return mp.quad(f, [0, mp.pi], [0, mp.pi]) / mp.pi**2
gs = [gd(0), gd(1), gd(2)]
for hb, v in zip(arcs, row):
    print(hb, mp.pslq([v, 1, gs[0], gs[1], gs[2]], maxcoeff=10**8, maxsteps=10**7))
pickle.dump((arcs, row, gs), open("pyroQ45.pkl", "wb"))
