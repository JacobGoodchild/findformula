"""Anisotropic triangular Szegedy walk: bond directions e1, e2, e1-e2 with probabilities q1, q2, q3 per
direction (each direction has 2 arcs, P(arc) = q_i/2).  L_P(k) = q1(1-c1) + q2(1-c2) + q3(1-cos(k1-k2)),
Q_P(k) = q1(1+c1) + q2(1+c2) + q3(1+cos(k1-k2))."""
import sys, mpmath as mp
sys.path.insert(0, '../quantum2d')
from lgf2d import a as akern, b as bgreen
mp.mp.dps = 25
dirs = [(1, 0), (0, 1), (1, -1)]
def channels(q, start=0):
    q = [mp.mpf(x) for x in q]
    # L_P = A(k1) - Re(C(k1) e^{ik2}), A = q1(1-c1) + q2 + q3, C = q2 + q3 e^{-ik1}
    A = lambda k: q[0] * (1 - mp.cos(k)) + q[1] + q[2]; C = lambda k: q[1] + q[2] * mp.exp(-1j * k)
    Ap = lambda k: q[0] * (1 + mp.cos(k)) + q[1] + q[2]; Cp = lambda k: -(q[1] + q[2] * mp.exp(-1j * k))
    arcs = [(s, p) for p in range(3) for s in (1, -1)]
    P = lambda arc: q[arc[1]] / 2
    ca, cg = {}, {}
    def aa(y):
        if y == (0, 0): return mp.mpf(0)
        if y not in ca: ca[y] = akern(y, A, C)
        return ca[y]
    def gp(y):
        if y not in cg: cg[y] = bgreen(y, Ap, Cp)
        return cg[y]
    a0 = arcs[start]
    Yr, Qr = [], []
    for b in arcs:
        dp = tuple(a0[0] * v for v in dirs[a0[1]]); dq = tuple(b[0] * v for v in dirs[b[1]])
        y = (dq[0] - dp[0], dq[1] - dp[1])
        w = mp.sqrt(P(a0) * P(b))
        Yr.append(w * (aa(dp) + aa(dq) - aa(y)))                 # G0 - G(dq) - G(dp) + G(dq - dp)
        Qr.append(w * (gp((0, 0)) + gp(dp) + gp(dq) + gp(y)))
    pp = mp.fsum(((1 if b == a0 else 0) - y)**2 for b, y in zip(arcs, Yr)) / 4
    pm = mp.fsum(((1 if b == a0 else 0) - y)**2 for b, y in zip(arcs, Qr)) / 4
    return pp, pm, Yr, Qr
if __name__ == "__main__":
    import math
    for q in [(mp.mpf(1)/3,)*3, (mp.mpf(1)/2, mp.mpf(1)/3, mp.mpf(1)/6), (mp.mpf(2)/5, mp.mpf(7)/20, mp.mpf(1)/4)]:
        pp, pm, Yr, Qr = channels(q)
        lam = mp.findroot(lambda l: sum(mp.atan(l * x) for x in q) - mp.pi / 2, 1)
        th = [mp.atan(lam * x) for x in q]
        print("q", [mp.nstr(x, 5) for x in q], "+1", pp, "-1", pm, "total", pp + pm)
        print("   Y(a,a) =", Yr[0], " 2theta1/pi =", 2 * th[0] / mp.pi, " lam^2 =", mp.identify(lam**2))
