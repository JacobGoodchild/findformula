import mpmath as mp
from tri_w import channels
mp.mp.dps = 25
def law(q):
    q1, q2, q3 = [mp.mpf(x) for x in q]
    lam = 1 / mp.sqrt(q1*q2 + q2*q3 + q3*q1)
    th = [mp.atan(lam * x) for x in (q1, q2, q3)]
    t1, t2 = th[0] / mp.pi, th[1] / mp.pi
    c1 = 1 / (lam * q1)
    T = {}
    T[(1, 0)] = 2 * t1
    T[(-1, 0)] = 2 / mp.pi * (c1 - th[0] * c1**2)
    T[(1, 1)] = -q2 / (2 * q3) + q2 * (q1 + q3) / (q1 * q3) * t1 + (q2 + q3) / q3 * t2
    T[(-1, 1)] = q2 * (q1 + q3) / q1**2 * t1 + (q2 + q3) / q2 * t2 - c1 / mp.pi
    T[(1, 2)] = mp.mpf(1) / 2 - (q1 - q3) / q1 * t1 - (q2 + q3) / q2 * t2
    T[(-1, 2)] = 1 - sum(T.values())
    return T, lam, th
for q in [(mp.mpf('0.37'), mp.mpf('0.41'), mp.mpf('0.22')), (mp.mpf('0.2'), mp.mpf('0.3'), mp.mpf('0.5')), (mp.mpf(1)/3,)*3]:
    T, lam, th = law(q)
    pp, pm, Yr, Qr = channels(q)
    arcs = [(s, p) for p in range(3) for s in (1, -1)]
    err = max(abs(T[b] - mp.sqrt(q[b[1]] / q[0]) * y) for b, y in zip(arcs, Yr))
    print("q", [mp.nstr(x, 3) for x in q], " sum theta = pi/2?", mp.nstr(sum(th) - mp.pi / 2, 5), " max |law - engine| =", mp.nstr(err, 5))
