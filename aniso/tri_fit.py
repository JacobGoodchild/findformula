import mpmath as mp
from tri_w import channels
mp.mp.dps = 30
arcs = [(s, p) for p in range(3) for s in (1, -1)]
for q in [(mp.mpf(2)/5, mp.mpf(7)/20, mp.mpf(1)/4), (mp.mpf(9)/20, mp.mpf(7)/20, mp.mpf(1)/5), (mp.mpf(11)/25, mp.mpf(8)/25, mp.mpf(6)/25)]:
    pp, pm, Yr, Qr = channels(q)
    lam = mp.findroot(lambda l: sum(mp.atan(l * x) for x in q) - mp.pi / 2, 1)
    t = [mp.atan(lam * x) / mp.pi for x in q]
    print("q", [mp.nstr(x, 4) for x in q], "lam^2", mp.identify(lam**2), [mp.nstr(x, 6) for x in [1/(q[0]*q[1]+q[1]*q[2]+q[2]*q[0])]])
    for b, y in zip(arcs, Yr):
        T = mp.sqrt(q[b[1]] / q[0]) * y
        c1 = 1 / (lam * q[0]); r = mp.pslq([T, 1, t[0], t[1], c1 / mp.pi], maxcoeff=10**5, maxsteps=10**6)
        print("   ", b, mp.nstr(T, 18), r)
