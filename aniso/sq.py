"""Szegedy walk (flip-flop, weighted Grover coin) on the anisotropic square lattice:
P(+-x) = alpha/2, P(+-y) = beta/2, alpha + beta = 1.  L_P(k) = alpha(1-c1) + beta(1-c2).
a(x) = E[(1 - cos k.x)/L_P(k)]; k2 integrated exactly."""
import mpmath as mp
mp.mp.dps = 30
def a(x, al):
    be = 1 - al; m1, m2 = x
    def f(k1):
        A = al * (1 - mp.cos(k1)) + be          # L = A - be*cos k2
        root = mp.sqrt(A * A - be * be)
        if root == 0: return mp.mpf(0)
        r = (A - root) / be
        return (1 - mp.cos(m1 * k1) * r**abs(m2)) / root
    return mp.quad(f, [0, mp.pi]) / mp.pi
def trap(al):
    be = 1 - al
    a10, a01, a20, a11 = a((1, 0), al), a((0, 1), al), a((2, 0), al), a((1, 1), al)
    Yxx = al * a10; Yxmx = al / 2 * (2 * a10 - a20); Yxy = mp.sqrt(al * be) / 2 * (a10 + a01 - a11)
    rowsum = Yxx + Yxmx + 2 * Yxy
    pbar = 2 * ((1 - Yxx)**2 + Yxmx**2 + 2 * Yxy**2) / 4
    return pbar, (Yxx, Yxmx, Yxy, rowsum)
for al in [mp.mpf(1)/2, mp.mpf(1)/3, mp.mpf(1)/4, mp.mpf(1)/5, mp.mpf(2)/3]:
    p, Y = trap(al)
    th = mp.asin(mp.sqrt(al))
    print("alpha", al, "pbar", p, "rowsum", Y[3])
    for nm, v in zip(["Yxx", "Yx-x", "Yxy"], Y[:3]):
        B = [1, th / mp.pi, mp.sqrt(al * (1 - al)) / mp.pi, mp.sqrt(al / (1 - al)) / mp.pi, mp.sqrt((1 - al) / al) / mp.pi, th * mp.sqrt(al * (1 - al)) / mp.pi]
        print("   ", nm, v, mp.pslq([v] + B, maxcoeff=10**4, maxsteps=10**6))
