"""Abelian sandpile height-1 probability via Majumdar-Dhar: P1 = det(I - M), M_jk = D(x_j) + D(x_k) - D(x_j - x_k),
D(x) = G(0) - G(x) = R(x)/2, over the neighbours x_j whose bonds to the origin are cut (all but one)."""
import mpmath as mp, itertools, sys
def R_bravais(bonds, x, dps=30):
    """effective resistance to x on a one-site lattice with bonds [(vec, conductance)]."""
    mp.mp.dps = dps
    def f(k1, k2):
        D = mp.fsum(4 * c * mp.sin((k1 * v[0] + k2 * v[1]) / 2)**2 for v, c in bonds)
        if D == 0: return mp.mpf(0)
        return 4 * mp.sin((k1 * x[0] + k2 * x[1]) / 2)**2 / D
    return mp.quad(f, [-mp.pi, 0, mp.pi], [-mp.pi, 0, mp.pi]) / (4 * mp.pi**2)
def P1(bonds, keep=0, dps=30, Rfun=None):
    nbrs = []
    for v, c in bonds:
        nbrs += [tuple(v), tuple(-a for a in v)]
    cut = [n for i, n in enumerate(nbrs) if i != keep]
    cache = {}
    def D(x):
        key = tuple(sorted([abs(x[0]), abs(x[1])])) if Rfun is None else x
        if x == (0, 0): return mp.mpf(0)
        if key not in cache: cache[key] = (Rfun(x) if Rfun else R_bravais(bonds, x, dps)) / 2
        return cache[key]
    n = len(cut)
    M = mp.matrix(n, n)
    for j in range(n):
        for k in range(n):
            d = (cut[j][0] - cut[k][0], cut[j][1] - cut[k][1])
            M[j, k] = D(cut[j]) + D(cut[k]) - D(d)
    return mp.det(mp.eye(n) - M), cache
if __name__ == "__main__":
    sq = [((1, 0), 1), ((0, 1), 1)]
    v, c = P1(sq, dps=20); mp.mp.dps = 20
    print("square P1 =", v, " expected 2/pi^2 - 4/pi^3 =", 2 / mp.pi**2 - 4 / mp.pi**3)
