"""Flip-flop Grover walk trapping on a 2D Bravais lattice with neighbour pairs +-delta_p.
pbar(a) = ||M+ e_a||^2 + ||M- e_a||^2, M = I/2 - N,
N+ = (1/4)[a(dp) + a(dq) - a(y)],  N- = (1/4)[G'(0)+G'(dp)+G'(dq)+G'(y)],  y = t dq - s dp."""
import mpmath as mp
from lgf2d import a as akern, b as bgreen

def trap(deltas, A, C, Ap, Cp, bipartite=False):
    P = len(deltas)
    states = [(s, p) for p in range(P) for s in (1, -1)]
    cacheA, cacheG = {}, {}
    def aa(y):
        if y == (0, 0): return mp.mpf(0)
        if y not in cacheA: cacheA[y] = akern(y, A, C)
        return cacheA[y]
    def gp(y):
        if y not in cacheG: cacheG[y] = bgreen(y, Ap, Cp)
        return cacheG[y]
    n = 2 * P
    Np, Nm = mp.matrix(n, n), mp.matrix(n, n)
    for i, (s, p) in enumerate(states):
        for j, (t, q) in enumerate(states):
            dp, dq = deltas[p], deltas[q]
            y = (t * dq[0] - s * dp[0], t * dq[1] - s * dp[1])
            Np[i, j] = (aa(dp) + aa(dq) - aa(y)) / 4
            if not bipartite:
                Nm[i, j] = (gp((0, 0)) + gp(dp) + gp(dq) + gp(y)) / 4
    Mp = mp.eye(n) / 2 - Np; Mm = mp.eye(n) / 2 - Nm
    col = lambda M: sum(M[k, 0]**2 for k in range(n))
    if bipartite:
        return 2 * col(Mp), col(Mp), col(Mp), cacheA, cacheG
    return col(Mp) + col(Mm), col(Mp), col(Mm), cacheA, cacheG

if __name__ == "__main__":
    mp.mp.dps = 30
    sq = trap([(1, 0), (0, 1)], lambda k: 2 - mp.cos(k), lambda k: mp.mpf(1), None, None, bipartite=True)
    print("square flip-flop pbar:", sq[0], sq[1], sq[2])
    tri = trap([(1, 0), (0, 1), (1, -1)], lambda k: 3 - mp.cos(k), lambda k: 1 + mp.exp(-1j * k),
               lambda k: 3 + mp.cos(k), lambda k: -(1 + mp.exp(-1j * k)))
    print("triangular flip-flop pbar:", tri[0], " (+1:", tri[1], " -1:", tri[2], ")")
    print("a values:", tri[3]); print("G' values:", tri[4])
