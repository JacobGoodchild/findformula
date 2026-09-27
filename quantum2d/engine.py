"""Flip-flop Grover walk trapping on 2D periodic graphs whose Bloch determinant reduces to a
triangular-lattice dispersion  kappa * (alpha - c1 - c2 - cos(k1-k2)).
Uses:  <b|Pi_+|a> = (delta - Y(a,b))/2,  Y = d^* L^{-1} d (transfer currents),
       <b|Pi_-|a> = (delta - YQ(a,b))/2, YQ = |d|^* Q^{-1} |d| (signless Laplacian)."""
import sympy as sp, mpmath as mp
from lgf2d import a as akern, b as bgreen

z1, z2 = sp.symbols('z1 z2')
def bloch(nv, edges, signless):
    M = sp.zeros(nv, nv)
    for (u, v, (o1, o2)) in edges:
        ph = z1**o1 * z2**o2          # v sits in cell R + o
        M[u, u] += 1; M[v, v] += 1
        s = 1 if signless else -1
        M[u, v] += s * ph; M[v, u] += s / ph
    return M

def laurent_terms(expr):
    expr = sp.expand(expr)
    poly = sp.Poly(sp.expand(expr * z1**20 * z2**20), z1, z2)
    return {(m[0] - 20, m[1] - 20): c for m, c in poly.terms()}

def match_triangular(det):
    """det = kappa*(alpha - (z1+1/z1+z2+1/z2+z1/z2+z2/z1)/2)"""
    T = laurent_terms(det)
    kappa = -2 * T.get((1, 0), 0)
    alpha = T.get((0, 0), 0) / kappa
    ref = laurent_terms(kappa * (alpha - (z1 + 1/z1 + z2 + 1/z2 + z1/z2 + z2/z1) / 2))
    assert all(sp.simplify(T.get(k, 0) - ref.get(k, 0)) == 0 for k in set(T) | set(ref)), (T, ref)
    return sp.nsimplify(kappa), sp.nsimplify(alpha)

def expect(numer, kappa, alpha, cache):
    """E[numer/det] with det = kappa*(alpha - ...)."""
    T = laurent_terms(numer)
    sgn = 1
    if alpha < 0:       # det = kappa*(alpha - S) = (-kappa)*((-alpha) + S): flip
        A = lambda k: -alpha + mp.cos(k); C = lambda k: -(1 + mp.exp(-1j * k)); sgn = -1
    else:
        A = lambda k: alpha - mp.cos(k); C = lambda k: 1 + mp.exp(-1j * k)
    tot = mp.mpf(0)
    if alpha == 3:     # singular: use sum c_x (G(x)-G(0)) = -sum c_x a(x); needs sum c_x = 0
        assert sp.simplify(sum(T.values())) == 0
        for x, c in T.items():
            if x == (0, 0): continue
            key = ("a", x)
            if key not in cache: cache[key] = akern(x, A, C)
            tot += -mp.mpf(sp.N(c, 50)) * cache[key]
    else:
        for x, c in T.items():
            key = ("g", float(alpha), x)
            if key not in cache: cache[key] = bgreen(x, A, C)
            tot += mp.mpf(sp.N(c, 50)) * cache[key]
    return sgn * tot / mp.mpf(sp.N(kappa, 50))

def transfer(nv, edges, origin, signless, cache):
    L = bloch(nv, edges, signless)
    det = sp.factor(L.det())
    kappa, alpha = match_triangular(sp.expand(L.det()))
    adj = L.adjugate()
    # arcs out of origin vertex (in cell 0): list of (neighbour vertex, offset)
    arcs = []
    for (u, v, o) in edges:
        if u == origin: arcs.append((v, o))
        if v == origin: arcs.append((u, (-o[0], -o[1])))
    n = len(arcs)
    Y = [[None] * n for _ in range(n)]
    for i, (vi, oi) in enumerate(arcs):
        for j, (vj, oj) in enumerate(arcs):
            # Y = [G(o,o) -/+ G(o,vj) -/+ G(vi,o) + G(vi,vj)] with Bloch phases
            s = 1 if signless else -1
            g = lambda p, q, op, oq: adj[p, q] * z1**(op[0] - oq[0]) * z2**(op[1] - oq[1])
            num = g(origin, origin, (0,0), (0,0)) + s * g(origin, vj, (0,0), oj) + s * g(vi, origin, oi, (0,0)) + g(vi, vj, oi, oj)
            Y[i][j] = expect(num, kappa, alpha, cache)
    return Y, arcs, (kappa, alpha)

def trap(nv, edges, origin, bipartite, cache=None):
    cache = {} if cache is None else cache
    Y, arcs, ka = transfer(nv, edges, origin, False, cache)
    n = len(arcs)
    pp = sum(((1 if j == 0 else 0) - Y[0][j])**2 for j in range(n)) / 4
    if bipartite:
        return 2 * pp, pp, pp, Y, None, ka
    YQ, _, kaq = transfer(nv, edges, origin, True, cache)
    pm = sum(((1 if j == 0 else 0) - YQ[0][j])**2 for j in range(n)) / 4
    return pp + pm, pp, pm, Y, YQ, (ka, kaq)

if __name__ == "__main__":
    mp.mp.dps = 30
    tri = [(0, 0, (1, 0)), (0, 0, (0, 1)), (0, 0, (1, -1))]
    r = trap(1, tri, 0, False)
    print("triangular check:", r[0], "(expect 0.2679557025370805512940)")
