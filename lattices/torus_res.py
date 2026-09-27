"""Independent check of adjacent-node resistances: solve Kirchhoff's laws on an N x N torus of unit cells."""
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl, sys
from builders import LATTICES
def resistances(nv, E, N):
    idx = lambda v, x, y: v + nv * ((x % N) + N * (y % N))
    n = nv * N * N; rows = []; cols = []
    for x in range(N):
        for y in range(N):
            for (u, v, o) in E:
                a, b = idx(u, x, y), idx(v, x + o[0], y + o[1]); rows += [a, b]; cols += [b, a]
    A = sps.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    L = sps.diags(np.asarray(A.sum(1)).ravel()) - A
    L = L.tolil(); L[0, 0] += 1.0; L = L.tocsc()           # ground node 0 weakly (removes zero mode)
    lu = spl.splu(L)
    out = {}
    for (u, v, o) in E:
        a, b = idx(u, N // 2, N // 2), idx(v, N // 2 + o[0], N // 2 + o[1])
        rhs = np.zeros(n); rhs[a] = 1; rhs[b] = -1
        phi = lu.solve(rhs); out[(u, v, o)] = phi[a] - phi[b]
    return out
if __name__ == "__main__":
    nv, E = LATTICES[sys.argv[1]]()
    for N in map(int, sys.argv[2:]):
        r = resistances(nv, E, N)
        print(N, sorted(set(round(x, 7) for x in r.values())))
