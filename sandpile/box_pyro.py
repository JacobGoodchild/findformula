"""Exact finite-box P1 for the pyrochlore lattice (4-site FCC cell; up tetrahedron {0,1,2,3}@R,
down tetrahedron {0@R, 1@R-e1, 2@R-e2, 3@R-e3}), dissipative boundary."""
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl, itertools, sys
E = [(0, 0, 0)]
def P1(N):
    idx = lambda s, x, y, w: s + 4 * (x + N * (y + N * w)); n = 4 * N**3; R = []; C = []
    e = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    for x, y, w in itertools.product(range(N), repeat=3):
        up = [(s, (x, y, w)) for s in range(4)]
        dn = [(0, (x, y, w))] + [(s, (x - e[s][0], y - e[s][1], w - e[s][2])) for s in (1, 2, 3)]
        for T in (up, dn):
            for (s1, p1), (s2, p2) in itertools.combinations(T, 2):
                if all(0 <= c < N for c in p1 + p2):
                    a, b = idx(s1, *p1), idx(s2, *p2); R += [a, b]; C += [b, a]
    A = sps.csr_matrix((np.ones(len(R)), (R, C)), shape=(n, n))
    D = (6 * sps.identity(n) - A).tocsc(); lu = spl.splu(D)
    c = N // 2; o = idx(0, c, c, c)
    nb = [idx(s, c, c, c) for s in (1, 2, 3)] + [idx(s, c - e[s][0], c - e[s][1], c - e[s][2]) for s in (1, 2, 3)]
    cut = nb[1:]; B = np.zeros((n, len(cut)))
    for j, v in enumerate(cut): B[o, j] = 1; B[v, j] = -1
    GB = np.column_stack([lu.solve(B[:, j]) for j in range(len(cut))])
    return np.linalg.det(np.eye(len(cut)) - B.T @ GB)
for N in map(int, sys.argv[1:]): print(N, P1(N))
