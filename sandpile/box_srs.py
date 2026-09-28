"""Exact finite-box sandpile P1 on the Laves graph (srs / K4 crystal), dissipative boundary. Expect 1/12."""
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl, itertools, sys
E = [(0, 1, (0, 0, 0)), (0, 2, (0, 0, 0)), (0, 3, (0, 0, 0)), (1, 2, (1, 0, 0)), (2, 3, (0, 1, 0)), (3, 1, (0, 0, 1))]
def P1(N):
    idx = lambda s, x, y, w: s + 4 * (x + N * (y + N * w)); n = 4 * N**3; R = []; C = []
    for x, y, w in itertools.product(range(N), repeat=3):
        for (u, v, o) in E:
            X, Y, Z = x + o[0], y + o[1], w + o[2]
            if 0 <= X < N and 0 <= Y < N and 0 <= Z < N:
                a, b = idx(u, x, y, w), idx(v, X, Y, Z); R += [a, b]; C += [b, a]
    A = sps.csr_matrix((np.ones(len(R)), (R, C)), shape=(n, n))
    lu = spl.splu((3 * sps.identity(n) - A).tocsc())
    c = N // 2; o = idx(0, c, c, c); nb = list(A[o].indices); cut = nb[1:]
    B = np.zeros((n, len(cut)))
    for j, v in enumerate(cut): B[o, j] = 1; B[v, j] = -1
    GB = np.column_stack([lu.solve(B[:, j]) for j in range(len(cut))])
    return np.linalg.det(np.eye(len(cut)) - B.T @ GB)
for N in map(int, sys.argv[1:]): print(N, P1(N), " (1/12 =", 1 / 12, ")")
