"""Exact finite-box sandpile P1 (dissipative boundary) at a central site of any 2D builder lattice."""
import sys, numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl
sys.path.insert(0, '/home/user/findformula/lattices')
from builders import LATTICES
def P1(name, origin, N):
    nv, E = LATTICES[name]()
    deg = [0] * nv
    for (u, v, o) in E: deg[u] += 1; deg[v] += 1
    idx = lambda s, x, y: s + nv * (x + N * y); n = nv * N * N; R = []; C = []
    for x in range(N):
        for y in range(N):
            for (u, v, o) in E:
                X, Y = x + o[0], y + o[1]
                if 0 <= X < N and 0 <= Y < N:
                    a, b = idx(u, x, y), idx(v, X, Y); R += [a, b]; C += [b, a]
    A = sps.csr_matrix((np.ones(len(R)), (R, C)), shape=(n, n))
    th = np.array([deg[i % nv] for i in range(n)], float)
    lu = spl.splu((sps.diags(th) - A).tocsc())
    c = N // 2; o0 = idx(origin, c, c)
    nb = [j for j in A[o0].indices]
    cut = nb[1:]; B = np.zeros((n, len(cut)))
    for j, v in enumerate(cut): B[o0, j] = 1; B[v, j] = -1
    GB = np.column_stack([lu.solve(B[:, j]) for j in range(len(cut))])
    return np.linalg.det(np.eye(len(cut)) - B.T @ GB)
if __name__ == "__main__":
    name, origin = sys.argv[1], int(sys.argv[2])
    for N in map(int, sys.argv[3:]): print(name, origin, N, P1(name, origin, N))
