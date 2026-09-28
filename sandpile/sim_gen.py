"""Generic abelian sandpile simulation on an N x N patch of a periodic lattice (open boundary: grains that
would leave the patch are lost). Heights 1..deg; topple when h > deg. Reports bulk density of h = 1 per sublattice."""
import sys, numpy as np, scipy.sparse as sps
sys.path.insert(0, '/home/user/findformula/lattices')
from builders import LATTICES
def build(name, N):
    nv, E = LATTICES[name]()
    idx = lambda s, x, y: s + nv * (x + N * y)
    deg = np.zeros(nv, int)
    for (u, v, o) in E: deg[u] += 1; deg[v] += 1
    R = []; C = []
    for x in range(N):
        for y in range(N):
            for (u, v, o) in E:
                X, Y = x + o[0], y + o[1]
                if 0 <= X < N and 0 <= Y < N:
                    a, b = idx(u, x, y), idx(v, X, Y); R += [a, b]; C += [b, a]
    n = nv * N * N
    A = sps.csr_matrix((np.ones(len(R), int), (R, C)), shape=(n, n))
    th = np.array([deg[s] for y in range(N) for x in range(N) for s in range(nv)])
    sub = np.array([s for y in range(N) for x in range(N) for s in range(nv)])
    xy = np.array([(x, y) for y in range(N) for x in range(N) for s in range(nv)])
    return A, th, sub, xy, nv
def run(name, N=80, steps=150000, seed=1):
    rng = np.random.default_rng(seed)
    A, th, sub, xy, nv = build(name, N)
    h = rng.integers(1, th + 1)
    def relax(h):
        while True:
            t = (h > th).astype(np.int64)
            if not t.any(): return h
            h = h - th * t + A @ t
    h = relax(h + th)
    for _ in range(steps // 5):
        i = rng.integers(0, len(h)); h[i] += 1
        if h[i] > th[i]: h = relax(h)
    c = N // 4; bulk = (xy[:, 0] >= c) & (xy[:, 0] < N - c) & (xy[:, 1] >= c) & (xy[:, 1] < N - c)
    acc = []
    for it in range(steps):
        i = rng.integers(0, len(h)); h[i] += 1
        if h[i] > th[i]: h = relax(h)
        if it % 100 == 0: acc.append([((h == 1) & bulk & (sub == s)).sum() / (bulk & (sub == s)).sum() for s in range(nv)])
    a = np.array(acc); return a.mean(0), a.std(0) / np.sqrt(len(a) / 10)
if __name__ == "__main__":
    m, e = run(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 80)
    print(sys.argv[1], "P(h=1) per sublattice:", m, "+-", e)
