"""Projector-grid check of flip-flop Grover trapping on 3D periodic graphs (e.g. diamond).
Averages the eigenprojectors of the Bloch walk operator U(k) over a k-grid."""
import numpy as np, sys, itertools
def pbar(nv, E, origin, start, N):
    arcs = []
    for (u, v, o) in E:
        arcs.append((u, v, o)); arcs.append((v, u, tuple(-x for x in o)))
    n = len(arcs)
    out = {x: [i for i, a in enumerate(arcs) if a[0] == x] for x in range(nv)}
    rev = [arcs.index((a[1], a[0], tuple(-x for x in a[2]))) for a in arcs]
    C = np.zeros((n, n))
    for x, idx in out.items():
        for i in idx:
            for j in idx: C[i, j] = 2 / len(idx) - (i == j)
    a0 = out[origin][start]; acc = {}
    ks = (np.arange(N) + 0.5) * 2 * np.pi / N
    for k in itertools.product(ks, repeat=len(E[0][2])):
        S = np.zeros((n, n), complex)
        for i, (u, v, o) in enumerate(arcs):
            S[rev[i], i] = np.exp(-1j * np.dot(k, o))
        w, V = np.linalg.eig(S @ C)
        ang = np.round(np.mod(np.angle(w) + 1e-4, 2 * np.pi) - 1e-4, 5)
        for key in set(ang):
            Q, _ = np.linalg.qr(V[:, np.abs(ang - key) < 1e-5])
            acc.setdefault(key, np.zeros(n, complex))
            acc[key] += (Q @ Q.conj().T)[:, a0]
    Nd = N ** len(E[0][2])
    res = {key: float(np.sum(np.abs(col[out[origin]] / Nd)**2)) for key, col in acc.items()}
    return {k: v for k, v in res.items() if v > 1e-7}
if __name__ == "__main__":
    # diamond: FCC primitive vectors a1=(0,1,1)/2.., sublattice B at A + (1,1,1)/4; bonds A(R)-B(R-o) for o in {0, a1, a2, a3}
    diamond = [(0, 1, (0, 0, 0)), (0, 1, (-1, 0, 0)), (0, 1, (0, -1, 0)), (0, 1, (0, 0, -1))]
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    r = pbar(2, diamond, 0, 0, N); print(r, "sum", sum(r.values()))
