"""Independent check: Bloch-grid eigenprojectors of the lazy Szegedy walk on diamond (loop arcs included)."""
import numpy as np, itertools, sys
eps = float(sys.argv[1]); N = int(sys.argv[2])
E = [(0, 1, (0, 0, 0)), (0, 1, (-1, 0, 0)), (0, 1, (0, -1, 0)), (0, 1, (0, 0, -1))]
arcs = []
for (u, v, o) in E: arcs.append((u, v, o)); arcs.append((v, u, tuple(-x for x in o)))
arcs += [(0, 0, (0, 0, 0)), (1, 1, (0, 0, 0))]            # loop arcs
n = len(arcs); rev = []
for a in arcs:
    rev.append(arcs.index((a[1], a[0], tuple(-x for x in a[2]))))
C = np.zeros((n, n))
for x in (0, 1):
    idx = [i for i, a in enumerate(arcs) if a[0] == x]
    w = np.array([np.sqrt(eps) if arcs[i][1] == x and arcs[i][2] == (0, 0, 0) and arcs[i][0] == arcs[i][1] else np.sqrt((1 - eps) / 4) for i in idx])
    C[np.ix_(idx, idx)] = 2 * np.outer(w, w) - np.eye(len(idx))
a0 = 0; out0 = [i for i, a in enumerate(arcs) if a[0] == 0]
acc = {1: np.zeros(n, complex), -1: np.zeros(n, complex)}
ks = (np.arange(N) + 0.5) * 2 * np.pi / N
for k in itertools.product(ks, repeat=3):
    S = np.zeros((n, n), complex)
    for i, (u, v, o) in enumerate(arcs):
        S[rev[i], i] = np.exp(-1j * np.dot(k, o))
    w, V = np.linalg.eig(S @ C)
    for lam in (1, -1):
        sel = np.abs(w - lam) < 1e-6
        if sel.any():
            Q, _ = np.linalg.qr(V[:, sel]); acc[lam] += (Q @ Q.conj().T)[:, a0]
tot = sum(np.sum(np.abs(acc[l][out0] / N**3)**2) for l in (1, -1))
print(eps, N, "+1", np.sum(np.abs(acc[1][out0] / N**3)**2), "-1", np.sum(np.abs(acc[-1][out0] / N**3)**2), "total", tot)
