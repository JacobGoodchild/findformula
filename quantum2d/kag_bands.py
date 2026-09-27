"""Numerical check: Bloch diagonalisation of the kagome flip-flop Grover walk; average the
projectors of each flat eigenvalue (+1, -1, e^{+-2pi i/3}) over the Brillouin zone."""
import numpy as np
kag = [(0, 1, (0, 0)), (0, 2, (0, 0)), (1, 2, (0, 0)), (1, 0, (1, 0)), (2, 0, (0, 1)), (1, 2, (1, -1))]
arcs = []
for (u, v, o) in kag:
    arcs.append((u, v, o)); arcs.append((v, u, (-o[0], -o[1])))
n = len(arcs)
rev = [arcs.index((a[1], a[0], (-a[2][0], -a[2][1]))) for a in arcs]
out0 = [i for i, a in enumerate(arcs) if a[0] == 0]
C = np.zeros((n, n))
for x in range(3):
    idx = [i for i, a in enumerate(arcs) if a[0] == x]
    for i in idx:
        for j in idx:
            C[i, j] = 2 / len(idx) - (i == j)
targets = {"+1": 1, "-1": -1, "w": np.exp(2j*np.pi/3), "w*": np.exp(-2j*np.pi/3)}
avg = {k: np.zeros((n, n), complex) for k in targets}
N = 90
for k1 in (np.arange(N) + 0.5) * 2 * np.pi / N:
    for k2 in (np.arange(N) + 0.5) * 2 * np.pi / N:
        S = np.zeros((n, n), complex)
        for i, (u, v, o) in enumerate(arcs):
            S[rev[i], i] = np.exp(1j * (k1 * o[0] + k2 * o[1]))
        U = S @ C
        w, V = np.linalg.eig(U)
        for name, lam in targets.items():
            sel = np.abs(w - lam) < 1e-6
            if sel.any():
                Q, _ = np.linalg.qr(V[:, sel])
                avg[name] += Q @ Q.conj().T
a = out0[0]
tot = 0
for name in targets:
    M = avg[name] / N**2
    p = sum(abs(M[b, a])**2 for b in out0)
    print(name, "contribution", p, " rank/k ~", np.real(np.trace(M)))
    tot += p
print("total", tot, " (simulation 0.19099)")
