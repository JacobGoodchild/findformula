import numpy as np, sys
from gen3d import geometric3
def grid3(nv, edges, origin, N):
    arcs = []
    for (u, v, o) in edges:
        arcs.append((u, v, o)); arcs.append((v, u, tuple(-x for x in o)))
    n = len(arcs); rev = [arcs.index((a[1], a[0], tuple(-x for x in a[2]))) for a in arcs]
    C = np.zeros((n, n))
    for x in range(nv):
        idx = [i for i, a in enumerate(arcs) if a[0] == x]
        for i in idx:
            for j in idx: C[i, j] = 2 / len(idx) - (i == j)
    out0 = [i for i, a in enumerate(arcs) if a[0] == origin]
    O = np.array([a[2] for a in arcs], float)
    def U(k):
        S = np.zeros((n, n), complex)
        S[rev, range(n)] = np.exp(1j * O @ k)
        return S @ C
    rng = np.random.default_rng(3)
    ws = [np.linalg.eigvals(U(rng.uniform(-3, 3, 3))) for _ in range(3)]
    flat = [w for w in ws[0] if all(min(abs(x - w)) < 1e-8 for x in ws[1:])]
    uniq = []
    for w in flat:
        if all(abs(w - u) > 1e-6 for u in uniq): uniq.append(w)
    avg = [np.zeros((n, n), complex) for _ in uniq]
    ks = (np.arange(N) + 0.37) * 2 * np.pi / N
    for a in ks:
        for b in ks:
            for c in ks:
                w, V = np.linalg.eig(U(np.array([a, b, c])))
                for i, lam in enumerate(uniq):
                    sel = np.abs(w - lam) < 1e-7
                    if sel.any():
                        Q, _ = np.linalg.qr(V[:, sel]); avg[i] += Q @ Q.conj().T
    a0 = out0[0]
    per = [sum(abs(M[b, a0] / N**3)**2 for b in out0) for M in avg]
    return uniq, per
if __name__ == "__main__":
    fcc = [(0, 2, 2), (2, 0, 2), (2, 2, 0)]
    nv, E = geometric3(fcc, [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)], 2)
    N = int(sys.argv[1])
    uniq, per = grid3(nv, E, 0, N)
    for u, p in zip(uniq, per): print(np.round(u, 5), p)
    print("total", sum(per), " (engine 0.29388359496)")
