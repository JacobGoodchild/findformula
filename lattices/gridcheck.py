"""Independent check: diagonalise the Bloch flip-flop Grover operator on an N x N momentum grid,
group eigenvalues that are the same at every k (flat bands), average their projectors, and
return the trapping probability for each start arc at the origin vertex."""
import numpy as np, sys
from builders import LATTICES
def gridcheck(nv, edges, origin, N=120):
    arcs = []
    for (u, v, o) in edges:
        arcs.append((u, v, o)); arcs.append((v, u, (-o[0], -o[1])))
    n = len(arcs)
    rev = [arcs.index((a[1], a[0], (-a[2][0], -a[2][1]))) for a in arcs]
    C = np.zeros((n, n))
    for x in range(nv):
        idx = [i for i, a in enumerate(arcs) if a[0] == x]
        for i in idx:
            for j in idx: C[i, j] = 2 / len(idx) - (i == j)
    out0 = [i for i, a in enumerate(arcs) if a[0] == origin]
    ks = (np.arange(N) + 0.37) * 2 * np.pi / N
    # find flat eigenvalues from two random k points
    def U(k1, k2):
        S = np.zeros((n, n), complex)
        for i, (u, v, o) in enumerate(arcs):
            S[rev[i], i] = np.exp(1j * (k1 * o[0] + k2 * o[1]))
        return S @ C
    w1 = np.linalg.eigvals(U(0.3, 1.1)); w2 = np.linalg.eigvals(U(2.2, -0.7)); w3 = np.linalg.eigvals(U(-1.3, 2.9))
    flat = []
    for w in w1:
        if min(abs(w2 - w)) < 1e-8 and min(abs(w3 - w)) < 1e-8 and all(abs(w - f) > 1e-6 for f in flat):
            flat.append(w)
    avg = {i: np.zeros((n, n), complex) for i in range(len(flat))}
    for k1 in ks:
        for k2 in ks:
            w, V = np.linalg.eig(U(k1, k2))
            for i, lam in enumerate(flat):
                sel = np.abs(w - lam) < 1e-7
                if sel.any():
                    Q, _ = np.linalg.qr(V[:, sel]); avg[i] += Q @ Q.conj().T
    res = {}
    for s, a in enumerate(out0):
        res[(arcs[a][1], arcs[a][2])] = sum(sum(abs(avg[i][b, a] / N**2)**2 for b in out0) for i in avg)
    return flat, res
if __name__ == "__main__":
    name, origin = sys.argv[1], int(sys.argv[2])
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 120
    nv, E = LATTICES[name]()
    flat, res = gridcheck(nv, E, origin, N)
    print(name, "flat eigenvalues:", np.round(flat, 6))
    for k, v in res.items(): print("  start arc to", k, ":", round(v, 8))
