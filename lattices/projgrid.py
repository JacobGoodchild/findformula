"""Independent check of flip-flop trapping: diagonalise the Bloch walk operator U(k) on a grid,
average the flat-band eigenprojectors over the Brillouin zone, sum squared return amplitudes."""
import numpy as np, sys
from builders import LATTICES
def pbar(nv, E, origin, start, N, tol=1e-7):
    arcs = []
    for (u, v, o) in E:
        arcs.append((u, v, o)); arcs.append((v, u, (-o[0], -o[1])))
    n = len(arcs)
    out = {x: [i for i, a in enumerate(arcs) if a[0] == x] for x in range(nv)}
    rev = [arcs.index((a[1], a[0], (-a[2][0], -a[2][1]))) for a in arcs]
    C = np.zeros((n, n))
    for x, idx in out.items():
        d = len(idx)
        for i in idx:
            for j in idx: C[i, j] = 2 / d - (i == j)
    a0 = out[origin][start]
    acc = {}
    ks = (np.arange(N) + 0.5) * 2 * np.pi / N
    for k1 in ks:
        for k2 in ks:
            S = np.zeros((n, n), complex)
            for i, (u, v, o) in enumerate(arcs):   # arc i in cell R -> arc rev[i] in cell R + o
                S[rev[i], i] = np.exp(-1j * (k1 * o[0] + k2 * o[1]))
            U = S @ C
            w, V = np.linalg.eig(U)
            # group eigenvalues; flat bands are the same at every k, so key by rounded angle
            ang = np.round(np.mod(np.angle(w) + 1e-4, 2 * np.pi) - 1e-4, 5)
            for key in set(ang):
                sel = np.abs(ang - key) < 1e-5
                Q, _ = np.linalg.qr(V[:, sel])
                P = Q @ Q.conj().T
                acc.setdefault(key, np.zeros(n, complex))
                acc[key] += P[:, a0]
    res = {}
    for key, col in acc.items():
        col /= N * N
        res[key] = np.sum(np.abs(col[out[origin]])**2)
    # only flat bands contribute a nonzero BZ average of order 1; dispersive keys are sampled once each
    return res
if __name__ == "__main__":
    name, origin, start, N = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    nv, E = LATTICES[name]()
    r = pbar(nv, E, origin, start, N)
    big = {k: v for k, v in r.items() if v > 1e-6}
    print(big, "sum", sum(big.values()))
