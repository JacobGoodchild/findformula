"""Moving-shift Grover walk on 2D lattices: grid estimate of flat-band trapping."""
import numpy as np, sys
def trap(deltas, N=200, start=0):
    n = 2 * len(deltas)
    dirs = []
    for d in deltas: dirs += [d, tuple(-x for x in d)]
    D = np.array(dirs, float)
    G = 2 / n * np.ones((n, n)) - np.eye(n)
    ks = (np.arange(N) + 0.37) * 2 * np.pi / N; ks2 = (np.arange(N) + 0.61) * 2 * np.pi / N
    res = {}
    for lam in (1, -1):
        M = np.zeros(n, complex)
        for k1 in ks:
            for k2 in ks2:
                th = D @ np.array([k1, k2])
                v = 1 / (1 + lam * np.exp(-1j * th))
                M += v * np.conj(v[start]) / np.sum(np.abs(v)**2)
        M /= N * N
        res[lam] = np.sum(np.abs(M)**2)
    return res
if __name__ == "__main__":
    r = trap([(1, 0), (0, 1)]); print("square (check vs Formula 3 d=2: 2A^2+B^2... ):", r, sum(r.values()))
    r = trap([(1, 0), (0, 1), (1, -1)]); print("triangular:", r, sum(r.values()))
