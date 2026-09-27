"""Weighted Grover coin (moving shift) on Z^d: check pbar = 2[pi1^2/4 + (p1 A - pi1/2)^2 + p1(1-p1)A^2/2],
A = (1/2) E[min_k Z_k^2/p_k] = int_0^inf prod_k erfc(sqrt(p_k u)) du,  pi1 = P(argmin = 1)."""
import numpy as np, mpmath as mp, itertools
def formula(p):
    mp.mp.dps = 20
    p = [mp.mpf(x) for x in p]
    A = mp.quad(lambda u: mp.fprod(mp.erfc(mp.sqrt(pk * u)) for pk in p), [0, 1, 5, 20, mp.inf])
    # pi1 = int_0^inf e^{-s}/sqrt(pi s) prod_{j>1} erfc(sqrt(p_j s/p_1)) ds
    pi1 = mp.quad(lambda s: mp.exp(-s) / mp.sqrt(mp.pi * s) * mp.fprod(mp.erfc(mp.sqrt(p[j] * s / p[0])) for j in range(1, len(p))), [0, 1, 5, 20, mp.inf])
    return 2 * (pi1**2 / 4 + (p[0] * A - pi1 / 2)**2 + p[0] * (1 - p[0]) * A**2 / 2), A, pi1
def grid(p, N):
    d = len(p); n = 2 * d
    s = np.sqrt(np.repeat(np.array(p) / 2, 2))
    dirs = []
    for i in range(d):
        e = np.zeros(d); e[i] = 1; dirs += [e, -e]
    D = np.array(dirs)
    ks = [(np.arange(N) + 0.37 + 0.11 * i) * 2 * np.pi / N for i in range(d)]
    tot = 0
    for lam in (1, -1):
        M = np.zeros(n, complex)
        for k in itertools.product(*ks):
            th = D @ np.array(k)
            v = s / (1 + lam * np.exp(-1j * th))
            M += v * np.conj(v[0]) / np.sum(np.abs(v)**2)
        tot += np.sum(np.abs(M / N**d)**2)
    return tot
for p, N in [((0.3, 0.7), 300), ((0.5, 0.5), 300), ((0.2, 0.3, 0.5), 50), ((1/3, 1/3, 1/3), 50)]:
    f, A, pi1 = formula(p)
    print(p, "formula", mp.nstr(f, 12), " grid", round(grid(p, N), 8), "  A", mp.nstr(A, 10), " pi1", mp.nstr(pi1, 10))
