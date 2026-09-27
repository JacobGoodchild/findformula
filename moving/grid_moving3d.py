import numpy as np, itertools, sys
def trap3(deltas, N):
    dirs = []
    for d in deltas: dirs += [d, tuple(-x for x in d)]
    D = np.array(dirs, float); n = len(dirs)
    ks = [(np.arange(N) + off) * 2 * np.pi / N for off in (0.37, 0.61, 0.83)]
    out = {}
    for lam in (1, -1):
        M = np.zeros(n, complex)
        for k in itertools.product(*ks):
            th = D @ np.array(k)
            v = 1 / (1 + lam * np.exp(-1j * th))
            M += v * np.conj(v[0]) / np.sum(np.abs(v)**2)
        out[lam] = (M / N**3)
    return out
bcc = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]
fcc = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
for name, dl in [("BCC", bcc), ("FCC", fcc)]:
    for N in (24, 40):
        o = trap3(dl, N)
        print(name, N, {l: round(float(np.sum(np.abs(M)**2)), 7) for l, M in o.items()}, "total", round(float(sum(np.sum(np.abs(M)**2) for M in o.values())), 7))
        print("    +1 amps", np.round(o[1].real, 5), "\n    -1 amps", np.round(o[-1].real, 5))
