import numpy as np
from honey_formula import pbar
def sim(p, T=400):
    # honeycomb, A(R) -- B(R + o) for o in offs with weights p (same weight seen from both ends)
    offs = [(0, 0), (-1, 0), (0, -1)]
    arcs = []   # (tail sublattice, head sublattice, offset, weight index)
    for t, o in enumerate(offs):
        arcs.append((0, 1, o, t)); arcs.append((1, 0, (-o[0], -o[1]), t))
    n = len(arcs); L = 2 * T + 5; c = L // 2
    rev = [arcs.index((a[1], a[0], (-a[2][0], -a[2][1]), a[3])) for a in arcs]
    psi = np.zeros((n, L, L), complex)
    start = [i for i, a in enumerate(arcs) if a[0] == 0 and a[3] == 0][0]
    psi[start, c, c] = 1
    rec = []
    for _ in range(T):
        new = np.zeros_like(psi)
        for x in (0, 1):
            idx = [i for i, a in enumerate(arcs) if a[0] == x]
            w = np.sqrt([p[arcs[i][3]] for i in idx])
            s = sum(w[k] * psi[i] for k, i in enumerate(idx))
            for k, i in enumerate(idx): new[i] = 2 * w[k] * s - psi[i]
        psi2 = np.zeros_like(psi)
        for i, a in enumerate(arcs): psi2[rev[i]] = np.roll(new[i], a[2], axis=(0, 1))
        psi = psi2
        rec.append(sum(np.abs(psi[i, c, c])**2 for i, a in enumerate(arcs) if a[0] == 0))
    return np.mean(rec[T // 2:])
for p in [(0.5, 0.25, 0.25), (0.6, 0.25, 0.15)]:
    print(p, "sim", sim(p), " formula", pbar(*p)[0])
