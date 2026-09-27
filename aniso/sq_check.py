import mpmath as mp, numpy as np
from sq import trap
mp.mp.dps = 25
def closed(al):
    al = mp.mpf(al); be = 1 - al; th = mp.asin(mp.sqrt(al))
    Yxx = 2 * th / mp.pi
    Yxmx = 2 / mp.pi * (mp.sqrt(be / al) - th * be / al)
    Yxy = mp.sqrt(al / be) * (1 - Yxx - Yxmx) / 2
    return ((1 - Yxx)**2 + Yxmx**2 + 2 * Yxy**2) / 2
for al in ["0.5", "0.25", "0.2", "0.7"]:
    print(al, "closed", closed(al), " numeric", trap(mp.mpf(al))[0])
def sim(al, T=400):
    be = 1 - al
    w = np.sqrt([al / 2, al / 2, be / 2, be / 2]); C = 2 * np.outer(w, w) - np.eye(4)
    L = 2 * T + 3; c = L // 2
    psi = np.zeros((4, L, L), complex); psi[0, c, c] = 1
    rec = []
    for t in range(T):
        psi = np.tensordot(C, psi, axes=(1, 0))
        new = np.empty_like(psi)
        new[1] = np.roll(psi[0], 1, axis=0); new[0] = np.roll(psi[1], -1, axis=0)   # +x hop lands on -x arc
        new[3] = np.roll(psi[2], 1, axis=1); new[2] = np.roll(psi[3], -1, axis=1)
        psi = new; rec.append(np.sum(np.abs(psi[:, c, c])**2))
    return np.mean(rec[T // 2:])
for al in [0.25, 0.7]:
    print("sim alpha", al, sim(al), " closed", closed(al))
