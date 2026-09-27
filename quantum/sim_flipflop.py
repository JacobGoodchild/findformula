import numpy as np, mpmath as mp
T = 80; d = 3; N = 6
L = 2 * T + 3; c = L // 2
G = np.ones((N, N)) / d - np.eye(N)
psi = np.zeros((N, L, L, L), np.complex64); psi[0, c, c, c] = 1     # start: coin +x
rec = []
for t in range(T):
    psi = np.tensordot(G.astype(np.complex64), psi, axes=(1, 0))
    new = np.empty_like(psi)
    for i in range(d):
        new[2*i+1] = np.roll(psi[2*i], 1, axis=i)    # hop +e_i, coin +i -> -i (points back)
        new[2*i] = np.roll(psi[2*i+1], -1, axis=i)   # hop -e_i, coin -i -> +i
    psi = new
    rec.append(float(np.sum(np.abs(psi[:, c, c, c])**2)))
mp.mp.dps = 30
W = mp.sqrt(6) / (32 * mp.pi**3) * mp.gamma(mp.mpf(1)/24) * mp.gamma(mp.mpf(5)/24) * mp.gamma(mp.mpf(7)/24) * mp.gamma(mp.mpf(11)/24)
X = 54 / (mp.pi**2 * W)
al = (7 * W + X - 12) / 36; be = (24 - 7 * W - X) / 144
Mphi2 = (mp.mpf(1)/3)**2 + al**2 + 4 * be**2
print("predicted time-average:", 2 * Mphi2, "  simulated avg t=21..80:", np.mean(rec[20:]))
print("last few P(t):", np.round(rec[-4:], 5))
