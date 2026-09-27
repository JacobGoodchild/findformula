import mpmath as mp, numpy as np
mp.mp.dps = 30
def a2(l):
    l = mp.mpf(l); q = mp.sqrt(1 + l / 2)
    return (1 - 4 / mp.pi * mp.atan(q) / q) / l
for l in [0.5, 1, 4]:
    num = mp.quad(lambda s: mp.exp(-l * s / 2) * mp.erfc(mp.sqrt(s))**2, [0, 1, 5, 20, mp.inf]) / 2
    print("a2 check", l, num - a2(l))
T = 600; l = 1.0; N = 5
w = np.array([1.0] * 4 + [np.sqrt(l)]) / np.sqrt(4 + l)
G = 2 * np.outer(w, w) - np.eye(N)
L = 2 * T + 3; c = L // 2
a = float(a2(l)); b = (1 - (4 + l) * a) / 4
for start, pred in [(4, l * a * a * (4 + l)), (0, 4 * a * a + 2 * b * b + l * a * a)]:
    psi = np.zeros((N, L, L), complex); psi[start, c, c] = 1; rec = []
    for t in range(T):
        psi = np.tensordot(G, psi, axes=(1, 0))
        for j in range(4):
            psi[j] = np.roll(psi[j], 1 if j % 2 == 0 else -1, axis=j // 2)
        rec.append(np.sum(np.abs(psi[:, c, c])**2))
    print(f"2D l=1 start={'loop' if start==4 else '+x'}: avg t=301..600 = {np.mean(rec[300:]):.6f}  predicted {pred:.6f}")
