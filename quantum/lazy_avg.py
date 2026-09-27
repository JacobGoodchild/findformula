import numpy as np, sys
from lazy import a3
l, T = float(sys.argv[1]), 90
N = 7
w = np.array([1.0] * 6 + [np.sqrt(l)]) / np.sqrt(6 + l)
G = (2 * np.outer(w, w) - np.eye(N)).astype(np.complex64)
L = 2 * T + 3; c = L // 2
a = float(a3(l)); b = (1 - (6 + l) * a) / 6
for start, pred in [(6, l * a * a * (6 + l)), (0, 6 * a * a + 2 * b * b + l * a * a)]:
    psi = np.zeros((N, L, L, L), np.complex64); psi[start, c, c, c] = 1
    rec = []
    for t in range(T):
        psi = np.tensordot(G, psi, axes=(1, 0))
        for j in range(6):
            psi[j] = np.roll(psi[j], 1 if j % 2 == 0 else -1, axis=j // 2)
        rec.append(float(np.sum(np.abs(psi[:, c, c, c])**2)))
    print(f"l={l} start={'loop' if start==6 else '+x'}: time-average t=31..90 = {np.mean(rec[30:]):.6f}   predicted = {pred:.6f}")
