import numpy as np
from lazy import a3
l, T = 1.0, 90
N = 7
w = np.array([1.0] * 6 + [np.sqrt(l)]) / np.sqrt(6 + l)
G = 2 * np.outer(w, w) - np.eye(N)
L = 2 * T + 3; c = L // 2
psi = np.zeros((N, L, L, L), np.complex64); psi[6, c, c, c] = 1
a = float(a3(l)); pred = l * a * a * (6 + l)
rec = []
for t in range(1, T + 1):
    psi = np.tensordot(G.astype(np.complex64), psi, axes=(1, 0))
    for j in range(6):
        psi[j] = np.roll(psi[j], 1 if j % 2 == 0 else -1, axis=j // 2)
    rec.append(float(np.sum(np.abs(psi[:, c, c, c])**2)))
for t in [10, 20, 30, 45, 60, 75, 90]:
    print(t, round(rec[t - 1], 6), " error vs prediction:", round(rec[t - 1] - pred, 6), " error*t^1.5:", round((rec[t-1]-pred)*t**1.5, 4))
print("prediction", pred)
