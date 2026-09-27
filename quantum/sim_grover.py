"""Brute-force simulation of the d-dimensional Grover walk (moving shift) started at the
origin, compared with the predicted long-time return probability from the exact formula."""
import numpy as np, sys
d, T = int(sys.argv[1]), int(sys.argv[2])
N = 2 * d
L = 2 * T + 3
psi = np.zeros((N,) + (L,) * d, complex)
c = (L // 2,) * d
phi = np.zeros(N); phi[0] = 1.0                       # start in coin state '+x'
psi[(slice(None),) + c] = phi
G = 2.0 / N * np.ones((N, N)) - np.eye(N)
out = []
for t in range(1, T + 1):
    psi = np.tensordot(G, psi, axes=(1, 0))
    for j in range(N):
        axis, step = j // 2, (1 if j % 2 == 0 else -1)
        psi[j] = np.roll(psi[j], step, axis=axis)
    out.append(np.sum(np.abs(psi[(slice(None),) + c])**2))
A = {2: 0.5 - 1 / np.pi, 3: 0.5 - (3 - np.sqrt(3)) / np.pi}[d]
B = (1 - d * A) / d
M = A / 2 * np.ones((N, N))
for i in range(d):
    M[2*i:2*i+2, 2*i:2*i+2] += B / 2 * np.array([[1, -1], [-1, 1]])
pbar = 2 * np.sum((M @ phi)**2)
print(f"d={d}: predicted P(even t -> inf) = {2*pbar:.6f}, time-average = {pbar:.6f}")
for t in [T - 3, T - 2, T - 1, T]:
    print(f"  t={t}: P(return) = {out[t-1]:.6f}")
