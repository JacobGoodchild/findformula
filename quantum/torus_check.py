"""Independent check: average the projector onto the eigenvalue +1 eigenvector of the
Fourier-space 3D Grover walk over a k-grid (numerical diagonalisation), compare with
the closed-form matrix M = (A/2) J + (B/2) blockdiag([[1,-1],[-1,1]])."""
import numpy as np
d, n = 3, 64
N = 2 * d
G = 2 / N * np.ones((N, N)) - np.eye(N)
ks = (np.arange(n) + 0.5) * 2 * np.pi / n - np.pi          # midpoint grid, avoids k=0,pi exactly
Msum = np.zeros((N, N), complex)
for kx in ks:
    for ky in ks:
        for kz in ks:
            k = (kx, ky, kz)
            s = np.array([np.exp(1j * k[j // 2] * (1 if j % 2 == 0 else -1)) for j in range(N)])
            w, V = np.linalg.eig(np.diag(s) @ G)
            i = np.argmin(np.abs(w - 1))
            v = V[:, i] / np.linalg.norm(V[:, i])
            Msum += np.outer(v, v.conj())
M = Msum / n**3
A = 0.5 - (3 - np.sqrt(3)) / np.pi; B = (1 - 3 * A) / 3
Mf = A / 2 * np.ones((N, N))
for i in range(d):
    Mf[2*i:2*i+2, 2*i:2*i+2] += B / 2 * np.array([[1, -1], [-1, 1]])
print("numeric M row 0:", np.round(M[0].real, 6))
print("formula M row 0:", np.round(Mf[0], 6))
print("max |difference|:", np.abs(M - Mf).max())
