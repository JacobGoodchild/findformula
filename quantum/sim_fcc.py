import numpy as np
T = 56; P = 6; N = 12
deltas = [(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
L = 2*T + 3; c = L // 2
G = (np.ones((N, N)) / P - np.eye(N)).astype(np.complex64)
psi = np.zeros((N, L, L, L), np.complex64); psi[0, c, c, c] = 1
rec = []
for t in range(T):
    psi = np.tensordot(G, psi, axes=(1, 0))
    new = np.empty_like(psi)
    for p, d in enumerate(deltas):
        new[2*p+1] = np.roll(psi[2*p], d, axis=(0, 1, 2))
        new[2*p] = np.roll(psi[2*p+1], tuple(-v for v in d), axis=(0, 1, 2))
    psi = new
    rec.append(float(np.sum(np.abs(psi[:, c, c, c])**2)))
print("simulated time-average t=17..56:", np.mean(rec[16:]), "  predicted 0.3789146311")
print("last P(t):", np.round(rec[-4:], 5))
