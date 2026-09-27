import numpy as np, sys
def sim(deltas, T, start=0):
    P = len(deltas); N = 2 * P
    L = 2 * T + 3; c = L // 2
    G = np.ones((N, N)) / P - np.eye(N)
    psi = np.zeros((N, L, L), complex); psi[start, c, c] = 1
    rec = []
    for t in range(T):
        psi = np.tensordot(G, psi, axes=(1, 0))
        new = np.empty_like(psi)
        for p, d in enumerate(deltas):
            new[2*p+1] = np.roll(psi[2*p], d, axis=(0, 1))
            new[2*p] = np.roll(psi[2*p+1], (-d[0], -d[1]), axis=(0, 1))
        psi = new
        rec.append(np.sum(np.abs(psi[:, c, c])**2))
    return np.array(rec)
if __name__ == "__main__":
    r = sim([(1, 0), (0, 1), (1, -1)], 500)
    print("triangular: time-average t=201..500:", r[200:].mean(), " predicted 0.2679557025")
    r = sim([(1, 0), (0, 1)], 500)
    print("square: time-average t=201..500:", r[200:].mean(), " predicted 0.1673437786")
