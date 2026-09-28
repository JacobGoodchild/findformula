"""Diamond resistor network + next-nearest (FCC) bonds of conductance w: R_NN reduces to the FCC Green function.
R_NN(w) = ((4w+1)/(8w^2)) P(t_w),  P(t) = E[1/(t - s)], s = c1c2 + c2c3 + c3c1 (cubic half-angle cosines),  t_w = 3 + (8w+1)/(4w^2)."""
import mpmath as mp, numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl, itertools, sys
def P(t, dps=20):
    mp.mp.dps = dps; t = mp.mpf(t)
    def f(k1, k2):   # k3 exactly: s = c1c2 + c3(c1 + c2): E_k3[1/(A - B c3)] = 1/sqrt(A^2 - B^2)
        c1, c2 = mp.cos(k1), mp.cos(k2); A = t - c1 * c2; B = c1 + c2
        return 1 / mp.sqrt(A * A - B * B)
    return mp.quad(f, [0, mp.pi / 2, mp.pi], [0, mp.pi / 2, mp.pi]) / mp.pi**2
def R_formula(w):
    w = mp.mpf(w); return (4 * w + 1) / (8 * w * w) * P(3 + (8 * w + 1) / (4 * w * w))
def R_torus(w, N):
    # FCC primitive vectors a1=(0,1,1)/2, a2=(1,0,1)/2, a3=(1,1,0)/2 in lattice coords; A at R, B at R + (1,1,1)/4.
    # A(R)-B(R - a_j) for a_0 = 0 and a_1..a_3; NNN: same-sublattice bonds along +-a_i, +-(a_i - a_j)
    idx = lambda s, x, y, z: s + 2 * ((x % N) + N * ((y % N) + N * (z % N)))
    n = 2 * N**3; R = []; C = []; V = []
    nn = [(0, 0, 0), (-1, 0, 0), (0, -1, 0), (0, 0, -1)]
    nnn = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, 0), (0, 1, -1), (1, 0, -1)]
    for x, y, z in itertools.product(range(N), repeat=3):
        for o in nn:
            a, b = idx(0, x, y, z), idx(1, x + o[0], y + o[1], z + o[2]); R += [a, b]; C += [b, a]; V += [1.0, 1.0]
        for s in (0, 1):
            for o in nnn:
                a, b = idx(s, x, y, z), idx(s, x + o[0], y + o[1], z + o[2]); R += [a, b]; C += [b, a]; V += [w, w]
    A = sps.csr_matrix((V, (R, C)), shape=(n, n)); L = (sps.diags(np.asarray(A.sum(1)).ravel()) - A).tolil(); L[0, 0] += 1
    lu = spl.splu(L.tocsc()); a, b = idx(0, N // 2, N // 2, N // 2), idx(1, N // 2, N // 2, N // 2)
    r = np.zeros(n); r[a] = 1; r[b] = -1; p = lu.solve(r); return p[a] - p[b]
if __name__ == "__main__":
    for w in map(float, sys.argv[1:]):
        print(w, "formula", R_formula(w), " torus N=12,16,20:", [round(R_torus(w, N), 7) for N in (12, 16, 20)])
