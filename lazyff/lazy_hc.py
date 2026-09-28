"""Lazy flip-flop Grover (Szegedy) walk on the honeycomb lattice. Stay with prob eps, else move to one of
3 neighbours. -1 band via GQ = ((I + P_walk)/2)^(-1); +1 band unchanged (= 1/24)."""
import mpmath as mp, numpy as np
mp.mp.dps = 20
# sublattice A at 0, B neighbours of A(0,0): B(0,0), B(-1,0), B(0,-1); Bloch f(k) = 1 + e^{-ik1} + e^{-ik2}
NB = [(0, 0), (-1, 0), (0, -1)]
def GQ_entries(eps):
    eps = mp.mpf(eps); P = (1 - eps) / 3; al = (1 + eps) / 2; be = P / 2
    # (I + Pw)/2 = [[al, be f],[be f*, al]] ; inverse = [[al, -be f],[-be f*, al]] / (al^2 - be^2 |f|^2)
    def E(fun):
        g = lambda k1, k2: fun(k1, k2)
        return mp.quad(lambda k1: mp.quad(lambda k2: g(k1, k2), [-mp.pi, 0, mp.pi]), [-mp.pi, 0, mp.pi]) / (4 * mp.pi**2)
    def den(k1, k2):
        f = 1 + mp.expj(-k1) + mp.expj(-k2); return al * al - be * be * abs(f)**2, f
    # G(A0, A0)
    def gAA(k1, k2): d, f = den(k1, k2); return al / d
    # G(B_R, A0) = E[-be f* e^{ik.R}/d]... use real parts by symmetry
    def gBA(R):
        def h(k1, k2):
            d, f = den(k1, k2); return mp.re(-be * mp.conj(f) * mp.expj(k1 * R[0] + k2 * R[1]) / d)
        return E(h)
    def gBB(dR):
        def h(k1, k2):
            d, f = den(k1, k2); return mp.re(al * mp.expj(k1 * dR[0] + k2 * dR[1]) / d)
        return E(h)
    return P, E(gAA), gBA, gBB
def pbar(eps):
    P, g00, gBA, gBB = GQ_entries(eps)
    ha = NB[0]
    gAa = gBA(ha)
    amps = []
    for hb in NB:
        d = 1 if hb == ha else 0
        amps.append((d - P / 2 * (g00 + gBA(hb) + gAa + gBB((hb[0] - ha[0], hb[1] - ha[1])))) / 2)
    aloop = -mp.sqrt(eps) * mp.sqrt(P / 2) * (g00 + gAa) / mp.sqrt(2)
    pm = mp.fsum(a * a for a in amps) + aloop**2
    return mp.mpf(1) / 24 + pm, pm, amps, aloop
def sim(eps, T=300):
    # arcs 0,1,2 at A to B(R+NB), arcs 3,4,5 at B to A; coin 4 states per site (3 arcs + loop)
    w = np.sqrt([(1 - eps) / 3] * 3 + [eps]); C = 2 * np.outer(w, w) - np.eye(4)
    L = 2 * T + 5; c = L // 2
    A = np.zeros((4, L, L), complex); B = np.zeros((4, L, L), complex); A[0, c, c] = 1; rec = []
    for _ in range(T):
        A = np.tensordot(C, A, axes=(1, 0)); B = np.tensordot(C, B, axes=(1, 0))
        nA = np.zeros_like(A); nB = np.zeros_like(B)
        for i, o in enumerate(NB):
            nB[i] = np.roll(A[i], o, axis=(0, 1))                  # A(R) arc i -> B(R+o), lands on B's arc i (back to A)
            nA[i] = np.roll(B[i], (-o[0], -o[1]), axis=(0, 1))
        nA[3] = A[3]; nB[3] = B[3]; A, B = nA, nB
        rec.append(np.sum(np.abs(A[:, c, c])**2))
    return np.mean(rec[T // 2:])
if __name__ == "__main__":
    print("eps->0:", pbar(mp.mpf('1e-25'))[0], "(expect 1/12)")
    for e in (0.3, 0.6):
        print(e, pbar(e)[0], "sim", sim(e))
