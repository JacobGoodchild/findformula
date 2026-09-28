"""Lazy flip-flop Grover (Szegedy) walk on the triangular lattice: P = (1-eps) P_simple + eps (self-loop).
Mirrors lazy_sq.py. -1 band uses the signless Green function E[cos(k.x)/(t + S)], t = 3(1+eps)/(1-eps),
S = cos k1 + cos k2 + cos(k1 - k2)."""
import mpmath as mp, numpy as np
mp.mp.dps = 25
HEADS = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
def g(x, t):
    m1, m2 = x
    def f(k1):
        # E_k2[cos(m1 k1 + m2 k2)/(t + cos k1 + cos k2 + cos(k1 - k2))]; cos k2 + cos(k1-k2) = R cos(k2 - k1/2), R = 2cos(k1/2)
        A = t + mp.cos(k1); R = 2 * mp.cos(k1 / 2)
        root = mp.sqrt(A * A - R * R)
        r = (A - root) / R if R != 0 else mp.mpf(0)
        # E[cos(m2 (phi + k1/2))/(A + R cos phi)] = (-r)^|m2| cos(m2 k1/2)/root
        return mp.cos(m1 * k1 + m2 * k1 / 2) * (-r)**abs(m2) / root
    return mp.quad(f, [0, mp.pi / 2, mp.pi]) / mp.pi
def pbar(eps):
    eps = mp.mpf(eps); P = (1 - eps) / 6; t = 3 * (1 + eps) / (1 - eps)
    GQ = lambda x: 6 / (1 - eps) * g(x, t)
    ha = HEADS[0]; amps = []
    for hb in HEADS:
        d = 1 if hb == ha else 0
        amps.append((d - P / 2 * (GQ((0, 0)) + GQ(ha) + GQ(hb) + GQ((hb[0] - ha[0], hb[1] - ha[1])))) / 2)
    aloop = -mp.sqrt(eps) * mp.sqrt(P / 2) * (GQ((0, 0)) + GQ(ha)) / mp.sqrt(2)
    pm = mp.fsum(a * a for a in amps) + aloop**2
    pp = mp.mpf("0.13428605624308877787")          # triangular +1 channel (Formula 7), unchanged by laziness
    return pp + pm, pp, pm, amps, aloop
def sim(eps, T=300):
    w = np.sqrt([(1 - eps) / 6] * 6 + [eps]); C = 2 * np.outer(w, w) - np.eye(7)
    L = 2 * T + 3; c = L // 2
    psi = np.zeros((7, L, L), complex); psi[0, c, c] = 1; rec = []
    rev = [1, 0, 3, 2, 5, 4]
    for _ in range(T):
        psi = np.tensordot(C, psi, axes=(1, 0)); new = np.empty_like(psi)
        for i, h in enumerate(HEADS):   # amplitude on arc i at x moves to x+h and becomes the reverse arc
            new[rev[i]] = np.roll(psi[i], h, axis=(0, 1))
        new[6] = psi[6]; psi = new; rec.append(np.sum(np.abs(psi[:, c, c])**2))
    return np.mean(rec[T // 2:])
if __name__ == "__main__":
    print("eps=0 check (expect triangular total 0.26795570253708055129):", pbar(1e-30)[0])
    for e in [0.2, 0.5]:
        v = pbar(e); print("eps", e, "formula", v[0], "(+1", v[1], "-1", v[2], ")  sim", sim(e))
