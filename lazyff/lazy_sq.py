"""Lazy flip-flop Grover (Szegedy) walk on Z^2: P = (1-eps) P_simple + eps (self-loop)."""
import mpmath as mp, numpy as np
mp.mp.dps = 25
def g(x, t):   # E[cos(k.x)/(t + c1 + c2)], k2 integrated exactly
    m1, m2 = x
    def f(k1):
        A = t + mp.cos(k1); root = mp.sqrt(A * A - 1)
        return mp.cos(m1 * k1) * (-(A - root))**abs(m2) / root     # E_k2[cos(m k2)/(A + cos k2)] = (-r)^|m|/root
    return mp.quad(f, [0, mp.pi]) / mp.pi
def pbar(eps):
    eps = mp.mpf(eps); P = (1 - eps) / 4; t = 2 + 4 * eps / (1 - eps)
    GQ = lambda x: 4 / (1 - eps) * g(x, t)
    ha = (1, 0)
    heads = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    amps = []
    for hb in heads:
        d = (1 if hb == ha else 0)
        amps.append((d - P / 2 * (GQ((0, 0)) + GQ(ha) + GQ(hb) + GQ((hb[0] - ha[0], hb[1] - ha[1])))) / 2)
    aloop = -mp.sqrt(eps) * mp.sqrt(P / 2) * (GQ((0, 0)) + GQ(ha)) / mp.sqrt(2)
    pm = mp.fsum(a * a for a in amps) + aloop**2
    pp = ((mp.mpf(1)/2)**2 + (2/mp.pi - mp.mpf(1)/2)**2 + 2 * (mp.mpf(1)/2 - 1/mp.pi)**2) / 4
    return pp + pm, pp, pm
def sim(eps, T=400):
    w = np.sqrt([(1 - eps) / 4] * 4 + [eps]); C = 2 * np.outer(w, w) - np.eye(5)
    L = 2 * T + 3; c = L // 2
    psi = np.zeros((5, L, L), complex); psi[0, c, c] = 1; rec = []
    for _ in range(T):
        psi = np.tensordot(C, psi, axes=(1, 0)); new = np.empty_like(psi)
        new[1] = np.roll(psi[0], 1, axis=0); new[0] = np.roll(psi[1], -1, axis=0)
        new[3] = np.roll(psi[2], 1, axis=1); new[2] = np.roll(psi[3], -1, axis=1); new[4] = psi[4]
        psi = new; rec.append(np.sum(np.abs(psi[:, c, c])**2))
    return np.mean(rec[T // 2:])
if __name__ == "__main__":
    for e in [0.2, 0.5]:
        v = pbar(e); print("eps", e, "formula", v[0], "(+1", v[1], "-1", v[2], ")  sim", sim(e))
