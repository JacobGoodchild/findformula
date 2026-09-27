"""Lackadaisical (lazy) Grover walk on Z^3 with self-loop weight l.
Closed form for the trapped (localized) part, checked against the integral and a simulation."""
import mpmath as mp, numpy as np, sys

def K(b):   # int_0^inf e^{-b x^2} erfc(x)^2 dx, closed form (trivariate orthant probability)
    return mp.sqrt(mp.pi / b) / 2 * (1 + 2 / mp.pi * (mp.asin(1 / (b + 1)) - 2 * mp.asin(1 / mp.sqrt(b + 1))))

def a3(l):  # a(l) = E[1/(6 + l + 2T)] for the 3D lazy walk
    c = mp.mpf(l) / 2
    J = (1 - 3 / mp.sqrt(c + 1) * (1 + 2 / mp.pi * (mp.asin(1 / (c + 2)) - 2 * mp.asin(1 / mp.sqrt(c + 2))))) / c
    return J / 2

if __name__ == "__main__":
    mp.mp.dps = 40
    for b in [mp.mpf(1) / 3, 2, 7]:
        num = mp.quad(lambda x: mp.exp(-b * x * x) * mp.erfc(x)**2, [0, 1, 3, 8, mp.inf])
        print("K check", b, num - K(b))
    for l in [mp.mpf(1) / 2, 1, 3]:
        num = mp.quad(lambda s: mp.exp(-l * s / 2) * mp.erfc(mp.sqrt(s))**3, [0, 1, 5, 20, mp.inf]) / 2
        print("a(l) check", l, num - a3(l))
    # simulation
    l, T = float(sys.argv[1]), int(sys.argv[2])
    d, N = 3, 7
    w = np.array([1.0] * 6 + [np.sqrt(l)]) / np.sqrt(6 + l)
    G = 2 * np.outer(w, w) - np.eye(N)
    L = 2 * T + 3; c = L // 2
    for start in ["loop", "+x"]:
        psi = np.zeros((N, L, L, L), complex)
        psi[6 if start == "loop" else 0, c, c, c] = 1
        for t in range(T):
            psi = np.tensordot(G, psi, axes=(1, 0))
            for j in range(6):
                psi[j] = np.roll(psi[j], 1 if j % 2 == 0 else -1, axis=j // 2)
        P = np.sum(np.abs(psi[:, c, c, c])**2)
        a = float(a3(l)); bb = (1 - (6 + l) * a) / 6
        pred = l * a * a * (6 + l) if start == "loop" else 6 * a * a + 2 * bb * bb + l * a * a
        print(f"l={l} start={start}: simulated P(return, t={T}) = {P:.6f}   predicted limit = {pred:.6f}")
