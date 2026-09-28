"""Lazy flip-flop Grover walk on the DIAMOND lattice: closed form via the FCC Green function (Joyce)."""
import mpmath as mp, sympy as sp, sys
def g_joyce(t):
    x = 3 / t
    k2 = mp.mpf(1) / 2 - (mp.sqrt(3) / 2) * (3 + x)**(-1.5) * (4 * x + (3 - x) * mp.sqrt(1 - x))
    LGF = 3**1.5 * (3 + x)**(-1.5) * (2 - mp.sqrt(1 - x)) * (2 * mp.ellipk(k2) / mp.pi)**2
    return LGF / t
def g_num(t):
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2); A = t - c1 * c2; B = c1 + c2
        return 1 / mp.sqrt(A * A - B * B)
    return mp.quad(f, [0, mp.pi / 2, mp.pi], [0, mp.pi / 2, mp.pi]) / mp.pi**2
def pbar(eps, g=g_joyce):
    eps = mp.mpf(eps); P = (1 - eps) / 4; al = (1 + eps) / 2; be = P / 2
    r = al / be; tau = (r * r - 4) / 4; G = g(tau)
    G00 = al / (4 * be * be) * G
    Gn = -((1 + tau) * G - 1) / (4 * be)
    Gt = al / (4 * be * be) * (tau * G - 1) / 3
    As = (1 - P * (G00 + Gn)) / 2; Ao = -(P / 4) * (G00 + 2 * Gn + Gt); Al2 = eps * P * (G00 + Gn)**2 / 4
    return mp.mpf(1) / 12 + As**2 + 3 * Ao**2 + Al2
if __name__ == "__main__":
    mp.mp.dps = 20
    print("check g:", g_joyce(mp.mpf(5)) - g_num(mp.mpf(5)))
    print("eps->0:", pbar(mp.mpf('1e-15')), "(expect 1/6)")
    for e in ['0.2', '0.5']: print(e, pbar(e), pbar(e, g_num))
    r, gs = sp.symbols('r g')
    eps = (r - 4) / (r + 4); P = (1 - eps) / 4; al = (1 + eps) / 2; be = P / 2; tau = (r * r - 4) / 4
    G00 = al / (4 * be * be) * gs; Gn = -((1 + tau) * gs - 1) / (4 * be); Gt = al / (4 * be * be) * (tau * gs - 1) / 3
    As = (1 - P * (G00 + Gn)) / 2; Ao = -(P / 4) * (G00 + 2 * Gn + Gt); Al2 = eps * P * (G00 + Gn)**2 / 4
    print("As =", sp.factor(sp.simplify(As)), "\nAo =", sp.factor(sp.simplify(Ao)), "\nAloop^2 =", sp.factor(sp.simplify(Al2)))
    print("pbar =", sp.factor(sp.expand(sp.simplify(sp.Rational(1, 12) + As**2 + 3 * Ao**2 + Al2))))
