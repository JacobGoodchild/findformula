"""Independent check: sphere radius a touching inside a container radius b, potential flow,
solved as a boundary-value problem in tangent-sphere/Hankel variables (no images used)."""
import sys
import numpy as np
from scipy.integrate import solve_bvp, quad
from scipy.special import zeta

def added_mass(a, b, s0, S=80.0):
    m1, m2 = 1 / (2 * a), 1 / (2 * b)
    d = m1 - m2
    U = 1.0
    def f(s, y):
        X, Xp, Y, Yp = y
        cth, csh = 1 / np.tanh(s * d), 1 / np.sinh(s * d)
        Xpp = Xp / s - X / s**2 + m1**2 * X + m1 * (X * cth - Y * csh) / s - U * (1 - 2 * m1 * s) * np.exp(-s * m1)
        Ypp = Yp / s - Y / s**2 + m2**2 * Y + m2 * (X * csh - Y * cth) / s
        return np.vstack([Xp, Xpp, Yp, Ypp])
    def bc(ya, yb):
        X, Xp, Y, Yp = ya
        al, alp = (X - Y) / (m1 - m2), (Xp - Yp) / (m1 - m2)
        ga, gap = (m1 * Y - m2 * X) / (m1 - m2), (m1 * Yp - m2 * Xp) / (m1 - m2)
        return np.array([al - s0 * alp / 2, ga - s0 * gap, yb[1] + m1 * yb[0], yb[3] + m2 * yb[2]])
    s = np.geomspace(s0, S, 6000)
    sol = solve_bvp(f, bc, s, np.zeros((4, s.size)), tol=1e-10, max_nodes=2 * 10**6)
    X = lambda t: sol.sol(t)[0]; Y = lambda t: sol.sol(t)[2]
    P1 = lambda t: (X(t) * np.cosh(t * d) - Y(t)) / (t * np.sinh(t * d))
    pts = list(np.geomspace(s0, S, 40))
    I = sum(quad(lambda t: P1(t) * (1 - 2 * t * m1) * np.exp(-t * m1), pts[i], pts[i + 1], limit=200)[0] for i in range(len(pts) - 1))
    flux = 2 * np.pi / (3 * m1) * I
    return abs(flux) / (4 / 3 * np.pi * a**3), sol.status

def formula(r):
    q = 1 / (1 - r)
    return 1.5 * q**3 * zeta(3, q) - 1

if __name__ == "__main__":
    a, b = float(sys.argv[1]), float(sys.argv[2])
    for s0 in [1e-3, 1e-4, 1e-5, 1e-6]:
        print(s0, added_mass(a, b, s0), " formula:", formula(a / b))
