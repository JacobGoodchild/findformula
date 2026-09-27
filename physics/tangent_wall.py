"""Sphere (radius a=1) touching a plane wall, moving normal to it, potential flow.
Tangent-sphere coordinates; reduces to an ODE in the Hankel variable s. Validation case."""
import numpy as np
from scipy.integrate import solve_bvp, quad

m = 0.5  # mu-coordinate of sphere of radius 1
U = 1.0

def rhs(s, y):
    X, Xp = y
    coth = 1 / np.tanh(s * m)
    Xpp = Xp / s - X / s**2 + m**2 * X + m * X * coth / s - U * (1 - 2 * m * s) * np.exp(-s * m)
    return np.vstack([Xp, Xpp])

s0, S = 1e-3, 60.0
def bc(ya, yb):
    # near 0: X ~ c s^2 (+ s^2 log s)  -> use X - s X'/2 ~ 0 approx; at S: decaying ~ exp(-m s)
    return np.array([ya[0] - s0 * ya[1] / 2, yb[1] + m * yb[0]])

s = np.linspace(s0, S, 4000)
sol = solve_bvp(rhs, bc, s, np.zeros((2, s.size)), tol=1e-10, max_nodes=10**6)
print(sol.status, sol.message)
X = lambda t: sol.sol(t)[0]
P1 = lambda t: X(t) / np.tanh(t * m) / t
I = quad(lambda t: P1(t) * (1 - 2 * t * m) * np.exp(-t * m), s0, S, limit=500)[0]
flux = 2 * np.pi / (3 * m) * I
print("C =", abs(flux) / (4 / 3 * np.pi), " (expect ~0.803 = 1.5*zeta(3)-1 =", 1.5 * 1.2020569031595942 - 1, ")")
