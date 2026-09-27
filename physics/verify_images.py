"""Brute-force check of the image solution: build phi from N image doublets, test the
boundary conditions at random surface points, and integrate the added mass directly."""
import numpy as np, sys
from scipy.special import zeta

a, b, N = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
U = 1.0
beta = b - a
# positions along axis x (contact point at x=0, centres at x=a and x=b)
pos, strength = [], []
A0 = -U * a**3 / 2
for k in range(1, N + 1):          # inner images (k=1 is the sphere centre)
    Ak = A0 * b**3 / (a + k * beta)**3
    wk = 1 / a + (k - 1) * (1 / a - 1 / b)
    pos.append(1 / wk); strength.append(Ak)
    # image of this doublet in the container (outside container, x<0)
    f = b - 1 / wk
    pos.append(b - b**2 / f); strength.append(-Ak * (b / f)**3)
pos, strength = np.array(pos), np.array(strength)

def grad(P):
    """gradient of sum of x-doublets  A (x-x0)/|r-r0|^3  at points P (n,3)"""
    R = P[:, None, :] - np.stack([pos, 0 * pos, 0 * pos], 1)[None]
    r2 = (R**2).sum(-1); r3 = r2**1.5; r5 = r2**2.5
    gx = strength * (1 / r3 - 3 * R[..., 0]**2 / r5)
    gy = strength * (-3 * R[..., 0] * R[..., 1] / r5)
    gz = strength * (-3 * R[..., 0] * R[..., 2] / r5)
    return np.stack([gx.sum(1), gy.sum(1), gz.sum(1)], 1)

def phi(P):
    R = P[:, None, :] - np.stack([pos, 0 * pos, 0 * pos], 1)[None]
    r3 = ((R**2).sum(-1))**1.5
    return (strength * R[..., 0] / r3).sum(1)

rng = np.random.default_rng(1)
def sphere_pts(c, R, n, min_angle):
    # points on sphere (centre (c,0,0), radius R), excluding a cap of half-angle min_angle round contact point
    th = np.arccos(rng.uniform(-1, np.cos(min_angle), n)); ph = rng.uniform(0, 2 * np.pi, n)
    nrm = np.stack([-np.cos(th), np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph)], 1)  # angle measured from -x (contact side)
    return np.array([c, 0, 0]) + R * nrm, nrm

for cap in [0.3, 0.1]:
    P, n_ = sphere_pts(a, a, 2000, cap)
    res_in = np.abs((grad(P) * n_).sum(1) - U * n_[:, 0]).max()
    P, n_ = sphere_pts(b, b, 2000, cap * a / b)
    res_out = np.abs((grad(P) * n_).sum(1)).max()
    print(f"N={N} cap={cap}: max BC residual inner={res_in:.2e}  outer={res_out:.2e}")

# added mass by direct surface quadrature on inner sphere: M = -(rho/U) * surface integral of phi*n_x
from numpy.polynomial.legendre import leggauss
t, w = leggauss(4000)                      # t = cos(theta) with theta from the far pole
th = np.arccos(t)
P = np.stack([a + a * np.cos(th), a * np.sin(th), 0 * th], 1)   # axisymmetric: phi independent of azimuth
integral = 2 * np.pi * a**2 * np.sum(w * phi(P) * np.cos(th))
C_num = -integral / U / (4 / 3 * np.pi * a**3)
q = b / (b - a)
print("C direct quadrature:", C_num, "  formula:", 1.5 * q**3 * zeta(3, q) - 1)
