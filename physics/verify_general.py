"""Generic image iteration (no closed form assumed) for a moving sphere (radius a, centre x=a)
touching a fixed sphere of signed radius b at the origin:
  b>0: container (centre x=b, fluid inside it);  b<0: external sphere radius |b| (centre x=b<0).
Checks boundary conditions and integrates the added mass directly."""
import numpy as np, sys
from scipy.special import zeta
from numpy.polynomial.legendre import leggauss

a, b, N = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
U = 1.0
cb, Rb = b, abs(b)
pos, st = [a], [-U * a**3 / 2]
x, A = a, -U * a**3 / 2
for _ in range(N):
    f = x - cb; x = cb + Rb**2 / f; A = -A * (Rb / abs(f))**3; pos.append(x); st.append(A)   # image in fixed sphere
    f = x - a;  x = a + a**2 / f;   A = -A * (a / abs(f))**3;  pos.append(x); st.append(A)    # image in moving sphere
pos, st = np.array(pos), np.array(st)

def field(P, deriv):
    R = P[:, None, :] - np.stack([pos, 0 * pos, 0 * pos], 1)[None]
    r2 = (R**2).sum(-1)
    if not deriv:
        return (st * R[..., 0] / r2**1.5).sum(1)
    g = [st * (1 / r2**1.5 - 3 * R[..., 0]**2 / r2**2.5), st * (-3 * R[..., 0] * R[..., 1] / r2**2.5), st * (-3 * R[..., 0] * R[..., 2] / r2**2.5)]
    return np.stack([gi.sum(1) for gi in g], 1)

rng = np.random.default_rng(0)
def pts(c, R, n, cap, sgn):
    th = np.arccos(rng.uniform(-1, np.cos(cap), n)); ph = rng.uniform(0, 2 * np.pi, n)
    nr = np.stack([sgn * np.cos(th), np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph)], 1)
    return np.array([c, 0, 0]) + R * nr, nr
P, nr = pts(a, a, 2000, 0.2, -1)
print("inner BC residual:", np.abs((field(P, 1) * nr).sum(1) - U * nr[:, 0]).max())
P, nr = pts(cb, Rb, 2000, 0.2 * a / Rb, 1 if b < 0 else -1)
print("fixed-sphere BC residual:", np.abs((field(P, 1) * nr).sum(1)).max())
t, w = leggauss(4000); th = np.arccos(t)
P = np.stack([a + a * np.cos(th), a * np.sin(th), 0 * th], 1)
C = -2 * np.pi * a**2 * np.sum(w * field(P, 0) * np.cos(th)) / U / (4 / 3 * np.pi * a**3)
lam = 1 - a / b
print(f"a={a} b={b}: C numeric = {C:.9f}   formula (3/2) lam^-3 zeta(3,1/lam) - 1 = {1.5 / lam**3 * zeta(3, 1 / lam) - 1:.9f}")
