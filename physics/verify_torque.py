"""Stokes flow: ball (radius a) spinning with rate Omega about the line of centres while
touching a fixed sphere of signed radius b (b>0 container, b<0 external ball).
Build the flow from point rotlets (images), check no-slip on both surfaces, and compute the
torque by integrating the viscous shear stress over the ball's surface numerically."""
import numpy as np, sys
from scipy.special import zeta
from numpy.polynomial.legendre import leggauss

a, b, N = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
Om, mu = 1.0, 1.0
cb, Rb = b, abs(b)
pos, G = [a], [Om * a**3]
x, g = a, Om * a**3
for _ in range(N):
    f = x - cb; x = cb + Rb**2 / f; g = -g * (Rb / abs(f))**3; pos.append(x); G.append(g)
    f = x - a;  x = a + a**2 / f;   g = -g * (a / abs(f))**3;  pos.append(x); G.append(g)
pos, G = np.array(pos), np.array(G)

def uphi(xx, rho):
    """azimuthal velocity of axial rotlets: sum G * rho / |r - r0|^3"""
    d2 = (xx[..., None] - pos)**2 + rho[..., None]**2
    return (G * rho[..., None] / d2**1.5).sum(-1)

rng = np.random.default_rng(0)
th = np.arccos(rng.uniform(-1, np.cos(0.2), 3000))
res_ball = np.abs(uphi(a - a * np.cos(th), a * np.sin(th)) - Om * a * np.sin(th)).max()
thb = np.arccos(rng.uniform(-1, np.cos(0.2 * a / Rb), 3000))
xb = cb + (-np.cos(thb) if b > 0 else np.cos(thb)) * Rb
res_fix = np.abs(uphi(xb, Rb * np.sin(thb))).max()
print(f"no-slip residual: ball {res_ball:.1e}, fixed sphere {res_fix:.1e}")

# torque on ball = integral over surface of rho * sigma_{phi n} dS,  sigma_{phi n} = mu*(dU/dn - U n_rho/rho)
t, w = leggauss(3000); th = np.arccos(t)          # theta from the far pole (x=2a side)
nx, nr = np.cos(th), np.sin(th)
X, R = a + a * nx, a * nr
h = 1e-6
dUdn = (uphi(X + h * nx, R + h * nr) - uphi(X - h * nx, R - h * nr)) / (2 * h)
sig = mu * (dUdn - uphi(X, R) * nr / R)
T = -2 * np.pi * a**2 * np.sum(w * R * sig)       # sign: torque resisting the spin reported positive
lam = 1 - a / b
print(f"a={a} b={b}: torque/(8 pi mu Omega a^3) numeric = {T / (8 * np.pi * mu * Om * a**3):.8f}   formula lam^-3 zeta(3,1/lam) = {zeta(3, 1 / lam) / lam**3:.8f}")
