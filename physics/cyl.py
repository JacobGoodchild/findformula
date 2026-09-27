"""Wall-correction coefficient k for a sphere moving along the axis of a circular pipe
(Stokes flow): drag = 6*pi*mu*a*U*(1 + k a/R + ...). Derived by reflecting a Stokeslet
off the no-slip wall using the axisymmetric stream function, Fourier transformed in z."""
import mpmath as mp

def integrand(l):
    extra = int(3 * max(0, -mp.log10(l))) + 20
    with mp.extradps(extra):
        return _integrand(mp.mpf(l))

def _integrand(l):
    I0, I1 = mp.besseli(0, l), mp.besseli(1, l)
    K0, K1 = mp.besselk(0, l), mp.besselk(1, l)
    c = 2 / mp.pi / (8 * mp.pi)          # Stokeslet: psi = F/(8 pi mu) rho^2/r, F = mu = 1
    # unknowns A (rho I1), B (rho^2 I0);  R = 1
    r1, r2 = -c * K0, -c * (2 * K0 - l * K1)
    det = I1 * (2 * I0 + l * I1) - l * I0 ** 2
    A = (r1 * (2 * I0 + l * I1) - I0 * r2) / det
    B = (I1 * r2 - l * I0 * r1) / det
    return l * A + 2 * B

def k(dps=30):
    mp.mp.dps = dps
    w = mp.quad(integrand, [mp.mpf(10) ** (-dps - 10), mp.mpf(10) ** -6, 0.01, 1, 5, 20, 60, 3 * dps])
    return -6 * mp.pi * w

if __name__ == "__main__":
    print(k(30))
