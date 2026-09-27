"""High-precision lattice Green-function differences for 2D Bravais lattices whose
dispersion is  D(k) = A(k1) - Re(C(k1) e^{i k2}) (only e^{+-i k2} harmonics).
a(x) = E[(1 - cos k.x)/D(k)]  (potential kernel; finite even when E[1/D] diverges)
b(x) = E[cos(k.x)/D(k)]       (only for gapped D)"""
import mpmath as mp

def _k2avg(A, C, m):
    """E_{k2}[ e^{i m k2} / (A - Re(C e^{ik2})) ]"""
    R = abs(C)
    root = mp.sqrt(A * A - R * R)
    if R == 0:
        return (1 if m == 0 else 0) / A
    r = (A - root) / R
    return r**abs(m) * mp.exp(-1j * m * mp.arg(C)) / root

def a(x, A, C):
    x1, x2 = x
    def f(k1):
      with mp.extradps(mp.mp.dps):
        Ak, Ck = A(k1), C(k1)
        if abs(Ak * Ak - abs(Ck)**2) < mp.mpf(10)**(-2 * mp.mp.dps):
            return mp.mpf(0)
        v = _k2avg(Ak, Ck, 0) - (mp.exp(1j * x1 * k1) * _k2avg(Ak, Ck, x2)).real
        return +mp.re(v)
    return mp.quad(f, [-mp.pi, 0, mp.pi]) / (2 * mp.pi)

def b(x, A, C):
    x1, x2 = x
    f = lambda k1: mp.re(mp.exp(1j * x1 * k1) * _k2avg(A(k1), C(k1), x2))
    return mp.quad(f, [-mp.pi, 0, mp.pi]) / (2 * mp.pi)

if __name__ == "__main__":
    mp.mp.dps = 30
    # square lattice D = 2 - cos k1 - cos k2
    A = lambda k1: 2 - mp.cos(k1); C = lambda k1: mp.mpf(1)
    print("square a(1,0) =", a((1,0), A, C), "(1/2)")
    print("square a(1,1) =", a((1,1), A, C), "vs 2/pi =", 2/mp.pi)
    print("square a(2,0) =", a((2,0), A, C), "vs 2-4/pi =", 2 - 4/mp.pi)
