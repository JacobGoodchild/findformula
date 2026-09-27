"""FCC flip-flop Grover walk. Pairs delta: (1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1).
G(x)  = E[cos(k.x)/(6 - 2s)],  G'(x) = E[cos(k.x)/(6 + 2s)],  s = c1c2 + c2c3 + c3c1.
k3 integrated exactly:  E_k3[cos(m k3)/(A - B cos k3)] = r^m / sqrt(A^2 - B^2), r = (A - sqrt(A^2-B^2))/B."""
import mpmath as mp, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 25

def lgf(x, sign):
    m1, m2, m3 = x
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2)
        A = 6 - sign * 2 * c1 * c2          # sign=+1: 6-2s ; sign=-1: 6+2s
        B = sign * 2 * (c1 + c2)
        root = mp.sqrt(A * A - B * B)
        if B == 0:
            e = (1 if m3 == 0 else 0) / A
        else:
            e = ((A - root) / B)**abs(m3) / root
        return mp.cos(m1 * k1) * mp.cos(m2 * k2) * e
    return mp.quad(f, [0, mp.pi], [0, mp.pi]) / mp.pi**2

if __name__ == "__main__":
    # needed separations: 0, delta, and t*d_q - s*d_p
    pts = [(0,0,0), (1,1,0), (2,2,0), (2,0,0), (1,1,2), (2,1,1), (0,2,0), (1,2,1), (0,0,2)]
    for sgn in (-1,):
        for x in pts:
            print(sgn, x, lgf(x, sgn))

def lgf_diff(x):
    """Delta(x) = G(0) - G(x) = E[(1 - cos k.x)/(6 - 2s)] (regular integrand)."""
    m1, m2, m3 = x
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2)
        A = 6 - 2 * c1 * c2; B = 2 * (c1 + c2)
        root = mp.sqrt(A * A - B * B)
        if root == 0:
            return mp.mpf(0)
        r = (A - root) / B if B != 0 else mp.mpf(0)
        return (1 - mp.cos(m1 * k1) * mp.cos(m2 * k2) * (r**abs(m3) if m3 else 1)) / root
    return mp.quad(f, [0, mp.pi], [0, mp.pi]) / mp.pi**2
