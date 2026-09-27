"""BCC lattice Green function via series:
G(x) = E[cos(k.x)/(4(1 - c1 c2 c3))] = (1/4) sum_n prod_i E[cos(x_i k) cos(k)^n],
E[cos(m k) cos^n k] = C(n,(n-m)/2)/2^n  (n>=|m|, n-m even)."""
import mpmath as mp
mp.mp.dps = 40

def G(x):
    m = [abs(v) for v in x]
    par = m[0] % 2
    assert all(v % 2 == par for v in m)
    j0 = (max(m) - par) // 2          # first j with n >= max(m)
    def term(j):
        n = 2 * (int(j) + j0) + par
        p = mp.mpf(1)
        for mi in m:
            p *= mp.binomial(n, (n - mi) // 2) / mp.mpf(2)**n
        return p
    return mp.nsum(term, [0, mp.inf], method="levin") / 4

if __name__ == "__main__":
    Wb = mp.gamma(mp.mpf(1) / 4)**4 / (4 * mp.pi**3)
    print("4*G(0) =", 4 * G((0, 0, 0)), " Watson BCC =", Wb)
