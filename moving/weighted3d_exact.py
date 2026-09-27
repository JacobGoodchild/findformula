"""Exact A(p) = (1/2) E[min_k Z_k^2/p_k] for d = 3 via spherical triangles:
 int_T x_i^2 dsigma = Omega/3 - (1/3) sum_edges m_{e,i} (a_e + b_e)_i tan(Theta_e/2)   (m_e outward edge normal).
 E[Z_i^2 1_C] = (3/(4 pi)) int_{C cap S^2} x_i^2 dsigma ; each cone C_i = 8 congruent octant pieces."""
import mpmath as mp
mp.mp.dps = 50
def tri_integral(V, i):
    """V: 3 unit vertices (counterclockwise seen from outside). returns (Omega, int x_i^2)."""
    a, b, c = [mp.matrix(v) for v in V]
    # solid angle (Van Oosterom-Strackee)
    num = a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])
    den = 1 + (a.T*b)[0] + (b.T*c)[0] + (c.T*a)[0]
    Om = 2 * mp.atan2(abs(num), den)
    s = Om / 3
    centroid = (a + b + c) / 3
    for (u, w) in [(a, b), (b, c), (c, a)]:
        cr = mp.matrix([u[1]*w[2]-u[2]*w[1], u[2]*w[0]-u[0]*w[2], u[0]*w[1]-u[1]*w[0]])
        m = cr / mp.norm(cr)
        if (m.T * centroid)[0] > 0: m = -m        # outward
        Th = mp.acos((u.T * w)[0])
        s -= m[i] * (u[i] + w[i]) * mp.tan(Th / 2) / 3
    return Om, s
def A_exact(p):
    p = [mp.mpf(x) for x in p]
    tot = 0; pis = []
    for i in range(3):
        j, k = [x for x in range(3) if x != i]
        # octant piece of cone C_i: 0 < Z_i < a_j Z_j, Z_i < a_k Z_k, a = sqrt(p_i/p_*)
        ej = [0, 0, 0]; ej[j] = 1; ek = [0, 0, 0]; ek[k] = 1
        v = [0, 0, 0]; v[i] = 1; v[j] = mp.sqrt(p[j] / p[i]); v[k] = mp.sqrt(p[k] / p[i])
        nv = mp.sqrt(sum(t * t for t in v)); v = [t / nv for t in v]
        Om, s = tri_integral([ej, ek, v], i)
        pis.append(8 * Om / (4 * mp.pi))
        tot += (1 / p[i]) * 8 * (3 / (4 * mp.pi)) * s
    return tot / 2, pis
p = (mp.sqrt(2)/10, mp.sqrt(7)/10, 1 - mp.sqrt(2)/10 - mp.sqrt(7)/10)
A, pis = A_exact(p)
num = mp.quad(lambda u: mp.fprod(mp.erfc(mp.sqrt(pk * u)) for pk in p), [0, 0.5, 1, 2, 5, 10, 20, 40, 80, mp.inf])
print("A exact", A, "\nA num  ", num, "\npis", pis, sum(pis))
p = (mp.mpf(1)/3,)*3; A, pis = A_exact(p); print("isotropic A/3 =", A/3, " vs 1/2 - (3-sqrt3)/pi =", mp.mpf(1)/2 - (3 - mp.sqrt(3))/mp.pi)
