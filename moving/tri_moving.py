"""Moving-shift Grover walk on the triangular lattice: exact flat-band projector averages.
k uniform <-> x = tan(k1/2), y = tan(k2/2) independent standard Cauchy.  tan((k1-k2)/2) = (x-y)/(1+xy)."""
import sympy as sp, mpmath as mp, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 25
x, y = sp.symbols('x y')
t = [x, y, (x - y) / (1 + x * y)]            # tan(theta_j / 2) for directions e1, e2, e1-e2
u = [1 / x, 1 / y, (1 + x * y) / (x - y)]    # cot(theta_j / 2)
arcs = [(s, p) for p in range(3) for s in (1, -1)]
def entry(lam, b, a):
    if lam == 1:
        vb = 1 + sp.I * b[0] * t[b[1]]; va = 1 + sp.I * a[0] * t[a[1]]; den = 2 * (3 + sum(z**2 for z in t))
    else:
        vb = 1 - sp.I * b[0] * u[b[1]]; va = 1 - sp.I * a[0] * u[a[1]]; den = 2 * (3 + sum(z**2 for z in u))
    return sp.cancel(sp.together(vb * sp.conjugate(va).subs({sp.conjugate(x): x, sp.conjugate(y): y}) / den))
def avg(expr):
    """E over Cauchy x, y:  (1/pi^2) int int expr dx dy/((1+x^2)(1+y^2)); inner y-integral by residues."""
    R = sp.cancel(sp.together(expr / (1 + y**2)))
    N, D = sp.fraction(R)
    Np = sp.Poly(sp.expand(N), y); Dp = sp.Poly(sp.expand(D), y)
    Nc = [sp.lambdify(x, c, "mpmath") for c in Np.all_coeffs()]
    Dc = [sp.lambdify(x, c, "mpmath") for c in Dp.all_coeffs()]
    def inner(xv):
        with mp.extradps(mp.mp.dps):
            dc = [mp.mpc(f(xv)) for f in Dc]
            while abs(dc[0]) < mp.mpf(10)**(-3 * mp.mp.dps): dc = dc[1:]
            nc = [mp.mpc(f(xv)) for f in Nc]
            roots = mp.polyroots(dc, maxsteps=400, extraprec=6 * mp.mp.prec)
            tot = mp.mpc(0)
            for i, r in enumerate(roots):
                if mp.im(r) > 0:
                    dP = dc[0]
                    for j, s in enumerate(roots):
                        if j != i: dP *= (r - s)
                    tot += mp.polyval(nc, r) / dP
            return 2j * mp.pi * tot / mp.pi          # (1/pi) int R dy
    f = lambda xv: mp.re(inner(xv)) / (mp.pi * (1 + xv**2))
    fi = lambda xv: mp.im(inner(xv)) / (mp.pi * (1 + xv**2))
    return mp.quad(f, [-mp.inf, -1, 0, 1, mp.inf]) + 1j * mp.quad(fi, [-mp.inf, -1, 0, 1, mp.inf])
if __name__ == "__main__":
    a = arcs[0]
    tot = 0; vals = {}
    for lam in (1, -1):
        s = 0
        for b in arcs:
            v = avg(entry(lam, b, a)); vals[(lam, b)] = v; s += abs(v)**2
            print(lam, b, mp.nstr(v, 20), flush=True)
        print("channel", lam, s, flush=True); tot += s
    print("TOTAL", tot)
    import pickle; pickle.dump(vals, open(f"tri_moving_{mp.mp.dps}.pkl", "wb"))
