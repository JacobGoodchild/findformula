"""General engine: flip-flop Grover walk trapping on ANY 2D periodic graph.

Graph: nv vertices per cell, edges (u, v, (o1, o2)) meaning u in cell R -- v in cell R+o.
Grover coin at each vertex over its out-arcs, flip-flop shift.
Flat-band channels:
  +1  : <b|Pi|a> = (delta_ab - Y(a,b))/2,   Y  = d^* L^{-1} d  (L = D - A)
  -1  : <b|Pi|a> = (delta_ab - YQ(a,b))/2,  YQ = |d|^* Q^{-1} |d| (Q = D + A)
  lambda = mu +- i sqrt(1-mu^2) for each flat band mu of H = D^{-1/2} A D^{-1/2}:
        <b|Pi|a> = [F(o,o)/deg_o - lam F(h_b,o)/sqrt(deg_hb deg_o)
                    - conj(lam) F(o,h_a)/sqrt(deg_o deg_ha) + F(h_b,h_a)/sqrt(deg_hb deg_ha)] / (2 (1-mu^2))
Brillouin-zone averages E_k[N(z)/Det(z)] are done as: exact residue sum over z2 on the
unit circle, then tanh-sinh quadrature over k1.
"""
import sympy as sp, mpmath as mp

z1, z2, mu = sp.symbols('z1 z2 mu')

def bloch_mats(nv, edges):
    A = sp.zeros(nv, nv); deg = [0] * nv
    for (u, v, (o1, o2)) in edges:
        ph = z1**o1 * z2**o2
        A[u, v] += ph; A[v, u] += 1 / ph
        deg[u] += 1; deg[v] += 1
    D = sp.diag(*deg)
    return A, D, deg

class Averager:
    """E_k[N/Det] for Laurent polynomials N, Det in z1, z2."""
    def __init__(self, det):
        self.det = sp.expand(det)
        self.cache = {}
        self._prep(self.det)
    def _prep(self, det):
        # Det as polynomial in z2 with coefficients in z1:  z2^{-dmin} * P(z2)
        det = sp.expand(det)
        pw = [t.as_powers_dict().get(z2, 0) for t in sp.Add.make_args(det)]
        self.dmin = min(pw); dmax = max(pw)
        P = sp.Poly(sp.expand(det * z2**(-self.dmin)), z2)
        self.Pcoef = [sp.lambdify(z1, c, "mpmath") for c in P.all_coeffs()]   # highest first
    def inner(self, Ncoefs, nmin, zz1):
        """E_{k2}[N/Det] at fixed z1: N = z2^{nmin} * sum Ncoefs[j] z2^j (low->high)."""
        pc = [mp.mpc(f(zz1)) for f in self.Pcoef]
        # strip leading zeros
        # drop (near-)vanishing leading coefficients: the corresponding roots run off to infinity,
        # lie outside the unit circle and do not contribute; keeping them stalls polyroots
        big = max(abs(c) for c in pc)
        while abs(pc[0]) <= big * mp.mpf(10)**(-mp.mp.dps): pc = pc[1:]
        roots, err = mp.polyroots(pc, maxsteps=400, extraprec=4 * mp.mp.prec, error=True)
        if err > mp.mpf(10)**(-mp.mp.dps // 2):
            raise ValueError("polyroots did not converge: err=%s" % mp.nstr(err, 5))
        # Newton-polish the roots near/inside the unit circle (polyroots' accuracy is relative
        # to the largest root, which is poor when the leading coefficient is tiny)
        dpc0 = [c * (len(pc) - 1 - j) for j, c in enumerate(pc[:-1])]
        pol = []
        for r in roots:
            if abs(r) < 2:
                for _ in range(8):
                    r = r - mp.polyval(pc, r) / mp.polyval(dpc0, r)
            pol.append(r)
        roots = pol
        lead = pc[0]
        Nc = [mp.mpc(f(zz1)) for f in Ncoefs]
        # R(z)/z = z^{e} * PN(z) / PD(z),  e = nmin - dmin - 1
        e = nmin - self.dmin - 1
        PN = lambda z: mp.fsum(c * z**j for j, c in enumerate(Nc))
        res = mp.mpc(0)
        dpc = [c * (len(pc) - 1 - j) for j, c in enumerate(pc[:-1])]     # derivative coefficients
        tiny = [r for r in roots if abs(r) < mp.mpf('1e-3')]
        if e < 0 and tiny:
            # roots clustered at the pole z = 0: residues cancel catastrophically, so take the
            # combined residue of {0} U cluster as a contour integral on |z| = rho (trapezoid,
            # geometrically convergent), plus ordinary residues of the remaining inside roots
            rest = [r for r in roots if abs(r) >= mp.mpf('1e-3')]
            rc = max(abs(r) for r in tiny); ro = min(abs(r) for r in rest) if rest else mp.mpf(10)
            rho = mp.sqrt(rc * ro); q = mp.sqrt(rc / ro)
            M = int(mp.ceil(mp.mp.prec * mp.log(2) / -mp.log(q))) + 4
            g = lambda z: z**e * PN(z) / mp.polyval(pc, z)
            res += mp.fsum(g(rho * mp.expj(2 * mp.pi * j / M)) * rho * mp.expj(2 * mp.pi * j / M) for j in range(M)) / M
            for r in rest:
                if abs(r) < 1:
                    res += r**e * PN(r) / mp.polyval(dpc, r)
            return res
        for i, r in enumerate(roots):
            if abs(r) < 1:
                dP = mp.polyval(dpc, r)
                res += r**e * PN(r) / dP
        if e < 0:   # pole at z = 0 of order -e: residue = coeff of z^{-e-1} in PN/PD (PD(0) != 0)
            n = -e - 1
            pd = list(reversed(pc))            # low->high
            q = []                             # power series of PN/PD
            for k in range(n + 1):
                s = (Nc[k] if k < len(Nc) else 0) - mp.fsum(q[j] * (pd[k - j] if k - j < len(pd) else 0) for j in range(k))
                q.append(s / pd[0])
            res += q[n]
        return res
    def avg(self, N, breaks=None):
        N = sp.expand(N)
        if N == 0: return mp.mpf(0)
        pw = [t.as_powers_dict().get(z2, 0) for t in sp.Add.make_args(N)]
        nmin, nmax = min(pw), max(pw)
        Pn = sp.Poly(sp.expand(N * z2**(-nmin)), z2)
        co = list(reversed(Pn.all_coeffs()))
        Ncoefs = [sp.lambdify(z1, c, "mpmath") for c in co]
        def f(k1):
            with mp.extradps(2 * mp.mp.dps):
                return self.inner(Ncoefs, nmin, mp.exp(1j * k1))
        if breaks is None:
            pts = [j * mp.pi / 6 for j in range(-6, 7)]
        else:
            pts = sorted(set([-mp.pi, mp.pi] + [mp.mpf(b) for b in breaks]))
        v = mp.quad(lambda k: mp.re(f(k)), pts) / (2 * mp.pi)
        return v

def out_arcs(edges, x):
    arcs = []
    for (u, v, o) in edges:
        if u == x: arcs.append((v, o))
        if v == x: arcs.append((u, (-o[0], -o[1])))
    return arcs

def phase(p, q):   # Fourier kernel factor for G(p at cell op, q at cell oq): z^(op - oq)
    return z1**(p[0] - q[0]) * z2**(p[1] - q[1])

def transfer_matrix(nv, edges, origin, signless):
    A, D, deg = bloch_mats(nv, edges)
    M = D + A if signless else D - A
    det = sp.expand(M.det()); adj = M.adjugate()
    av = Averager(det)
    arcs = out_arcs(edges, origin)
    s = 1 if signless else -1
    O = (origin, (0, 0))
    def G(p, q): return adj[p[0], q[0]] * phase(p[1], q[1])
    n = len(arcs); Y = [[None] * n for _ in range(n)]
    for i, hi in enumerate(arcs):
        for j, hj in enumerate(arcs):
            if j < i: Y[i][j] = Y[j][i]; continue
            N = G(O, O) + s * G(O, hj) + s * G(hi, O) + G(hi, hj)
            Y[i][j] = av.avg(N)
    return Y, arcs, deg

def flat_bands(nv, edges):
    """flat eigenvalues mu of H = D^{-1/2} A D^{-1/2} (constant in k), with multiplicities."""
    A, D, deg = bloch_mats(nv, edges)
    Dm = sp.diag(*[1 / sp.sqrt(d) for d in deg])
    H = Dm * A * Dm
    p = sp.expand((mu * sp.eye(nv) - H).det())
    num = sp.numer(sp.together(p))
    # factors independent of z1, z2 give flat bands
    fl = []
    for fac, m in sp.factor_list(num)[1]:
        if not fac.has(z1) and not fac.has(z2) and fac.has(mu):
            for r in sp.solve(fac, mu):
                fl.append((r, m))
    return fl, H, deg

def flat_projector_channel(nv, edges, origin, mu0, H, deg, start):
    """contribution of lambda = mu0 +- i sqrt(1-mu0^2) channels (simple flat band)."""
    X = mu * sp.eye(nv) - H
    p = sp.expand(X.det())
    dp = sp.diff(p, mu).subs(mu, mu0)
    adjX = X.adjugate().subs(mu, mu0)
    denom = sp.expand(sp.numer(sp.together(dp)))
    scale = sp.simplify(sp.together(dp) / denom)          # dp = denom * scale (scale indep of mu)
    av = Averager(denom)
    arcs = out_arcs(edges, origin)
    O = (origin, (0, 0))
    def F(p_, q_):
        N = sp.expand(sp.together(adjX[p_[0], q_[0]] / scale) * phase(p_[1], q_[1]))
        N = sp.expand(sp.numer(sp.together(N)) / sp.denom(sp.together(N)))
        return av.avg(N)
    ha = arcs[start]
    Foo = F(O, O); Foa = F(O, ha)
    th2 = 1 - mu0**2
    out = mp.mpf(0)
    for lam in (mp.mpc(float(mu0), mp.sqrt(th2)), mp.mpc(float(mu0), -mp.sqrt(th2))):
        lam = mp.mpc(mp.mpf(sp.N(mu0, 60)), (1 if lam.imag > 0 else -1) * mp.sqrt(mp.mpf(sp.N(th2, 60))))
        s = mp.mpf(0)
        for hb in arcs:
            val = (Foo / deg[origin] - lam * F(hb, O) / mp.sqrt(deg[hb[0]] * deg[origin])
                   - mp.conj(lam) * Foa / mp.sqrt(deg[origin] * deg[ha[0]]) + F(hb, ha) / mp.sqrt(deg[hb[0]] * deg[ha[0]])) / (2 * mp.mpf(sp.N(th2, 60)))
            s += abs(val)**2
        out += s
    return out

def trapping(nv, edges, origin=0, start=0, bipartite=None, verbose=True):
    Y, arcs, deg = transfer_matrix(nv, edges, origin, False)
    n = len(arcs)
    pp = mp.fsum(((1 if j == start else 0) - Y[start][j])**2 for j in range(n)) / 4
    YQ, _, _ = transfer_matrix(nv, edges, origin, True) if not bipartite else (Y, None, None)
    pm = mp.fsum(((1 if j == start else 0) - YQ[start][j])**2 for j in range(n)) / 4
    fl, H, _ = flat_bands(nv, edges)
    extra = {}
    for (m0, mult) in fl:
        if mult != 1: raise NotImplementedError("degenerate flat band")
        if abs(sp.N(m0)) >= 1: continue
        extra[str(m0)] = flat_projector_channel(nv, edges, origin, m0, H, deg, start)
    tot = pp + pm + mp.fsum(extra.values())
    if verbose:
        print(" rowsum Y:", mp.nstr(mp.fsum(Y[start]), 12), " deg:", deg)
        print(" +1:", pp, "\n -1:", pm, "\n flat:", extra, "\n TOTAL:", tot)
    return tot, pp, pm, extra, Y, YQ

if __name__ == "__main__":
    mp.mp.dps = 20
    tri = [(0, 0, (1, 0)), (0, 0, (0, 1)), (0, 0, (1, -1))]
    print("triangular (expect 0.26795570253708055129):"); trapping(1, tri)
    kag = [(0, 1, (0, 0)), (0, 2, (0, 0)), (1, 2, (0, 0)), (1, 0, (1, 0)), (2, 0, (0, 1)), (1, 2, (1, -1))]
    print("kagome (expect 0.19146243000085471390):"); trapping(3, kag)
