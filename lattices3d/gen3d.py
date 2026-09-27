"""3D version of the flip-flop Grover trapping engine.
E_k[N/Det] : exact residue sum over z3, tanh-sinh quadrature over (k1, k2)."""
import sympy as sp, mpmath as mp, itertools
import numpy as np

z1, z2, z3 = sp.symbols('z1 z2 z3')

def geometric3(avecs, basis, bond2, rng=2):
    A = [np.array(a, float) for a in avecs]; B = [np.array(b, float) for b in basis]
    edges = set()
    for u, ru in enumerate(B):
        for v, rv in enumerate(B):
            for o in itertools.product(range(-rng, rng + 1), repeat=3):
                d = rv + sum(oi * ai for oi, ai in zip(o, A)) - ru
                if abs(d @ d - bond2) < 1e-9 and not (u == v and o == (0, 0, 0)):
                    key = (u, v, o); rk = (v, u, tuple(-x for x in o))
                    if rk not in edges: edges.add(key)
    return len(B), sorted(edges)

def bloch(nv, edges):
    A = sp.zeros(nv, nv); deg = [0] * nv
    for (u, v, o) in edges:
        ph = z1**o[0] * z2**o[1] * z3**o[2]
        A[u, v] += ph; A[v, u] += 1 / ph; deg[u] += 1; deg[v] += 1
    return A, sp.diag(*deg), deg

def zpow(t, z):
    return t.as_powers_dict().get(z, 0)

class Avg3:
    def __init__(self, det):
        det = sp.expand(det)
        self.dmin = min(zpow(t, z3) for t in sp.Add.make_args(det))
        P = sp.Poly(sp.expand(det * z3**(-self.dmin)), z3)
        self.Pc = [sp.lambdify((z1, z2), c, "mpmath") for c in P.all_coeffs()]
    def inner(self, Nc, nmin, a, b):
        pc = [mp.mpc(f(a, b)) for f in self.Pc]
        while abs(pc[0]) == 0: pc = pc[1:]
        roots = mp.polyroots(pc, maxsteps=300, extraprec=4 * mp.mp.prec)
        lead = pc[0]
        Ncv = [mp.mpc(f(a, b)) for f in Nc]
        e = nmin - self.dmin - 1
        PN = lambda z: mp.fsum(c * z**j for j, c in enumerate(Ncv))
        res = mp.mpc(0)
        for i, r in enumerate(roots):
            if abs(r) < 1:
                dP = lead
                for j, s in enumerate(roots):
                    if j != i: dP *= (r - s)
                res += r**e * PN(r) / dP
        if e < 0:
            n = -e - 1; pd = list(reversed(pc)); q = []
            for k in range(n + 1):
                s = (Ncv[k] if k < len(Ncv) else 0) - mp.fsum(q[j] * (pd[k - j] if k - j < len(pd) else 0) for j in range(k))
                q.append(s / pd[0])
            res += q[n]
        return res
    def avg(self, N):
        N = sp.expand(N)
        if N == 0: return mp.mpf(0)
        nmin = min(zpow(t, z3) for t in sp.Add.make_args(N))
        P = sp.Poly(sp.expand(N * z3**(-nmin)), z3)
        Nc = [sp.lambdify((z1, z2), c, "mpmath") for c in reversed(P.all_coeffs())]
        def f(k1, k2):
            with mp.extradps(mp.mp.dps):
                return mp.re(self.inner(Nc, nmin, mp.exp(1j * k1), mp.exp(1j * k2)))
        pts = [-mp.pi, 0, mp.pi]
        return mp.quad(f, pts, pts) / (4 * mp.pi**2)

def ph(p, q):
    return z1**(p[0] - q[0]) * z2**(p[1] - q[1]) * z3**(p[2] - q[2])

def out_arcs(edges, x):
    arcs = []
    for (u, v, o) in edges:
        if u == x: arcs.append((v, o))
        if v == x: arcs.append((u, tuple(-t for t in o)))
    return arcs

def transfer_row(nv, edges, origin, start, signless):
    A, D, deg = bloch(nv, edges)
    M = D + A if signless else D - A
    av = Avg3(sp.expand(M.det())); adj = M.adjugate()
    arcs = out_arcs(edges, origin); s = 1 if signless else -1
    O = (origin, (0, 0, 0)); G = lambda p, q: adj[p[0], q[0]] * ph(p[1], q[1])
    hi = arcs[start]
    row = []
    for hj in arcs:
        row.append(av.avg(G(O, O) + s * G(O, hj) + s * G(hi, O) + G(hi, hj)))
        print("   ", "YQ" if signless else "Y", hj, mp.nstr(row[-1], 20), flush=True)
    return row, arcs, deg
