"""Exact contribution of the kagome flat band (random-walk eigenvalue -1/2) to flip-flop
Grover-walk trapping. Szegedy eigenvalues lambda = e^{+-2 pi i/3};
<b|Pi_lambda|a> = (1/6)[F(o,o) - lam F(h_b,o) - conj(lam) F(o,h_a) + F(h_b,h_a)],
F = real-space projector onto the kagome flat band (adjacency eigenvalue -2)."""
import sympy as sp, mpmath as mp
from engine import bloch, laurent_terms, match_triangular, expect, z1, z2
mp.mp.dps = 35
kag = [(0, 1, (0, 0)), (0, 2, (0, 0)), (1, 2, (0, 0)), (1, 0, (1, 0)), (2, 0, (0, 1)), (1, 2, (1, -1))]
Q = bloch(3, kag, True)                 # 4 + A_K
M = Q - 2 * sp.eye(3)                   # A_K + 2  (rank 2 generically)
r0, r1 = M.row(0), M.row(1)
conjz = {z1: 1/z1, z2: 1/z2}
v = sp.Matrix([r0[1]*r1[2] - r0[2]*r1[1], r0[2]*r1[0] - r0[0]*r1[2], r0[0]*r1[1] - r0[1]*r1[0]])
v = sp.Matrix([sp.factor(x) for x in v])
g = sp.gcd(sp.gcd(sp.numer(sp.together(v[0])), sp.numer(sp.together(v[1]))), sp.numer(sp.together(v[2])))
print("common factor:", g)
v = sp.Matrix([sp.expand(sp.cancel(x / g)) for x in v])
print("v =", list(v))
print("check M v = 0:", [sp.simplify(sp.expand(x)) for x in (M * v)])
vdag = sp.Matrix([[sp.expand(x.subs(conjz, simultaneous=True)) for x in v]])
norm = sp.expand((vdag * v)[0])
print("|v|^2 =", sp.factor(norm))
kappa, alpha = match_triangular(norm)
print("kappa, alpha:", kappa, alpha)
# F(p at op, q at oq) = E[ v_p conj(v_q) / |v|^2 * z^(op - oq) ]
def F(p, op, q, oq):
    num = sp.expand(v[p] * vdag[q] * z1**(op[0] - oq[0]) * z2**(op[1] - oq[1]))
    return expect(num, kappa, alpha, cache)
cache = {}
arcs = []
for (u, w, o) in kag:
    if u == 0: arcs.append((w, o))
    if w == 0: arcs.append((u, (-o[0], -o[1])))
lam = mp.exp(2j * mp.pi / 3)
o0 = (0, (0, 0))
a = arcs[0]
total = 0
for L in (lam, mp.conj(lam)):
    s = 0
    for b in arcs:
        val = (F(0, (0,0), 0, (0,0)) - L * F(b[0], b[1], 0, (0,0)) - mp.conj(L) * F(0, (0,0), a[0], a[1]) + F(b[0], b[1], a[0], a[1])) / 6
        s += abs(val)**2
    print("channel", mp.nstr(L, 5), s)
    total += s
print("flat-band total:", total)
print("F(o,o) =", F(0,(0,0),0,(0,0)))
