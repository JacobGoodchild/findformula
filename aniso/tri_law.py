import mpmath as mp, sympy as sp, itertools
from tri_w import channels
mp.mp.dps = 30
arcs = [(s, p) for p in range(3) for s in (1, -1)]
Q = [(sp.Rational(2,5), sp.Rational(7,20), sp.Rational(1,4)), (sp.Rational(9,20), sp.Rational(7,20), sp.Rational(1,5)),
     (sp.Rational(11,25), sp.Rational(8,25), sp.Rational(6,25)), (sp.Rational(1,2), sp.Rational(3,10), sp.Rational(1,5)),
     (sp.Rational(3,10), sp.Rational(3,10), sp.Rational(2,5)), (sp.Rational(3,5), sp.Rational(1,4), sp.Rational(3,20))]
data = {}
for q in Q:
    qm = [mp.mpf(sp.N(x, 40)) for x in q]
    pp, pm, Yr, Qr = channels(qm)
    lam = mp.findroot(lambda l: sum(mp.atan(l * x) for x in qm) - mp.pi / 2, 1)
    t = [mp.atan(lam * x) / mp.pi for x in qm]; c1 = 1 / (lam * qm[0])
    rel = {}
    for b, y in zip(arcs, Yr):
        T = mp.sqrt(qm[b[1]] / qm[0]) * y
        r = mp.pslq([T, 1, t[0], t[1], c1 / mp.pi], maxcoeff=10**6, maxsteps=10**7)
        rel[b] = [-sp.Rational(x, r[0]) for x in r[1:]]      # T = r1 + r2 t1 + r3 t2 + r4 c1/pi
    data[q] = rel
q1, q2, q3 = sp.symbols('q1 q2 q3')
monos = [sp.Integer(1), q1, q2, q1**2, q1*q2, q2**2]   # (q3 = 1 - q1 - q2) ; coefficients times q1 (and q1^2) allowed
for b in arcs[1:]:
    print("arc", b)
    for idx, name in enumerate(["const", "t1", "t2", "cot1/pi"]):
        for scale_name, scale in [("", 1), ("*q1", q1), ("*q1^2", q1**2)]:
            cs = sp.symbols('k0:%d' % len(monos))
            eqs = []
            for q, rel in data.items():
                sub = {q1: q[0], q2: q[1], q3: q[2]}
                eqs.append(sp.Eq(sum(c * m.subs(sub) for c, m in zip(cs, monos)), rel[b][idx] * scale.subs(sub) if scale != 1 else rel[b][idx]))
            sol = sp.solve(eqs, cs, dict=True)
            if sol:
                expr = sp.factor(sum(sol[0][c] * m for c, m in zip(cs, monos)).subs(q2, q2) / scale)
                print("   ", name, "=", expr, "(fit with", scale_name or "x1", ")")
                break
        else:
            print("   ", name, ": no low-degree fit")
