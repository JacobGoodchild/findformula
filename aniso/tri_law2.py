import mpmath as mp, sympy as sp
exec(open('tri_law.py').read().split("q1, q2, q3 = sp.symbols")[0])
for q, rel in data.items():
    S = q[0]*q[1] + q[1]*q[2] + q[2]*q[0]
    c = {"c1^2": S/q[0]**2, "c2^2": S/q[1]**2, "c3^2": S/q[2]**2, "c1c2": S/(q[0]*q[1]), "c1c3": S/(q[0]*q[2]), "c2c3": S/(q[1]*q[2])}
    print("q", q, {k: v for k, v in c.items()})
    for b in [(1, 1), (-1, 1), (1, 2), (-1, 2)]:
        print("    ", b, "const, t1, t2, cot1/pi coeffs:", rel[b])
