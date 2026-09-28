import sympy as sp, mpmath as mp
from ustdeg import degree_distribution
mp.mp.dps = 20
out = {}
# pyrochlore (Formula 11): tetra A arcs 0-2, tetra B arcs 3-5 (arc 3+i collinear with arc i)
W, X = sp.symbols('W X'); g = (1 + 8*W + 3*X)/48; d = (15 - 8*W - 3*X)/96
Y = sp.Matrix(6, 6, lambda i, j: (sp.Rational(1,3) if i == j else sp.Rational(1,6)) if (i < 3) == (j < 3) else (g if i % 3 == j % 3 else d))
Wv = 3*mp.gamma(mp.mpf(1)/3)**6/(2**(mp.mpf(14)/3)*mp.pi**4)
out['pyrochlore'] = (degree_distribution(Y), {W: Wv, X: 1/(mp.pi**2*Wv)})
# diamond: Y = 1/2 diag, 1/6 off
Y = sp.Matrix(4, 4, lambda i, j: sp.Rational(1,2) if i == j else sp.Rational(1,6)); out['diamond'] = (degree_distribution(Y), {})
# king (Formula 29), arcs [+x,-x,+y,-y,(1,1),(-1,-1),(1,-1),(-1,1)], p = 1/pi, L = ln(2+sqrt3)/(pi sqrt3)
p, L = sp.symbols('p L')
vec = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)]
Rk = {(1,0): L, (1,1): sp.Rational(1,2) - L, (2,0): 4*p - 4*L, (2,1): 6*L - sp.Rational(1,2) - 2*p, (2,2): 5 - 4*p - 14*L}
def Rv(x):
    k = tuple(sorted(map(abs, x), reverse=True)); return 0 if k == (0,0) else Rk[(k[0], k[1])]
Y = sp.Matrix(8, 8, lambda i, j: (Rv(vec[i]) + Rv(vec[j]) - Rv((vec[i][0]-vec[j][0], vec[i][1]-vec[j][1])))/2)
out['king'] = (degree_distribution(Y), {p: 1/mp.pi, L: mp.log(2+mp.sqrt(3))/(mp.pi*mp.sqrt(3))})
# snub square / Shastry-Sutherland (Formula 25): arcs [a, b, d, a', b'] ; s3 = sqrt3
s3 = sp.Symbol('s')
Ya = [[sp.Rational(1,2) - s3/18, 1 - s3/2, s3/9, 5*s3/9 - 1 + p/2, sp.Rational(1,2) - s3/9 - p/2],
      [None, sp.Rational(1,2) - s3/18, s3/9, sp.Rational(1,2) - s3/9 - p/2, 5*s3/9 - 1 + p/2],
      [None, None, 2*s3/9, sp.Rational(1,2) - 2*s3/9, sp.Rational(1,2) - 2*s3/9],
      [None, None, None, sp.Rational(1,2) - s3/18, (3 - s3)/6],
      [None, None, None, None, sp.Rational(1,2) - s3/18]]
Y = sp.Matrix(5, 5, lambda i, j: Ya[min(i,j)][max(i,j)])
print('SS row sums', [sp.simplify(sum(Y.row(i))) for i in range(5)])
out['snub_square'] = ([sp.factor(sp.expand(v.subs(s3, sp.sqrt(3)))) for v in degree_distribution(Y)], {p: 1/mp.pi})
for n, (P, sub) in out.items():
    subs = {k: sp.Float(str(v), 25) for k, v in sub.items()}
    vals = [sp.N(v.subs(subs), 18) for v in P]
    print(n, 'sum', sp.N(sum(vals), 12), 'mean', sp.N(sum(k*v for k, v in enumerate(vals)), 12))
    for k, v in enumerate(P):
        if v != 0: print('   deg', k, ':', v, ' =', vals[k])
