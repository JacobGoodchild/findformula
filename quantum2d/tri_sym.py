import sympy as sp
X, pi = sp.symbols('X pi', positive=True)
s3 = sp.sqrt(3)
Gp = {0: X, 1: (1 - 3*X)/3, "r3": (9*X - 2 - s3/(pi*X))/3, 2: (3*X - 4 + 2*s3/(pi*X))/3}
aa = {0: 0, 1: sp.Rational(1,3), "r3": 2*s3/pi - sp.Rational(2,3), 2: sp.Rational(8,3) - 4*s3/pi}
def cls(y):
    n2 = y[0]**2 + y[1]**2 + y[0]*y[1]
    return {0: 0, 1: 1, 3: "r3", 4: 2}[n2]
deltas = [(1,0),(0,1),(1,-1)]
states = [(s,p) for p in range(3) for s in (1,-1)]
pp = 0; pm = 0
i = 0; (s, p) = states[0]
for j,(t,q) in enumerate(states):
    y = (t*deltas[q][0]-s*deltas[p][0], t*deltas[q][1]-s*deltas[p][1]); c = cls(y)
    mp_ = sp.Rational(1,2)*(1 if j==0 else 0) - (2*aa[1]-aa[c])/4
    mm_ = sp.Rational(1,2)*(1 if j==0 else 0) - (Gp[0]+2*Gp[1]+Gp[c])/4
    pp += mp_**2; pm += mm_**2
print("+1 channel:", sp.simplify(sp.expand(pp)))
print("-1 channel:", sp.simplify(sp.expand(pm)))
tot = sp.expand(pp + pm)
print("total:", tot)
Xv = sp.Rational(1)*sp.cbrt(4)*sp.gamma(sp.Rational(1,3))**3/(8*sp.pi**2)
print(sp.N(tot.subs({X: Xv, pi: sp.pi}), 30))
