import sympy as sp, mpmath as mp
K, E, pi = sp.symbols('K E pi', positive=True)
s2, s3 = sp.sqrt(2), sp.sqrt(3)
Y = [sp.Rational(1,2), sp.Rational(1,4), (1 + 6*s3/pi)/36, (4 - 3*s3/pi)/18]
YQ = [sp.Rational(1,2), (5 - 5*s2*K/pi)/12, (24*s2*E/pi - 4*s2*K/pi - 11)/12, (4 + 3*s2*K/pi - 8*s2*E/pi)/4]
pp = sp.expand(sum(((1 if j == 0 else 0) - Y[j])**2 for j in range(4)) / 4)
pm = sp.expand(sum(((1 if j == 0 else 0) - YQ[j])**2 for j in range(4)) / 4)
flat = (19 - 42*s3/pi + 108/pi**2) / 648
tot = sp.expand(pp + pm + 2*flat)
print("+1:", sp.simplify(pp)); print("-1:", sp.collect(pm, [K, E])); print("total:", sp.collect(tot, [K, E]))
mp.mp.dps = 30
m = mp.mpf(5)/32
sub = {K: sp.Float(str(mp.ellipk(m)), 30), E: sp.Float(str(mp.ellipe(m)), 30), pi: sp.pi}
print("numeric total:", sp.N(tot.subs(sub), 25), " -1:", sp.N(pm.subs(sub), 25))
