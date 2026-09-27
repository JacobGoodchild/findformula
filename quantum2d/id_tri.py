import mpmath as mp
from lgf2d import a, b
mp.mp.dps = 45
pi = mp.pi; s3 = mp.sqrt(3)
A = lambda k: 3 - mp.cos(k); C = lambda k: 1 + mp.exp(-1j * k)
Ap = lambda k: 3 + mp.cos(k); Cp = lambda k: -(1 + mp.exp(-1j * k))
vals = {"a(1,1)": a((1, 1), A, C), "a(2,0)": a((2, 0), A, C),
        "G'(0)": b((0, 0), Ap, Cp), "G'(1,0)": b((1, 0), Ap, Cp), "G'(1,1)": b((1, 1), Ap, Cp), "G'(2,0)": b((2, 0), Ap, Cp)}
for k, v in vals.items(): print(k, v)
g0 = vals["G'(0)"]
cands = {"G(1/3)^3/pi^2": mp.gamma(mp.mpf(1)/3)**3/pi**2, "G(1/3)^6/pi^4": mp.gamma(mp.mpf(1)/3)**6/pi**4,
         "G(1/4)^2/pi^(3/2)": mp.gamma(mp.mpf(1)/4)**2/pi**1.5, "G(1/4)^4/pi^3": mp.gamma(mp.mpf(1)/4)**4/pi**3,
         "G(1/6)G(1/3)/pi^(3/2)": mp.gamma(mp.mpf(1)/6)*mp.gamma(mp.mpf(1)/3)/pi**1.5,
         "G(1/8)G(3/8)/pi^(3/2)": mp.gamma(mp.mpf(1)/8)*mp.gamma(mp.mpf(3)/8)/pi**1.5,
         "G(1/24)..G(11/24)/pi^3": mp.gamma(mp.mpf(1)/24)*mp.gamma(mp.mpf(5)/24)*mp.gamma(mp.mpf(7)/24)*mp.gamma(mp.mpf(11)/24)/pi**3,
         "1/pi": 1/pi, "1": mp.mpf(1)}
algs = {"1": 1, "s2": mp.sqrt(2), "s3": s3, "s6": mp.sqrt(6), "2^(1/3)": mp.cbrt(2), "4^(1/3)": mp.cbrt(4), "3^(1/4)": mp.root(3,4), "s3*2^(1/3)": s3*mp.cbrt(2), "s3*4^(1/3)": s3*mp.cbrt(4), "27^(1/4)": mp.root(27,4)}
for cn, c in cands.items():
    for an, al in algs.items():
        r = mp.pslq([g0, c * al], maxcoeff=5000, maxsteps=10**5)
        if r: print("G'(0) HIT", cn, an, r)
for k in ["a(1,1)", "a(2,0)"]:
    print(k, mp.pslq([vals[k], 1, s3/pi, 1/pi], maxcoeff=10**4, maxsteps=10**5))
X = mp.cbrt(4) * mp.gamma(mp.mpf(1)/3)**3 / (8 * pi**2)
print("G'(0)-X:", g0 - X)
for k in ["G'(1,0)", "G'(1,1)", "G'(2,0)"]:
    print(k, mp.pslq([vals[k], 1, X, 1/(pi**2*X), s3/pi, 1/pi, s3/(pi**3*X), s3*X/pi], maxcoeff=10**4, maxsteps=10**6))
print("--- wider basis")
names = ["v","1","X","1/(pi^2X)","s3/(pi^2X)","s3X","s3/pi","1/pi","s3","1/(pi X)","s3/(pi X)", "X/pi", "s3 X/pi"]
for k in ["G'(1,1)", "G'(2,0)"]:
    B = [vals[k], 1, X, 1/(pi**2*X), s3/(pi**2*X), s3*X, s3/pi, 1/pi, s3, 1/(pi*X), s3/(pi*X), X/pi, s3*X/pi]
    r = mp.pslq(B, maxcoeff=3000, maxsteps=10**6)
    print(k, [(n, c) for n, c in zip(names, r) if c] if r else None)
