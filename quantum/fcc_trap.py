import mpmath as mp, pickle
from fcc import lgf, lgf_diff
mp.mp.dps = 32
pi = mp.pi
Wsc = mp.sqrt(6)/(32*pi**3)*mp.gamma(mp.mpf(1)/24)*mp.gamma(mp.mpf(5)/24)*mp.gamma(mp.mpf(7)/24)*mp.gamma(mp.mpf(11)/24)
Wf = 3*mp.gamma(mp.mpf(1)/3)**6/(2**(mp.mpf(14)/3)*pi**4)
types = [(1,1,0),(2,0,0),(2,1,1),(2,2,0)]
D  = {x: lgf_diff(x) for x in types}
Gp = {x: lgf(x, -1) for x in [(0,0,0)] + types}
print("Delta(delta) should be 1/6:", D[(1,1,0)], "  G'(0)+G'(delta) should be 1/6:", Gp[(0,0,0)] + Gp[(1,1,0)])
for x in types:
    r1 = mp.pslq([D[x], 1, Wf, 1/(pi**2*Wf)], maxcoeff=10**5, maxsteps=10**6)
    r2 = mp.pslq([Gp[x], 1, Wsc, 1/(pi**2*Wsc)], maxcoeff=10**5, maxsteps=10**6)
    print(x, "Delta:", D[x], r1, "\n        G':", Gp[x], r2)
pickle.dump((D, Gp), open("fcc_vals.pkl", "wb"))
