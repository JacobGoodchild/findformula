import mpmath as mp
from fcc import lgf
mp.mp.dps = 40
g = lgf((0,0,0), -1)
print("G'(0) =", g)
open("fcc_gp0.txt","w").write(str(g))
pi = mp.pi
cands = {
 "G(1/3)^6/pi^4": mp.gamma(mp.mpf(1)/3)**6/pi**4,
 "G(1/4)^4/pi^3": mp.gamma(mp.mpf(1)/4)**4/pi**3,
 "G(1/8)^2G(3/8)^2/pi^3": (mp.gamma(mp.mpf(1)/8)*mp.gamma(mp.mpf(3)/8))**2/pi**3,
 "G(1/24)G(5/24)G(7/24)G(11/24)/pi^3": mp.gamma(mp.mpf(1)/24)*mp.gamma(mp.mpf(5)/24)*mp.gamma(mp.mpf(7)/24)*mp.gamma(mp.mpf(11)/24)/pi**3,
}
algs = {"1":1, "2^(1/3)":mp.cbrt(2), "4^(1/3)":mp.cbrt(4), "sqrt2":mp.sqrt(2), "sqrt3":mp.sqrt(3), "sqrt6":mp.sqrt(6), "3^(1/4)":mp.root(3,4), "sqrt3*2^(1/3)":mp.sqrt(3)*mp.cbrt(2), "sqrt3*4^(1/3)":mp.sqrt(3)*mp.cbrt(4), "2^(1/4)":mp.root(2,4), "sqrt(2+sqrt3)":mp.sqrt(2+mp.sqrt(3)), "sqrt(1+sqrt2)":mp.sqrt(1+mp.sqrt(2))}
for cn, c in cands.items():
    for an, a in algs.items():
        X = c * a
        r = mp.pslq([g, 1, X, 1/(pi**2*X)], maxcoeff=5000, maxsteps=10**5)
        if r: print("HIT", cn, an, r)
        r = mp.pslq([g, X], maxcoeff=10**4, maxsteps=10**5)
        if r: print("HIT ratio", cn, an, r)
print("done")
