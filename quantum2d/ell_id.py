import mpmath as mp
from lgf2d import b
mp.mp.dps = 40; pi = mp.pi
g = b((0,0), lambda k: 11 - mp.cos(k), lambda k: 1 + mp.exp(-1j*k))
print("g =", g)
cands = {}
for name, e in [("5", 5), ("1/5", mp.mpf(1)/5), ("25", 25), ("sqrt5", mp.sqrt(5)), ("3", 3), ("11", 11), ("11/3", mp.mpf(11)/3), ("5/3", mp.mpf(5)/3), ("3/5", mp.mpf(3)/5), ("9/5", mp.mpf(9)/5)]:
    for form, m in [("16e/((1+e)^3(3-e))", 16*e/((1+e)**3*(3-e))), ("16e^3/... inv", ((1+e)**3*(3-e))/(16*e) if e != 3 else None),
                    ("(e-1)^3(e+3)/(16 e)", (e-1)**3*(e+3)/(16*e)), ("16e/((e-1)^3(e+3))", 16*e/((e-1)**3*(e+3)) if e != 1 else None),
                    ("16 e/((e+1)^3 (e-3))", 16*e/((e+1)**3*(e-3)) if e != 3 else None)]:
        if m is None or m <= 0 or m >= 1: continue
        K = mp.ellipk(m)
        ratio = g * pi / K
        p = mp.findpoly(ratio, 4, maxcoeff=2000)
        if p: print("HIT eps=", name, form, "m=", m, "ratio poly", p)
        p = mp.findpoly(ratio**2, 4, maxcoeff=2000)
        if p: print("HIT(sq) eps=", name, form, "m=", m, "ratio^2 poly", p)
