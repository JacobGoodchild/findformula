import mpmath as mp
pi = mp.pi
def gamma_products():
    g = lambda a, b: mp.gamma(mp.mpf(a) / b)
    return {"G(1/3)^3/pi^2": g(1,3)**3/pi**2, "G(1/4)^2/pi^1.5": g(1,4)**2/pi**1.5,
            "G(1/8)G(3/8)/pi^1.5": g(1,8)*g(3,8)/pi**1.5, "G(1/24)G(5/24)G(7/24)G(11/24)/pi^3": g(1,24)*g(5,24)*g(7,24)*g(11,24)/pi**3,
            "G(1/3)^6/pi^4": g(1,3)**6/pi**4, "G(1/4)^4/pi^3": g(1,4)**4/pi**3, "(G(1/8)G(3/8))^2/pi^3": (g(1,8)*g(3,8))**2/pi**3,
            "G(1/5)G(2/5)/pi^1.5":g(1,5)*g(2,5)/pi**1.5, "G(1/12)G(5/12)/pi^1.5": g(1,12)*g(5,12)/pi**1.5,
            "G(1/20)G(9/20)/pi^1.5": g(1,20)*g(9,20)/pi**1.5, "G(1/7)G(2/7)G(4/7)/pi^2": g(1,7)*g(2,7)*g(4,7)/pi**2,
            "G(1/15)G(2/15)G(4/15)G(8/15)/pi^3": g(1,15)*g(2,15)*g(4,15)*g(8,15)/pi**3}
def algebraics():
    s = mp.sqrt
    base = {"1":1, "s2":s(2), "s3":s(3), "s5":s(5), "s6":s(6), "s7":s(7), "s10":s(10), "s15":s(15), "2^(1/3)":mp.cbrt(2), "4^(1/3)":mp.cbrt(4),
            "2^(1/4)":mp.root(2,4), "8^(1/4)":mp.root(8,4), "3^(1/4)":mp.root(3,4), "27^(1/4)":mp.root(27,4), "s(1+s2)":s(1+s(2)), "s(2+s3)":s(2+s(3)), "s(s2-1)": s(s(2)-1),
            "s3*2^(1/3)":s(3)*mp.cbrt(2), "s3*4^(1/3)":s(3)*mp.cbrt(4), "5^(1/4)":mp.root(5,4), "s((5+s5)/2)": s((5+s(5))/2)}
    return base
def find_gamma(v, maxc=5000):
    hits = []
    for cn, c in gamma_products().items():
        for an, al in algebraics().items():
            r = mp.pslq([v, c * al], maxcoeff=maxc, maxsteps=10**5)
            if r: hits.append((cn, an, r))
    return hits
def find_lin(v, names, vals, maxc=10**4):
    r = mp.pslq([v] + vals, maxcoeff=maxc, maxsteps=10**6)
    return [(n, c) for n, c in zip(["v"] + names, r) if c] if r else None
