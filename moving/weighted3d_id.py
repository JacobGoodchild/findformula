import mpmath as mp
mp.mp.dps = 35
def piA(p):
    p = [mp.mpf(x) for x in p]
    A = mp.quad(lambda u: mp.fprod(mp.erfc(mp.sqrt(pk * u)) for pk in p), [0, 1, 5, 20, 80, mp.inf])
    pis = []
    for i in range(3):
        o = [j for j in range(3) if j != i]
        pis.append(mp.quad(lambda s: mp.exp(-s) / mp.sqrt(mp.pi * s) * mp.fprod(mp.erfc(mp.sqrt(p[j] * s / p[i])) for j in o), [0, 1, 5, 20, 80, mp.inf]))
    return A, pis
for p in [(mp.mpf(1)/6, mp.mpf(1)/3, mp.mpf(1)/2), (mp.mpf(1)/5, mp.mpf(3)/10, mp.mpf(1)/2), (mp.mpf(3)/20, mp.mpf(7)/20, mp.mpf(1)/2), (mp.mpf(1)/3,)*3]:
    A, pis = piA(p)
    p1, p2, p3 = p
    at = lambda z: mp.atan(mp.sqrt(z)) / mp.pi
    B = {"1": 1, "at(p2/p1)": at(p2/p1), "at(p3/p1)": at(p3/p1), "at(p2p3/p1)": at(p2*p3/p1), "at(p1p3/p2)": at(p1*p3/p2)}
    names = list(B)
    r = mp.pslq([pis[0]] + [B[n] for n in names], maxcoeff=1000, maxsteps=10**6)
    print("p", [mp.nstr(x, 4) for x in p], "pi1", mp.nstr(pis[0], 20), [(n, c) for n, c in zip(["v"] + names, r) if c] if r else None)
    # A: allow algebraic/pi terms
    s = mp.sqrt
    BA = {"1": 1, "1/p1": 1/p1, "1/p2": 1/p2, "1/p3": 1/p3}
    for nm in ["at(p2p3/p1)", "at(p1p3/p2)"]:
        for q, qn in [(1/p1, "/p1"), (1/p2, "/p2"), (1/p3, "/p3")]:
            BA[nm + qn] = B[nm] * q
    BA["s(p1p2p3)/pi/(p1p2p3)"] = s(p1*p2*p3) / mp.pi / (p1*p2*p3)
    names = list(BA)
    r = mp.pslq([A] + [BA[n] for n in names], maxcoeff=10**4, maxsteps=10**7)
    print("    A", mp.nstr(A, 20), [(n, c) for n, c in zip(["v"] + names, r) if c] if r else None)
