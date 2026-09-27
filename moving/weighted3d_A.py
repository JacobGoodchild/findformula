import mpmath as mp, itertools
mp.mp.dps = 60
def A_of(p): return mp.quad(lambda u: mp.fprod(mp.erfc(mp.sqrt(pk * u)) for pk in p), [0, 1, 5, 20, 80, mp.inf])
def pi1(p):
    p1, p2, p3 = p
    return 2 / mp.pi * (mp.atan(mp.sqrt(p1 / p2)) + mp.atan(mp.sqrt(p1 / p3)) - mp.atan(mp.sqrt(p1 / (p2 * p3))))
for p in [(mp.sqrt(2)/10, mp.sqrt(7)/10, 1 - mp.sqrt(2)/10 - mp.sqrt(7)/10), (mp.sqrt(3)/9, mp.sqrt(5)/8, 1 - mp.sqrt(3)/9 - mp.sqrt(5)/8)]:
    A = A_of(p)
    # check pi1 formula
    num = mp.quad(lambda s: mp.exp(-s) / mp.sqrt(mp.pi * s) * mp.erfc(mp.sqrt(p[1] * s / p[0])) * mp.erfc(mp.sqrt(p[2] * s / p[0])), [0, 1, 5, 20, 80, mp.inf])
    print("pi1 formula check:", mp.nstr(pi1(p) - num, 5))
    names, vals = [], []
    for i in range(3):
        j, k = [x for x in range(3) if x != i]
        for nm, v in [("1", 1), ("at(pi/pj)", mp.atan(mp.sqrt(p[i] / p[j])) / mp.pi), ("at(pi/pk)", mp.atan(mp.sqrt(p[i] / p[k])) / mp.pi),
                      ("at(pi/(pjpk))", mp.atan(mp.sqrt(p[i] / (p[j] * p[k]))) / mp.pi)]:
            names.append(nm + "/p%d" % (i + 1)); vals.append(v / p[i])
    names.append("sqrt(p1p2p3)/pi"); vals.append(mp.sqrt(p[0] * p[1] * p[2]) / mp.pi)
    names.append("1/(pi sqrt(p1p2p3))"); vals.append(1 / (mp.pi * mp.sqrt(p[0] * p[1] * p[2])))
    r = mp.pslq([A] + vals, maxcoeff=10**4, maxsteps=10**7)
    print("A =", mp.nstr(A, 25), [(n, c) for n, c in zip(["A"] + names, r) if c] if r else None)
