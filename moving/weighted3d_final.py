import mpmath as mp
mp.mp.dps = 40
def pis_A(p):
    p = [mp.mpf(x) for x in p]; s = mp.sqrt; at = mp.atan
    pis = []; A = 0
    for i in range(3):
        j, k = [x for x in range(3) if x != i]
        pi_i = 2 / mp.pi * (at(s(p[i] / p[j])) + at(s(p[i] / p[k])) - at(s(p[i] / (p[j] * p[k]))))
        pis.append(pi_i)
        A += pi_i / p[i] - 2 / mp.pi * (s(p[k] / p[i]) / (1 + s(p[j])) + s(p[j] / p[i]) / (1 + s(p[k])))
    return pis, A / 2
def pbar(p):
    pis, A = pis_A(p); p1 = mp.mpf(p[0])
    return 2 * (pis[0]**2 / 4 + (p1 * A - pis[0] / 2)**2 + p1 * (1 - p1) * A**2 / 2)
for p in [(mp.sqrt(2)/10, mp.sqrt(7)/10, 1 - mp.sqrt(2)/10 - mp.sqrt(7)/10), (mp.mpf('0.2'), mp.mpf('0.3'), mp.mpf('0.5')), (mp.mpf(1)/3,)*3]:
    pis, A = pis_A(p)
    num = mp.quad(lambda u: mp.fprod(mp.erfc(mp.sqrt(pk * u)) for pk in p), [0, 0.5, 1, 2, 5, 10, 20, 40, 80, mp.inf])
    print([mp.nstr(x, 6) for x in p], " A formula - numeric:", mp.nstr(A - num, 3), " sum pi_i:", mp.nstr(sum(pis), 20), " pbar:", mp.nstr(pbar(p), 20))
