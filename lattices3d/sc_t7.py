import mpmath as mp, pickle
mp.mp.dps = 45
t = 7
def g(n):   # n-th derivative of g(t) = E[1/(t + sum cos)] = int e^{-ts} I0(s)^3 ds
    return mp.quad(lambda s: (-s)**n * mp.exp(-t * s) * mp.besseli(0, s)**3, [0, 1, 5, 20, 80, mp.inf])
gs = [g(0), g(1), g(2)]
def G3(x):  # E[cos(k.x)/(2(7 + sum cos))]
    return mp.quad(lambda s: mp.exp(-t * s) * mp.fprod((-1)**abs(xi) * mp.besseli(abs(xi), s) for xi in x), [0, 1, 5, 20, 80, mp.inf]) / 2
xs = [(0,0,0), (1,0,0), (2,0,0), (1,1,0), (2,1,0), (1,1,1), (3,0,0), (2,1,1), (2,2,0)]
res = {}
for x in xs:
    v = G3(x)
    r = mp.pslq([v, 1, gs[0], gs[1], gs[2]], maxcoeff=10**7, maxsteps=10**7)
    res[x] = (v, r); print(x, mp.nstr(v, 20), r, flush=True)
pickle.dump((gs, res), open("sc_t7.pkl", "wb"))
