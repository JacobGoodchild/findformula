import mpmath as mp, pickle
mp.mp.dps = 40
t = mp.mpf(15)
def gd(n):
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2)
        A = t - c1 * c2; B = c1 + c2; D = A * A - B * B
        return [1 / mp.sqrt(D), -A / D**1.5, (2 * A * A + B * B) / D**2.5][n]
    return mp.quad(f, [0, mp.pi], [0, mp.pi]) / mp.pi**2
gs = [gd(0), gd(1), gd(2)]
print(gs)
pickle.dump(gs, open("fcc_t15.pkl", "wb"))
arcs, YQ = pickle.load(open("pyroQ24.pkl", "rb"))
mp.mp.dps = 24
for a, y in zip(arcs, YQ):
    r = mp.pslq([y, 1, gs[0], gs[1], gs[2]], maxcoeff=10**5, maxsteps=10**6)
    print(a, mp.nstr(y, 20), r)
