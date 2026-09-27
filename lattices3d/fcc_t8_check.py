import mpmath as mp, sys, pickle
sys.argv = ['x', '60']
exec(open('lg_fcc.py').read().split('cache = {}')[0])
mp.mp.dps = 60
def hd(n, t=mp.mpf(8)):
    def f(k1, k2):
        c1, c2 = mp.cos(k1), mp.cos(k2)
        A = t + c1 * c2; B = c1 + c2; D = A * A - B * B
        return [1 / mp.sqrt(D), -A / D**1.5, (2 * A * A + B * B) / D**2.5][n]
    return mp.quad(f, [0, mp.pi / 2, mp.pi], [0, mp.pi / 2, mp.pi]) / mp.pi**2
hs = [hd(0), hd(1), hd(2)]
rels = {}
for x in [(0,0,0), (0,1,1), (0,2,2), (0,3,3), (0,0,2), (0,1,3), (1,1,2), (1,2,3), (2,2,2)]:
    v = base(x, 8, -1)
    r = mp.pslq([v, 1] + hs, maxcoeff=10**12, maxsteps=10**8)
    rels[x] = r; print(x, r, flush=True)
pickle.dump(rels, open('fcc_t8_rels.pkl', 'wb'))
