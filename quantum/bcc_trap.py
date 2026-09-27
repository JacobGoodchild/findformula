import mpmath as mp
from bcc_lgf import G
mp.mp.dps = 40
Wb = mp.gamma(mp.mpf(1) / 4)**4 / (4 * mp.pi**3)
vals = {x: G(x) for x in [(0,0,0), (1,1,1), (2,0,0), (2,2,0), (2,2,2)]}
g0, g1 = vals[(0,0,0)], vals[(1,1,1)]
alpha = (g0 - 2*g1 + vals[(2,2,2)]) / 4       # partner entry (same pair, opposite direction)
beta_s = (g0 - 2*g1 + vals[(2,0,0)]) / 4      # different pairs, difference vector type (2,0,0)
beta_d = (g0 - 2*g1 + vals[(2,2,0)]) / 4      # different pairs, type (2,2,0)
print("check G(0)-G(d) = 1/4:", g0 - g1)
for name, v in [("alpha", alpha), ("beta_s", beta_s), ("beta_d", beta_d)]:
    print(name, v, mp.pslq([v, 1, Wb, 1 / (mp.pi**2 * Wb)], maxcoeff=10**6, maxsteps=10**6))
# start in one coin state: M e = (1/2 - 1/8) on diagonal, -alpha partner, -beta_s x3, -beta_d x3
Me2 = (mp.mpf(3)/8)**2 + alpha**2 + 3*beta_s**2 + 3*beta_d**2
print("BCC flip-flop, start in one coin state: time-avg return prob =", 2 * Me2)
print("row-sum check (should be 1/2):", mp.mpf(1)/8 + alpha + 3*beta_s + 3*beta_d)
