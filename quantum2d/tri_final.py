import mpmath as mp, numpy as np
mp.mp.dps = 40; pi = mp.pi; s3 = mp.sqrt(3)
X = mp.cbrt(4) * mp.gamma(mp.mpf(1)/3)**3 / (8 * pi**2)
Gp = {0: X, 1: (1 - 3*X)/3, "r3": (9*X - 2 - s3/(pi*X))/3, 2: (3*X - 4 + 2*s3/(pi*X))/3}
aa = {0: 0, 1: mp.mpf(1)/3, "r3": 2*s3/pi - mp.mpf(2)/3, 2: mp.mpf(8)/3 - 4*s3/pi}
print("lattice eq at delta:", 3*Gp[1] + (Gp[2] + Gp[0] + 2*Gp["r3"] + 2*Gp[1])/2)
print("lattice eq at 0:", 3*Gp[0] + 6*Gp[1]/2 - 1)
def cls(y):   # classify axial vector by distance class: 0, 1, sqrt3, 2
    x1, x2 = y; n2 = x1*x1 + x2*x2 + x1*x2*(-1)   # |x1 e1 + x2 e2|^2 with e1.e2 = 1/2 -> x1^2+x2^2+x1x2 ; use e1.e2=+1/2
    n2 = x1*x1 + x2*x2 + x1*x2
    return {0: 0, 1: 1, 3: "r3", 4: 2}[n2]
deltas = [(1, 0), (0, 1), (1, -1)]
# NOTE: with neighbours e1, e2, e1-e2, the angle e1.e2 = 1/2
def n2(y): return y[0]**2 + y[1]**2 + y[0]*y[1]
states = [(s, p) for p in range(3) for s in (1, -1)]
Mp = mp.matrix(6, 6); Mm = mp.matrix(6, 6)
for i, (s, p) in enumerate(states):
    for j, (t, q) in enumerate(states):
        y = (t*deltas[q][0] - s*deltas[p][0], t*deltas[q][1] - s*deltas[p][1])
        c = cls(y)
        Mp[i, j] = (1 if i == j else 0)/mp.mpf(2) - (2*aa[1] - aa[c]) / 4
        Mm[i, j] = (1 if i == j else 0)/mp.mpf(2) - (Gp[0] + 2*Gp[1] + Gp[c]) / 4
pp = sum(Mp[k, 0]**2 for k in range(6)); pm = sum(Mm[k, 0]**2 for k in range(6))
print("triangular flip-flop pbar =", pp + pm, "(+1:", pp, " -1:", pm, ")")
