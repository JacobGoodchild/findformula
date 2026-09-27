import mpmath as mp
from flipflop2d import trap
from identify import find_gamma, find_lin
mp.mp.dps = 40; pi = mp.pi
deltas = [(1,0),(0,1),(1,1),(1,-1)]
A = lambda k: 4 - mp.cos(k); C = lambda k: 1 + 2*mp.cos(k)
Ap = lambda k: 4 + mp.cos(k); Cp = lambda k: -(1 + 2*mp.cos(k))
tot, pp, pm, ca, cg = trap(deltas, A, C, Ap, Cp)
print("king flip-flop pbar:", tot, "+1:", pp, "-1:", pm)
seen = {}
for y, v in sorted(ca.items()):
    k = mp.nstr(v, 25)
    if k in seen: continue
    seen[k] = y
    print("a", y, v, find_lin(v, ["1","1/pi","s2/pi","s3/pi","1/pi^2"], [1, 1/pi, mp.sqrt(2)/pi, mp.sqrt(3)/pi, 1/pi**2]))
print("G'(0) =", cg[(0,0)], find_gamma(cg[(0,0)]))
for y, v in sorted(cg.items()):
    print("G'", y, v)
