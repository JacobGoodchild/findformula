import mpmath as mp, numpy as np, itertools
mp.mp.dps = 40; pi = mp.pi; s3 = mp.sqrt(3)
Wsc = mp.sqrt(6)/(32*pi**3)*mp.gamma(mp.mpf(1)/24)*mp.gamma(mp.mpf(5)/24)*mp.gamma(mp.mpf(7)/24)*mp.gamma(mp.mpf(11)/24)
Y = Wsc * mp.sqrt(2)
Wf = 3*mp.gamma(mp.mpf(1)/3)**6/(2**(mp.mpf(14)/3)*pi**4)
Delta = {(0,0,0): 0, (0,1,1): mp.mpf(1)/6, (0,0,2): (4*Wf - 3/(pi**2*Wf))/6,
         (1,1,2): (2*Wf + 3/(pi**2*Wf) - 1)/3, (0,2,2): (8 - 12*Wf - 9/(pi**2*Wf))/3}
Gp = {(0,0,0): Y/12, (0,1,1): (2 - Y)/12, (0,0,2): (7*Y + 108/(pi**2*Y) - 36*s3/pi)/72,
      (1,1,2): (2*Y - 3 - 27/(pi**2*Y))/9, (0,2,2): (Y + 54/(pi**2*Y) + 6*s3/pi - 8)/6}
key = lambda y: tuple(sorted(abs(v) for v in y))
# lattice-equation check for G' at x = delta: 6 G'(d) + 1/2 sum_{12 nbrs} G'(d + e) = 0
nb = [v for v in itertools.product([-1,0,1], repeat=3) if sorted(map(abs, v)) == [0,1,1]]
print("G' lattice equation residual at delta:", 6*Gp[(0,1,1)] + sum(Gp[key(np.add((1,1,0), e))] for e in nb)/2)
print("G' lattice equation residual at 0:", 6*Gp[(0,0,0)] + sum(Gp[key(e)] for e in nb)/2 - 1)
deltas = [(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(0,1,1),(0,1,-1)]
states = [(s, p) for p in range(6) for s in (1, -1)]
Np = mp.matrix(12, 12); Nm = mp.matrix(12, 12)
for a, (s, p) in enumerate(states):
    for b, (t, q) in enumerate(states):
        y = tuple(t*deltas[q][i] - s*deltas[p][i] for i in range(3))
        Np[a, b] = (mp.mpf(1)/3 - Delta[key(y)]) / 4
        Nm[a, b] = (Gp[(0,0,0)] + 2*Gp[(0,1,1)] + Gp[key(y)]) / 4
e = mp.matrix(12, 1); e[0] = 1
Mp = mp.eye(12)/2 - Np; Mm = mp.eye(12)/2 - Nm
v1, v2 = Mp*e, Mm*e
pbar = sum(x**2 for x in v1) + sum(x**2 for x in v2)
print("FCC flip-flop, start in one coin state: pbar =", pbar)
print("  (+1 channel:", sum(x**2 for x in v1), ", -1 channel:", sum(x**2 for x in v2), ")")
