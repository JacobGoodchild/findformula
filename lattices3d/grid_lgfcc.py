import itertools, numpy as np
from grid3d import grid3
dirs = [(1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)]
def ends(t, x): return {x, tuple(a + b for a, b in zip(x, dirs[t]))}
edges = set()
for t0 in range(6):
    for t in range(6):
        for x in itertools.product(range(-2, 3), repeat=3):
            if (t, x) == (t0, (0, 0, 0)): continue
            if ends(t, x) & ends(t0, (0, 0, 0)):
                key = (t0, t, x); rk = (t, t0, tuple(-a for a in x))
                if rk not in edges: edges.add(key)
E = sorted(edges)
arcs = []
for (u, v, o) in E:
    if u == 0: arcs.append((v, o))
    if v == 0: arcs.append((u, tuple(-x for x in o)))
print("first arc", arcs[0], len(arcs))
for N in (10, 14):
    u, p = grid3(6, E, 0, N); print(N, sum(p))
