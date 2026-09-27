import itertools, numpy as np
from grid3d import grid3
dirs = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
def ends(t, x): return {x, tuple(a + b for a, b in zip(x, dirs[t]))}
edges = set()
for t0 in range(3):
    for t in range(3):
        for x in itertools.product(range(-2, 3), repeat=3):
            if (t, x) == (t0, (0, 0, 0)): continue
            if ends(t, x) & ends(t0, (0, 0, 0)):
                key = (t0, t, x); rk = (t, t0, tuple(-a for a in x))
                if rk not in edges: edges.add(key)
E = sorted(edges)
import grid3d
# grid3 returns per-flat-eigenvalue contributions for the FIRST out-arc; loop over arcs by reordering
arcs0 = [e for e in E if e[0] == 0] + [(e[1], e[0], tuple(-a for a in e[2])) for e in E if e[1] == 0 and e[0] != 0]
for N in (24, 36):
    u, p = grid3(3, E, 0, N); print(N, "flat eigs", np.round(u, 4), "first-arc total", sum(p))
