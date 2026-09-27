import sympy as sp, mpmath as mp, sys, pickle
from gen3d import *
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 15
fcc = [(0, 2, 2), (2, 0, 2), (2, 2, 0)]
nv, E = geometric3(fcc, [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)], 2)
deg = [0] * nv
for (u, v, o) in E: deg[u] += 1; deg[v] += 1
print("pyrochlore nv", nv, "edges", len(E), "deg", deg, flush=True)
Y, arcs, _ = transfer_row(nv, E, 0, 0, False)
print("Y row sum", mp.fsum(Y), flush=True)
YQ, _, _ = transfer_row(nv, E, 0, 0, True)
pickle.dump((arcs, Y, YQ), open("pyro_rows.pkl", "wb"))
pp = mp.fsum(((1 if j == 0 else 0) - Y[j])**2 for j in range(6)) / 4
pm = mp.fsum(((1 if j == 0 else 0) - YQ[j])**2 for j in range(6)) / 4
print("+1:", pp, " -1:", pm, flush=True)
