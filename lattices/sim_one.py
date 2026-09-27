import sys, numpy as np
from simgraph import simulate
from builders import LATTICES
name, origin, T = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
nv, E = LATTICES[name]()
r, a = simulate(nv, E, origin, T)
print(f"{name} origin {origin}: avg t={T//2}..{T} = {r[T//2:].mean():.7f}  avg last quarter = {r[3*T//4:].mean():.7f}  arc {a}")
