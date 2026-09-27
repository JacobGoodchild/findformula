import sys
from simgraph import simulate
from builders import LATTICES
from gen import out_arcs
name, origin, T = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
nv, E = LATTICES[name]()
print("engine arc order:", out_arcs(E, origin))
n = len(out_arcs(E, origin))
for s in range(n):
    r, a = simulate(nv, E, origin, T, start_arc=s)
    print(f"start {s}: arc {a}  avg t={T//2}..{T}: {r[T//2:].mean():.7f}", flush=True)
