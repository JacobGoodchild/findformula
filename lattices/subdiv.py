"""Test: does subdividing every edge halve flip-flop Grover trapping at original vertices,
and give 1/8 at the midpoints?"""
import numpy as np
from builders import LATTICES
from gridcheck import gridcheck
def subdivide(nv, edges):
    newE = []; k = nv
    for (u, v, o) in edges:
        m = k; k += 1
        newE.append((u, m, (0, 0)))        # midpoint in same cell as u
        newE.append((m, v, o))
    return k, newE
tri = (1, [(0, 0, (1, 0)), (0, 0, (0, 1)), (0, 0, (1, -1))])
honey = (2, [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))])
kag = (3, [(0, 1, (0, 0)), (0, 2, (0, 0)), (1, 2, (0, 0)), (1, 0, (1, 0)), (2, 0, (0, 1)), (1, 2, (1, -1))])
for name, (nv, E) in [("triangular", tri), ("honeycomb", honey), ("kagome", kag)]:
    _, base = gridcheck(nv, E, 0, 60)
    nv2, E2 = subdivide(nv, E)
    _, sub = gridcheck(nv2, E2, 0, 60)
    _, mid = gridcheck(nv2, E2, nv, 60)
    print(name, "base:", {k: round(v, 7) for k, v in list(base.items())[:2]}, "\n  subdivided vertex:", {k: round(v, 7) for k, v in list(sub.items())[:2]}, "\n  midpoint:", {k: round(v, 7) for k, v in mid.items()})
