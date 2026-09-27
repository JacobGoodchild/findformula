import mpmath as mp, pickle, sys
from gen3d import *
mp.mp.dps = 24
fcc = [(0, 2, 2), (2, 0, 2), (2, 2, 0)]
nv, E = geometric3(fcc, [(0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)], 2)
YQ, arcs, _ = transfer_row(nv, E, 0, 0, True)
pickle.dump((arcs, YQ), open("pyroQ24.pkl", "wb"))
print("done")
