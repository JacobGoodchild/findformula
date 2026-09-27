import mpmath as mp, pickle
from gen import trapping
from builders import LATTICES
from ident import lin, B_sq, B_tri, B_mix
mp.mp.dps = 30
for name, org, B in [("lieb", 0, B_sq), ("lieb", 1, B_sq), ("dice", 0, B_tri), ("dice", 1, B_tri), ("truncated_square(4.8.8)", 0, B_mix)]:
    nv, E = LATTICES[name]()
    tot, pp, pm, extra, Y, YQ = trapping(nv, E, origin=org, start=0, bipartite=True, verbose=False)
    print(name, org, "total", mp.nstr(tot, 25))
    print("   +1 channel:", lin(pp, B()))
    for k, v in extra.items(): print("   flat", k, ":", lin(v, B()))
    print("   total:", lin(tot, B()))
    print("   Y row:", [mp.nstr(y, 15) for y in Y[0]], [lin(y, B()) for y in Y[0]])
