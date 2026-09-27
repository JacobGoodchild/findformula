import mpmath as mp, sys, pickle, time
from gen import trapping
from builders import LATTICES
mp.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 25
name = sys.argv[1]
bip = {"lieb": True, "dice": True, "truncated_square(4.8.8)": True, "truncated_trihexagonal(4.6.12)": True}.get(name, False)
nv, E = LATTICES[name]()
deg = [0] * nv
for (u, v, o) in E: deg[u] += 1; deg[v] += 1
res = {}
seen = set()
for origin in range(nv):
    if deg[origin] in seen and name in ("lieb", "dice"): continue
    seen.add(deg[origin])
    t = time.time()
    print(f"== {name} origin {origin} (deg {deg[origin]})", flush=True)
    r = trapping(nv, E, origin=origin, start=0, bipartite=bip)
    res[origin] = r[:4]
    print("   time", round(time.time() - t), "s", flush=True)
    if name not in ("lieb", "dice"): break
pickle.dump(res, open(f"res_{name.split('(')[0]}.pkl", "wb"))
