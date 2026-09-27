import mpmath as mp, pickle
from engine import trap
mp.mp.dps = 35
honey = [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))]
r = trap(2, honey, 0, True)
print("honeycomb pbar:", r[0], " kappa/alpha:", r[5])
print("honeycomb Y row:", [mp.nstr(v, 20) for v in r[3][0]])
kag = [(0, 1, (0, 0)), (0, 2, (0, 0)), (1, 2, (0, 0)), (1, 0, (1, 0)), (2, 0, (0, 1)), (1, 2, (1, -1))]
r2 = trap(3, kag, 0, False)
print("kagome pbar:", r2[0], " +1:", r2[1], " -1:", r2[2], " kappa/alpha:", r2[5])
print("kagome Y row:", [mp.nstr(v, 20) for v in r2[3][0]])
print("kagome YQ row:", [mp.nstr(v, 20) for v in r2[4][0]])
pickle.dump((r[:3], [list(x) for x in r[3]], r2[:3], [list(x) for x in r2[3]], [list(x) for x in r2[4]]), open("hk.pkl", "wb"))
