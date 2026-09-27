import sys; from projgrid3d import pbar
# K4 crystal (Laves graph, srs net) = maximal abelian cover of K4: tree edges 0-1,0-2,0-3; others carry e1,e2,e3
srs = [(0, 1, (0, 0, 0)), (0, 2, (0, 0, 0)), (0, 3, (0, 0, 0)), (1, 2, (1, 0, 0)), (2, 3, (0, 1, 0)), (3, 1, (0, 0, 1))]
N = int(sys.argv[1])
for org in (0, 1):
    r = pbar(4, srs, org, 0, N); big = {k: v for k, v in r.items() if v > 1e-4}
    print(org, big, "sum(all)", sum(r.values()))
