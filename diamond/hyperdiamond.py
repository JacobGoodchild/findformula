import sys; from projgrid3d import pbar
d, N = int(sys.argv[1]), int(sys.argv[2])
E = [(0, 1, tuple([0] * d))] + [(0, 1, tuple(-1 if i == j else 0 for i in range(d))) for j in range(d)]
r = pbar(2, E, 0, 0, N); big = {k: v for k, v in r.items() if v > 1e-4}
D = d + 1
print(d, big, "sum", sum(big.values()), "predicted", (D - 2)**2 / (2 * D * (D - 1)))
