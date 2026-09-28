"""Abelian sandpile on an N x N king lattice (degree 8, open boundary): heights 1..8, topple when h > 8
(equivalently h >= 9 -> give 1 grain to each of 8 neighbours). Measures density of height 1 in the bulk."""
import numpy as np, sys
rng = np.random.default_rng(1)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 128
h = rng.integers(1, 9, size=(N, N))
SH = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)]
def relax(h):
    while True:
        t = h > 8
        if not t.any(): return h
        tt = t.astype(np.int64)
        h = h - 8 * tt
        for dx, dy in SH:
            s = np.zeros_like(tt)
            xs = slice(max(dx, 0), N + min(dx, 0)); xd = slice(max(-dx, 0), N + min(-dx, 0))
            ys = slice(max(dy, 0), N + min(dy, 0)); yd = slice(max(-dy, 0), N + min(-dy, 0))
            s[xs, ys] = tt[xd, yd]
            h = h + s
h = relax(h + 8)                      # drive into the recurrent class
for _ in range(20000):                # burn-in
    i, j = rng.integers(0, N, 2); h[i, j] += 1; h = relax(h) if h[i, j] > 8 else h
acc = []; c = N // 4
for it in range(200000):
    i, j = rng.integers(0, N, 2); h[i, j] += 1
    if h[i, j] > 8: h = relax(h)
    if it % 200 == 0: acc.append((h[c:-c, c:-c] == 1).mean())
a = np.array(acc); print(N, "P(h=1) bulk =", a.mean(), "+-", a.std() / np.sqrt(len(a) / 10))
