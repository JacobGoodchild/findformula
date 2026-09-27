"""Count linear extensions of posets built from lattice regions and flag families with
smooth ('round') counts, which signal hook-length-type product formulas."""
import sys, json
from functools import lru_cache
from sympy import primerange

SMALLP = list(primerange(2, 3000))
DIRS = {
    "sq": [(1, 0), (0, 1)],
    "tri": [(1, 0), (0, 1), (1, -1)],
    "sqd": [(1, 0), (0, 1), (1, 1)],
}


def rough(N):
    pm = 1
    for p in SMALLP:
        if N % p == 0:
            pm = p
            while N % p == 0:
                N //= p
    return N if N > 1 else pm


def count_linext(cells, dirs):
    cells = sorted(cells)
    idx = {c: i for i, c in enumerate(cells)}
    n = len(cells)
    preds = [0] * n
    for (x, y), i in idx.items():
        for dx, dy in dirs:
            q = (x + dx, y + dy)
            if q in idx:
                preds[idx[q]] |= 1 << i
    # forward DP over order ideals, layer by layer
    layer = {0: 1}
    for _ in range(n):
        new = {}
        for ideal, c in layer.items():
            for j in range(n):
                if not (ideal >> j) & 1 and (preds[j] & ~ideal) == 0:
                    k = ideal | (1 << j)
                    new[k] = new.get(k, 0) + c
        layer = new
    return layer[(1 << n) - 1]


def shape(name, n, k):
    S = set()
    R = range(-2 * n - 2 * k - 3, 2 * n + 2 * k + 4)
    for x in R:
        for y in R:
            ok = {
                "diamond": abs(x) + abs(y) <= n,
                "aztec": abs(x - 0.5) + abs(y - 0.5) <= n,
                "hex": abs(x) <= n and abs(y) <= n and abs(x + y) <= n,
                "tri": x >= 0 and y >= 0 and x + y <= n,
                "trap": 0 <= y < k and 0 <= x and x + y < n + k,
                "rect": 0 <= x < n and 0 <= y < k,
                "skewstair": x >= 0 and y >= 0 and n - 1 <= x + y <= n - 1 + k,  # staircase band
                "hexk": 0 <= x < n + k and 0 <= y < n + k and n - 1 <= x + y <= n + 2 * k - 1 + n - 1,
                "octa": abs(x) <= n and abs(y) <= n and abs(x) + abs(y) <= n + k,
                "cross": (abs(x) <= n and abs(y) <= k) or (abs(y) <= n and abs(x) <= k),
                "frame": 0 <= x < n + 2 * k and 0 <= y < n + 2 * k and not (k <= x < n + k and k <= y < n + k),
            }[name]
            if ok:
                S.add((x, y))
    return S


if __name__ == "__main__":
    out = open("linext_scan.jsonl", "w")
    for d in DIRS:
        for name in ["diamond", "aztec", "hex", "tri", "trap", "skewstair", "hexk", "octa", "cross", "frame", "rect"]:
            two = name in ("trap", "skewstair", "hexk", "octa", "cross", "frame", "rect")
            for k in (range(1, 6) if two else [0]):
                rows = []
                for n in range(1, 30):
                    S = shape(name, n, k)
                    if len(S) > 34 or len(S) == 0:
                        break
                    try:
                        c = count_linext(S, DIRS[d])
                    except KeyError:
                        c = 0
                    rows.append((n, len(S), c, rough(c) if c else 0))
                good = [r for r in rows if r[2] > 10 ** 3]
                smooth = len(good) >= 3 and all(r[3] < 3000 for r in good)
                out.write(json.dumps({"dirs": d, "shape": name, "k": k, "rows": [(a, b, str(c), str(e)) for a, b, c, e in rows], "smooth": smooth}) + "\n")
                out.flush()
                print(d, name, k, "SMOOTH" if smooth else "", [(r[1], r[3]) for r in rows], flush=True)
