"""Count perfect matchings of lattice regions exactly (frontier dynamic programming)
and look for families with 'round' (smooth) counts -> candidate product formulas."""
import json, math, sys
from sympy import primerange

LATT = {
    "square": [(1, 0), (0, 1)],
    "tri": [(1, 0), (0, 1), (1, -1)],
    "king": [(1, 0), (0, 1), (1, 1), (1, -1)],
    "honey": None,  # honeycomb handled as brick-wall: square lattice minus alternate vertical edges
}
SMALLP = list(primerange(2, 5000))


def nbrs(latt, p, pts):
    x, y = p
    out = []
    if latt == "honey":
        cand = [(x + 1, y), (x - 1, y)] + ([(x, y + 1)] if (x + y) % 2 == 0 else [(x, y - 1)])
    else:
        cand = []
        for dx, dy in LATT[latt]:
            cand += [(x + dx, y + dy), (x - dx, y - dy)]
    return [q for q in cand if q in pts]


def count_pm(latt, pts):
    order = sorted(pts)
    pos = {p: i for i, p in enumerate(order)}
    later = [[pos[q] for q in nbrs(latt, p, pts) if pos[q] > i] for i, p in enumerate(order)]
    states = {frozenset(): 1}
    for i in range(len(order)):
        new = {}
        for st, c in states.items():
            if i in st:
                ns = st - {i}
                new[ns] = new.get(ns, 0) + c
            else:
                for j in later[i]:
                    if j not in st:
                        ns = st | {j}
                        new[ns] = new.get(ns, 0) + c
        states = new
        if not states:
            return 0
    return states.get(frozenset(), 0)


def rough(N):
    if N == 0:
        return 0
    pm = 1
    for p in SMALLP:
        if N % p == 0:
            pm = p
            while N % p == 0:
                N //= p
    return N if N > 1 else pm


def region(name, n, k):
    P = set()
    R = range(-3 * n - 6, 3 * n + 7)
    for x in R:
        for y in R:
            ok = {
                "diamond": abs(x) + abs(y) <= n,
                "tri": x >= 0 and y >= 0 and x + y <= n,
                "rect": 0 <= x < n and 0 <= y < k,
                "hex": abs(x) <= n and abs(y) <= n and abs(x + y) <= n,
                "hexk": 0 <= x < n + k and 0 <= y < n + k and n - 1 <= x + y <= 2 * n + k - 2 + 0,
                "par": 0 <= x < n and 0 <= y < k,
                "trap": y >= 0 and y < k and x >= 0 and x + y < n,
                "stair": x >= 0 and y >= 0 and x + 2 * y < 2 * n,
                "adiam": abs(x - 0.5) + abs(y - 0.5) <= n,
                "octa": abs(x) <= n and abs(y) <= n and abs(x) + abs(y) <= n + k,
            }[name]
            if ok:
                P.add((x, y))
    return P


if __name__ == "__main__":
    latts = sys.argv[1].split(",")
    out = open(f"pm_scan_{'_'.join(latts)}.jsonl", "w")
    for latt in latts:
        for name in ["diamond", "adiam", "tri", "hex", "stair", "rect", "trap", "hexk", "octa"]:
            two = name in ("rect", "trap", "hexk", "octa")
            for k in (range(1, 7) if two else [0]):
                rows = []
                for n in range(1, 40):
                    pts = region(name, n, k)
                    if len(pts) > 260:
                        break
                    if len(pts) % 2:
                        continue
                    c = count_pm(latt, pts)
                    rows.append((n, len(pts), c, rough(c)))
                good = [r for r in rows if r[2] > 1000]
                smooth = all(r[3] < 5000 for r in good) and len(good) >= 3
                rec = {"latt": latt, "shape": name, "k": k, "rows": [(a, b, str(c), str(d)) for a, b, c, d in rows], "smooth": smooth}
                out.write(json.dumps(rec) + "\n"); out.flush()
                print(latt, name, k, "SMOOTH" if smooth else "", [r[2] for r in rows][:12], flush=True)
