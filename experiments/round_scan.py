"""Scan families of lattice-region graphs for 'round' spanning-tree counts.

For each (lattice, shape, boundary condition) family and size n, build the graph,
count spanning trees exactly (matrix-tree theorem, integer determinant via FLINT),
factor the count, and measure how smooth it is. Families whose counts are always
products of small primes are strong candidates for hidden product formulas.
"""
import sys, json, math
import flint
from sympy import factorint

# ---------- lattices: neighbour offsets on Z^2 coordinates ----------
LATT = {
    "square": [(1, 0), (0, 1)],
    "tri": [(1, 0), (0, 1), (1, -1)],          # triangular lattice in axial coords
    "king": [(1, 0), (0, 1), (1, 1), (1, -1)],
}


def hex_cells(a, b, c):
    """axial coordinates of a hexagon with side lengths a,b,c (triangular lattice points)."""
    pts = set()
    # hexagon: 0<=x<=a+b? use cube coords: x+y+z=0, -a<=... generic semi-regular
    for x in range(-60, 61):
        for y in range(-60, 61):
            z = -x - y
            if 0 <= x <= a + c and 0 <= y + c <= b + c and -a - b <= z <= 0 and True:
                pts.add((x, y))
    return pts


def shape(name, n):
    P = set()
    R = range(-2 * n - 3, 2 * n + 4)
    for x in R:
        for y in R:
            if name == "diamond" and abs(x) + abs(y) <= n: P.add((x, y))
            if name == "triangle" and x >= 0 and y >= 0 and x + y <= n: P.add((x, y))
            if name == "square" and 0 <= x <= n and 0 <= y <= n: P.add((x, y))
            if name == "rect2" and 0 <= x <= 2 * n and 0 <= y <= n: P.add((x, y))
            if name == "trihex" and x >= 0 and y >= 0 and x + y <= 2 * n and x <= 3 * n // 2 + 0 * y and y <= 3 * n // 2: P.add((x, y))
            if name == "hexreg" and abs(x) <= n and abs(y) <= n and abs(x + y) <= n: P.add((x, y))
            if name == "cross" and ((abs(x) <= n and abs(y) <= n // 2 + 0) or (abs(y) <= n and abs(x) <= n // 2)): P.add((x, y))
            if name == "octagon" and abs(x) <= n and abs(y) <= n and abs(x) + abs(y) <= 3 * n // 2: P.add((x, y))
            if name == "stair" and x >= 0 and y >= 0 and x + 2 * y <= 2 * n: P.add((x, y))
            if name == "halfdiamond" and y >= 0 and abs(x) + y <= n: P.add((x, y))
            if name == "annulus" and max(abs(x), abs(y)) <= n and max(abs(x), abs(y)) >= 1: P.add((x, y))
    return P


def graph(latt, pts, wired):
    idx = {p: i for i, p in enumerate(sorted(pts))}
    edges = []
    boundary = set()
    for (x, y) in pts:
        for dx, dy in LATT[latt]:
            for s in (1, -1):
                q = (x + s * dx, y + s * dy)
                if q in idx:
                    if s == 1:
                        edges.append((idx[(x, y)], idx[q]))
                else:
                    boundary.add(idx[(x, y)])
    n = len(idx)
    if wired:  # add one extra vertex joined to each missing neighbour (one edge per missing lattice edge)
        W = n
        for (x, y) in pts:
            for dx, dy in LATT[latt]:
                for s in (1, -1):
                    q = (x + s * dx, y + s * dy)
                    if q not in idx:
                        edges.append((idx[(x, y)], W))
        n += 1
    return n, edges


def spanning_trees(n, edges):
    if n == 1:
        return 1
    L = [[0] * n for _ in range(n)]
    for u, v in edges:
        L[u][u] += 1; L[v][v] += 1
        L[u][v] -= 1; L[v][u] -= 1
    M = flint.fmpz_mat([row[1:] for row in L[1:]])
    return int(M.det())


from sympy import primerange
SMALLP = list(primerange(2, 20000))
def smooth_info(N):
    """largest prime factor below 20000, or the rough cofactor if one remains"""
    N = int(N); pm = 1
    for p in SMALLP:
        if N % p == 0:
            pm = p
            while N % p == 0: N //= p
    if N > 1:
        return N, None
    return pm, None


if __name__ == "__main__":
    out = open("round_scan.jsonl", "w")
    for latt in LATT:
        for sh in ["diamond", "triangle", "square", "rect2", "hexreg", "cross", "octagon", "stair", "halfdiamond", "annulus", "trihex"]:
            for wired in (False, True):
                rows = []
                for n in range(1, 13):
                    pts = shape(sh, n)
                    if len(pts) > 700:
                        break
                    nv, E = graph(latt, pts, wired)
                    T = spanning_trees(nv, E)
                    pm, f = smooth_info(T)
                    rows.append((n, len(pts), pm, T.bit_length()))
                score = max((math.log(r[2]) / (r[3] * math.log(2)) for r in rows if r[3] > 20), default=1)
                rec = {"latt": latt, "shape": sh, "wired": wired, "rows": rows, "score": score}
                out.write(json.dumps(rec) + "\n"); out.flush()
                print(f"{latt:7s} {sh:12s} wired={wired!s:5s} score={score:.3f} maxp={[r[2] for r in rows]}", flush=True)
