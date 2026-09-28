"""Systematic search: square lattice, 2x2 supercell, each of the 4 plaquettes gets no diagonal, '/', '\\' or 'X'.
Flag decorations whose Laplacian spectral curve has genus 0 after removing the node (at most one quadratic
odd factor in the z2-discriminant): these have elementary resistances and Clausen-type tree entropies."""
import sympy as sp, itertools, json, sys
from gen import z1, z2
from genus_scan2 import lapdet, odd_deg
def build(conf):
    E = []
    def add(p, q):
        (x, y), (X, Y) = p, q
        E.append(((x % 2) + 2 * (y % 2), (X % 2) + 2 * (Y % 2), (X // 2 - x // 2, Y // 2 - y // 2), 1))
    for x in range(2):
        for y in range(2):
            add((x, y), (x + 1, y)); add((x, y), (x, y + 1))
            c = conf[x + 2 * y]
            if c in "/X": add((x, y), (x + 1, y + 1))
            if c in "\\X": add((x + 1, y), (x, y + 1))
    return 4, E
def canon(conf):
    # symmetry reduction: translations of the 2x2 pattern and the dihedral group (swap / and \ under reflections)
    grid = [[conf[x + 2 * y] for x in range(2)] for y in range(2)]
    best = None
    for rx in range(2):
        for ry in range(2):
            g = [[grid[(y + ry) % 2][(x + rx) % 2] for x in range(2)] for y in range(2)]
            for refl in (False, True):
                h = [[{"/": "\\", "\\": "/"}.get(c, c) if refl else c for c in (row[::-1] if refl else row)] for row in g]
                for tr in (False, True):
                    k = [[h[x][y] for x in range(2)] for y in range(2)] if tr else h
                    if tr: k = [[{"/": "/", "\\": "\\"}.get(c, c) for c in row] for row in k]
                    s = "".join(k[y][x] for y in range(2) for x in range(2))
                    best = s if best is None or s < best else best
    return best
if __name__ == "__main__":
    seen = set(); out = []
    for conf in itertools.product("0/\\X", repeat=4):
        c = canon("".join(conf))
        if c in seen: continue
        seen.add(c)
        nv, E = build(c)
        try:
            deg, odd = odd_deg(lapdet(nv, E))
        except Exception as e:
            deg, odd = None, str(e)
        flag = isinstance(odd, list) and sum(odd) <= 2
        print(c.replace("0", "."), deg, odd, "GENUS0" if flag else "", flush=True)
        out.append((c, deg, odd, flag))
    json.dump(out, open("plaquette_search.json", "w"))
