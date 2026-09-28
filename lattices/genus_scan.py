"""Scan weighted decorations of the square lattice (2x2 supercell) for Laplacian spectral curves that
become genus 0 after removing the node at k = 0: those have elementary resistances and Clausen-type
spanning-tree entropies. Prints the odd-multiplicity factors of the z2-discriminant."""
import sympy as sp
from gen import z1, z2
w = sp.symbols('w', positive=True)
def supercell(bondrule):
    """bondrule(x, y) -> list of (dx, dy, weight) bonds from site (x, y) (each bond listed once)."""
    E = []
    for x in range(2):
        for y in range(2):
            for dx, dy, c in bondrule(x, y):
                X, Y = x + dx, y + dy
                E.append((x + 2 * y, (X % 2) + 2 * (Y % 2), (X // 2, Y // 2), c))
    return 4, E
def lap_det(nv, E):
    A = sp.zeros(nv, nv); deg = [0] * nv
    for (u, v, o, c) in E:
        ph = z1**o[0] * z2**o[1]
        A[u, v] += c * ph; A[v, u] += c / ph; deg[u] += c; deg[v] += c
    return sp.expand((sp.diag(*deg) - A).det())
def odd_factors(d):
    pw = [t.as_powers_dict().get(z2, 0) for t in sp.Add.make_args(d)]
    P = sp.expand(d * z2**(-min(pw)))
    pw1 = [t.as_powers_dict().get(z1, 0) for t in sp.Add.make_args(P)]
    P = sp.expand(P * z1**(-min(pw1)))
    disc = sp.factor(sp.discriminant(P, z2))
    out = []
    for f, m in sp.factor_list(disc)[1]:
        if m % 2 == 1 and f.has(z1) and sp.simplify(f - z1) != 0:
            out.append(f)
    return sp.Poly(P, z2).degree(), out
RULES = {
 "SS (orthogonal dimers)": lambda x, y: [(1, 0, 1), (0, 1, 1)] + ([(1, 1, w)] if (x, y) == (0, 0) else [(-1, 1, w)] if (x, y) == (1, 0) and False else []) ,
 "checkerboard (both diagonals, alternate plaquettes)": lambda x, y: [(1, 0, 1), (0, 1, 1)] + ([(1, 1, w)] if (x + y) % 2 == 0 and x == y == 0 else []) + ([(1, 1, w)] if (x, y) == (1, 1) else []) + ([(-1, 1, w)] if (x, y) == (1, 0) else []) + ([(-1, 1, w)] if (x, y) == (0, 1) else []),
 "king (both diagonals everywhere)": lambda x, y: [(1, 0, 1), (0, 1, 1), (1, 1, w), (-1, 1, w)],
 "square + axis next-nearest": lambda x, y: [(1, 0, 1), (0, 1, 1), (2, 0, w), (0, 2, w)],
 "square + one diagonal on alternate plaquettes": lambda x, y: [(1, 0, 1), (0, 1, 1)] + ([(1, 1, w)] if (x + y) % 2 == 0 else []),
 "square + parallel diagonals on alternate columns": lambda x, y: [(1, 0, 1), (0, 1, 1)] + ([(1, 1, w)] if x == 0 else []),
 "square + diagonals on 1/4 plaquettes": lambda x, y: [(1, 0, 1), (0, 1, 1)] + ([(1, 1, w)] if (x, y) == (0, 0) else []),
}
if __name__ == "__main__":
    for name, rule in RULES.items():
        nv, E = supercell(rule)
        d = lap_det(nv, E)
        deg, fac = odd_factors(d)
        print(f"{name}: z2-degree {deg}, odd factors: {fac}", flush=True)
