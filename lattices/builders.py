"""Unit cells (vertices per cell, edges (u, v, offset)) for 2D lattices."""
import itertools, math
import numpy as np

def geometric(a1, a2, basis, bond=1.0, tol=1e-6):
    a1, a2 = np.array(a1, float), np.array(a2, float)
    B = [np.array(b, float) for b in basis]
    edges = set()
    for u, ru in enumerate(B):
        for v, rv in enumerate(B):
            for o1, o2 in itertools.product(range(-2, 3), repeat=2):
                d = rv + o1 * a1 + o2 * a2 - ru
                if abs(np.linalg.norm(d) - bond) < tol:
                    key = (u, v, (o1, o2)); rk = (v, u, (-o1, -o2))
                    if rk not in edges and not (u == v and (o1, o2) == (0, 0)):
                        edges.add(key)
    return len(B), sorted(edges)

s3 = math.sqrt(3); s2 = math.sqrt(2)

def lieb():   # A (corner, deg 4), B, C (edge centres, deg 2)
    return 3, [(0, 1, (0, 0)), (0, 2, (0, 0)), (1, 0, (1, 0)), (2, 0, (0, 1))]
def dice():   # hub (deg 6) = 0, rims 1, 2 (deg 3): hub connects to all 6 rim sites around it
    return 3, [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1)), (0, 2, (-1, 0)), (0, 2, (0, -1)), (0, 2, (-1, -1))]
def checkerboard():   # line graph of square lattice: 0 = horizontal edge (x,y)-(x+1,y), 1 = vertical edge (x,y)-(x,y+1)
    E = [(0, 0, (0, 1)), (0, 0, (1, 0)),     # (in line graph) h-edges sharing vertices? handled below
         ]
    # adjacency in line graph: edges sharing an endpoint
    h = lambda x, y: ("h", x, y); v = lambda x, y: ("v", x, y)
    def ends(e):
        t, x, y = e
        return {(x, y), (x + 1, y)} if t == "h" else {(x, y), (x, y + 1)}
    cands = [h(x, y) for x in range(-2, 3) for y in range(-2, 3)] + [v(x, y) for x in range(-2, 3) for y in range(-2, 3)]
    edges = set()
    for e0 in (h(0, 0), v(0, 0)):
        for e in cands:
            if e != e0 and ends(e) & ends(e0):
                u = 0 if e0[0] == "h" else 1; w = 0 if e[0] == "h" else 1
                off = (e[1] - e0[1], e[2] - e0[2])
                key = (u, w, off); rk = (w, u, (-off[0], -off[1]))
                if rk not in edges: edges.add(key)
    return 2, sorted(edges)
def shastry_sutherland():   # = snub square 3^2.4.3.4 topologically; 2x2 supercell of square lattice
    sites = [(0, 0), (1, 0), (0, 1), (1, 1)]
    idx = {s: i for i, s in enumerate(sites)}
    def cellsite(x, y):
        return idx[(x % 2, y % 2)], (x // 2, y // 2)
    edges = set()
    def add(p, q):
        u, ou = cellsite(*p); v, ov = cellsite(*q)
        off = (ov[0] - ou[0], ov[1] - ou[1])
        # normalise so that u's cell is the origin
        key = (u, v, off); rk = (v, u, (-off[0], -off[1]))
        if rk not in edges: edges.add(key)
    for (x, y) in sites:
        add((x, y), (x + 1, y)); add((x, y), (x, y + 1))
    add((0, 0), (1, 1)); add((0, 1), (-1, 2))     # dimers: (0,0)-(1,1) and (2,1)-(1,2) ~ (0,1)-(-1,2)
    return 4, sorted(edges)
def elongated_triangular():   # 3^3.4^2
    return geometric((1, 0), (0.5, 1 + s3 / 2), [(0, 0), (0, 1)])
def truncated_square():       # 4.8.8
    a = 1 + s2; r = s2 / 2
    return geometric((a, 0), (0, a), [(r, 0), (0, r), (-r, 0), (0, -r)])
def ruby():                   # 3.4.6.4 (rhombitrihexagonal)
    a = 1 + s3
    basis = [(math.cos(math.radians(30 + 60 * k)), math.sin(math.radians(30 + 60 * k))) for k in range(6)]
    return geometric((a, 0), (a / 2, a * s3 / 2), basis)
def star():                   # 3.12.12 (truncated hexagonal)
    h = 1 + 2 / s3; r = 1 / s3
    a1 = (s3 * h, 0); a2 = (s3 * h / 2, 1.5 * h)
    A = np.array([0, 0.0]); Bc = np.array([0, h])
    dirsA = [(0, 1), (-s3 / 2, -0.5), (s3 / 2, -0.5)]
    basis = [tuple(A + r * np.array(d)) for d in dirsA] + [tuple(Bc - r * np.array(d)) for d in dirsA]
    return geometric(a1, a2, basis)
def truncated_trihexagonal(): # 4.6.12
    h = 1 + s3
    a1 = (s3 * h, 0); a2 = (s3 * h / 2, 1.5 * h)
    A = np.array([0, 0.0]); Bc = np.array([0, h])
    basis = []
    for c, sgn in ((A, 1), (Bc, -1)):
        for ang in (90, 210, 330):
            for dphi in (-30, 30):
                t = math.radians(ang + dphi + (0 if sgn == 1 else 180))
                basis.append(tuple(c + np.array([math.cos(t), math.sin(t)])))
    return geometric(a1, a2, basis)
def snub_hexagonal():         # 3^4.6: triangular lattice minus a sqrt7 x sqrt7 superlattice
    A1, A2 = (3, -1), (1, 2)
    det = A1[0] * A2[1] - A1[1] * A2[0]     # 7
    def red(x, y):   # reduce (x,y) mod superlattice -> (rep, cell offset)
        # solve (x,y) = n1 A1 + n2 A2 + r with r in fundamental set
        n1 = math.floor((x * A2[1] - y * A2[0]) / det); n2 = math.floor((A1[0] * y - A1[1] * x) / det)
        return (x - n1 * A1[0] - n2 * A2[0], y - n1 * A1[1] - n2 * A2[1]), (n1, n2)
    reps = sorted(set(red(x, y)[0] for x in range(-6, 7) for y in range(-6, 7)))
    hole = red(0, 0)[0]
    sites = [r for r in reps if r != hole]
    idx = {s: i for i, s in enumerate(sites)}
    nbr = [(1, 0), (0, 1), (1, -1), (-1, 0), (0, -1), (-1, 1)]
    edges = set()
    for s in sites:
        for d in nbr:
            r, off = red(s[0] + d[0], s[1] + d[1])
            if r == hole: continue
            key = (idx[s], idx[r], off); rk = (idx[r], idx[s], (-off[0], -off[1]))
            if rk not in edges and key not in edges: edges.add(key)
    return len(sites), sorted(edges)

LATTICES = {"lieb": lieb, "dice": dice, "checkerboard": checkerboard, "shastry_sutherland(3^2.4.3.4)": shastry_sutherland,
            "elongated_triangular(3^3.4^2)": elongated_triangular, "truncated_square(4.8.8)": truncated_square,
            "ruby(3.4.6.4)": ruby, "star(3.12.12)": star, "truncated_trihexagonal(4.6.12)": truncated_trihexagonal,
            "snub_hexagonal(3^4.6)": snub_hexagonal}

if __name__ == "__main__":
    for name, f in LATTICES.items():
        nv, E = f()
        deg = [0] * nv
        for (u, v, o) in E: deg[u] += 1; deg[v] += 1
        print(f"{name:34s} nv={nv:2d} edges={len(E):2d} degrees={deg}")

def line_graph_of_bravais(dirs):
    """line graph of a 2D Bravais lattice with nearest-neighbour vectors dirs (one per pair)."""
    P = len(dirs)
    def ends(t, x, y):
        d = dirs[t]; return {(x, y), (x + d[0], y + d[1])}
    edges = set()
    for t0 in range(P):
        for t in range(P):
            for x in range(-3, 4):
                for y in range(-3, 4):
                    if (t, x, y) == (t0, 0, 0): continue
                    if ends(t, x, y) & ends(t0, 0, 0):
                        key = (t0, t, (x, y)); rk = (t, t0, (-x, -y))
                        if rk not in edges: edges.add(key)
    return P, sorted(edges)
LATTICES["linegraph_triangular"] = lambda: line_graph_of_bravais([(1, 0), (0, 1), (1, -1)])

def union_jack():      # tetrakis square: corners (0) deg 8, face centres (1) deg 4
    E = [(0, 0, (1, 0)), (0, 0, (0, 1)),
         (0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1)), (0, 1, (-1, -1))]
    return 2, E
def triakis_triangular():   # triangular lattice (0) + centroids of up (1) and down (2) triangles
    E = [(0, 0, (1, 0)), (0, 0, (0, 1)), (0, 0, (1, -1)),
         (0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1)),
         (0, 2, (-1, 0)), (0, 2, (0, -1)), (0, 2, (-1, -1))]
    return 3, E
LATTICES["union_jack"] = union_jack
LATTICES["triakis_triangular"] = triakis_triangular
# Laves tilings on the triangular lattice. Cell: 0 lattice point, 1-3 edge midpoints (edges e1, e2, e1-e2),
# 4 up-triangle centre (R, R+e1, R+e2), 5 down-triangle centre (R, R+e1, R+e1-e2)
def _deltoidal_edges():
    return [(0, 1, (0, 0)), (1, 0, (1, 0)), (0, 2, (0, 0)), (2, 0, (0, 1)), (0, 3, (0, 0)), (3, 0, (1, -1)),
            (4, 1, (0, 0)), (4, 2, (0, 0)), (4, 3, (0, 1)), (5, 1, (0, 0)), (5, 3, (0, 0)), (5, 2, (1, -1))]
def deltoidal_trihexagonal():   # dual of 3.4.6.4: degrees 6 (lattice), 4 (midpoints), 3 (centres)
    return 6, _deltoidal_edges()
def kisrhombille():              # dual of 4.6.12 (barycentric subdivision of triangular): degrees 12, 4, 6
    return 6, _deltoidal_edges() + [(4, 0, (0, 0)), (4, 0, (1, 0)), (4, 0, (0, 1)), (5, 0, (0, 0)), (5, 0, (1, 0)), (5, 0, (1, -1))]
LATTICES["deltoidal_trihexagonal"] = deltoidal_trihexagonal
LATTICES["kisrhombille"] = kisrhombille
def king():   # square lattice + both diagonals (degree 8), one site per cell
    return 1, [(0, 0, (1, 0)), (0, 0, (0, 1)), (0, 0, (1, 1)), (0, 0, (1, -1))]
LATTICES["king"] = king
def honeycomb_nnn():   # honeycomb + all next-nearest neighbours (unit weights), degree 9
    return 2, [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))] + [(s, s, o) for s in (0, 1) for o in [(1, 0), (0, 1), (1, -1)]]
LATTICES["honeycomb_nnn"] = honeycomb_nnn
LATTICES["honeycomb"] = lambda: (2, [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))])
