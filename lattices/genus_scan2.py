"""Fast genus screen with a generic numeric weight: odd-multiplicity z1-factors of the z2-discriminant
of the Laplacian determinant (after removing powers of z1 and the node factor (z1-1))."""
import sympy as sp, sys
from gen import z1, z2
W = sp.Rational(3, 7)
def lapdet(nv, E):
    A = sp.zeros(nv, nv); deg = [0] * nv
    for (u, v, o, c) in E:
        ph = z1**o[0] * z2**o[1]
        A[u, v] += c * ph; A[v, u] += c / ph; deg[u] += c; deg[v] += c
    return sp.expand((sp.diag(*deg) - A).det())
def odd_deg(d):
    pw = [t.as_powers_dict().get(z2, 0) for t in sp.Add.make_args(d)]
    P = sp.expand(d * z2**(-min(pw)))
    pw1 = [t.as_powers_dict().get(z1, 0) for t in sp.Add.make_args(P)]
    P = sp.expand(P * z1**(-min(pw1)))
    disc = sp.factor_list(sp.discriminant(P, z2))[1]
    odd = [f for f, m in disc if m % 2 and f.has(z1) and f != z1]
    return sp.Poly(P, z2).degree(), [sp.Poly(f, z1).degree() for f in odd]
FAM = {
 "square+axisNNN (1 site)": (1, [(0,0,(1,0),1),(0,0,(0,1),1),(0,0,(2,0),W),(0,0,(0,2),W)]),
 "triangular+NNN sqrt3 (1 site)": (1, [(0,0,(1,0),1),(0,0,(0,1),1),(0,0,(1,-1),1),(0,0,(1,1),W),(0,0,(2,-1),W),(0,0,(1,-2),W)]),
 "honeycomb+NNN (2 sites)": (2, [(0,1,(0,0),1),(0,1,(-1,0),1),(0,1,(0,-1),1)] + [(s,s,o,W) for s in (0,1) for o in [(1,0),(0,1),(1,-1)]]),
 "union jack, centre weight w": (2, [(0,0,(1,0),1),(0,0,(0,1),1),(0,1,(0,0),W),(0,1,(-1,0),W),(0,1,(0,-1),W),(0,1,(-1,-1),W)]),
 "kagome, one triangle orientation weight w": (3, [(0,1,(0,0),1),(0,2,(0,0),1),(1,2,(0,0),1),(1,0,(1,0),W),(2,0,(0,1),W),(1,2,(1,-1),W)]),
 "lieb + diagonal? (square with decorated edges, weight w on half-edges)": (3, [(0,1,(0,0),1),(1,0,(1,0),W),(0,2,(0,0),1),(2,0,(0,1),W)]),
 "anisotropic triangular (control)": (1, [(0,0,(1,0),1),(0,0,(0,1),W),(0,0,(1,-1),2*W)]),
 "square + one diagonal family, alt. rows (2 sites)": (2, [(0,1,(0,0),1),(1,0,(1,0),1),(0,0,(0,1),1),(1,1,(0,1),1),(0,1,(0,1),W)]),
}
if __name__ == "__main__":
    for n, (nv, E) in FAM.items():
        try:
            print(n, odd_deg(lapdet(nv, E)), flush=True)
        except Exception as e:
            print(n, "error", e, flush=True)
