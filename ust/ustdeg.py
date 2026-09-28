"""Exact degree distribution of the uniform spanning tree at a vertex, from the transfer-current matrix Y
of the star (Burton-Pemantle): P(exactly S present) = det K, K_ij = Y_ij (i in S), delta_ij - Y_ij (i not in S)."""
import sympy as sp, itertools
def degree_distribution(Y):
    """Generating function sum_k p_k x^k = det(I + (x - 1) Y) (multilinearity in the rows)."""
    from sympy.polys.matrices import DomainMatrix
    z = Y.shape[0]; x = sp.Symbol('x')
    M = sp.eye(z) + (x - 1) * Y
    gens = sorted(M.free_symbols, key=str)
    dM = DomainMatrix.from_Matrix(M).convert_to(sp.QQ[tuple(gens)])
    Phi = sp.Poly(dM.domain.to_sympy(dM.det()), x)
    return [sp.factor(Phi.coeff_monomial(x**k)) for k in range(z + 1)]
if __name__ == "__main__":
    q = sp.symbols('q')        # q = 1/pi
    # square lattice: Y(+x,+x) = 1/2, Y(+x,-x) = 2/pi - 1/2 (Formula 7 check), Y(+x,+y) = 1/2 - 1/pi
    a, b, c = sp.Rational(1, 2), 2 * q - sp.Rational(1, 2), sp.Rational(1, 2) - q
    arcs = ['+x', '-x', '+y', '-y']
    def Ysq(i, j):
        if i == j: return a
        if {arcs[i], arcs[j]} in ({'+x', '-x'}, {'+y', '-y'}): return b
        return c
    Y = sp.Matrix(4, 4, Ysq)
    p = degree_distribution(Y)
    for k, v in enumerate(p): print(k, v, sp.N(v.subs(q, 1 / sp.pi)))
