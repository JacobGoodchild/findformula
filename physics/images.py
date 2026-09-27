"""Potential-flow added mass of a sphere (radius a) touching the INSIDE of a fixed
spherical container (radius b), moving along the line of centres. Method of images
(axial doublets), contact point at x=0, both spheres on the x>0 side."""
import sympy as sp
a, b, n = sp.symbols('a b n', positive=True)
d = 1/a - 1/b
xin = lambda k: 1 / (1/a + k*d)        # inner-sphere images (k=0 is the sphere centre)
xex = lambda k: -1 / (k*d)             # images outside the container (k>=1)
# strength recursion: A_k = A_{k-1} * (b/(b - x_{k-1}))^3 * (a/(a - x'_k))^3
fac = sp.simplify((b/(b - xin(n-1)))**3 * (a/(a - xex(n)))**3)
print("ratio A_n/A_{n-1} =", sp.factor(fac))

k = sp.symbols('k', positive=True, integer=True)
U, rho = sp.symbols('U rho', positive=True)
A0 = -U*a**3/2
beta = b - a
A = lambda m: A0 * b**3 / (a + (m+1)*beta)**3           # inner images, m>=0
B = lambda m: -A(m-1) * (b/(b - xin(m-1)))**3             # images outside container, m>=1
inner_term = sp.simplify(A(k-1))
ext_term = sp.simplify(-2*B(k)/(a - xex(k))**3)
print("inner image strength A_(k-1):", sp.factor(inner_term))
print("ext contribution d(phi)/dz at centre, term k:", sp.factor(ext_term))
