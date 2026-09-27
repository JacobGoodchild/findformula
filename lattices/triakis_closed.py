"""Triakis triangular (kisrhombille-dual 3.3.3.3.3.3 + centroids) flip-flop Grover trapping: closed forms vs engine."""
import pickle, mpmath as mp, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 25
pi = mp.pi; r = mp.sqrt(3) / pi
m = mp.mpf(1) / 2 - 37 / (28 * mp.sqrt(7))
g = mp.ellipk(m) / (mp.mpf(3087)**0.25 * pi)                      # E[1/(15 + c1 + c2 + c12)]
num = lambda t: mp.quad(lambda k1: 1 / mp.sqrt((t + mp.cos(k1))**2 - (2 * mp.cos(k1 / 2))**2), [0, pi]) / pi
tol = mp.mpf(10)**(-mp.mp.dps + 2)
with mp.workdps(2 * mp.mp.dps):
    assert abs(num(mp.mpf(15)) - g) < tol
    gp = mp.diff(num, 15)                                           # g'(15)
d0 = pickle.load(open("full_triakis_triangular_0.pkl", "rb")); d1 = pickle.load(open("full_triakis_triangular_1.pkl", "rb"))
# corner (degree 12) arcs: 0 +e1, 1 -e1, 2 +e2, 3 -e2, 4 +(e1-e2), 5 -(e1-e2), 6-8 up-centroids, 9-11 down-centroids
Y = [mp.mpf(1) / 5, (6 * r - 3) / 5, mp.mpf(1) / 10, (2 - 3 * r) / 5, mp.mpf(1) / 10, (2 - 3 * r) / 5,
     mp.mpf(1) / 10, (3 * r - 1) / 15, (5 - 6 * r) / 30, (5 - 6 * r) / 30, mp.mpf(1) / 10, (3 * r - 1) / 15]
q1 = -13 - 264 * g - 6804 * gp; q3 = 231 * g + 3402 * gp; q7 = 5 + 3 * g + 1134 * gp; q8 = (1 - 468 * g - 6804 * gp) / 6
Q = [1 - 12 * g, q1, (3 - 42 * g) / 2, q3, (3 - 42 * g) / 2, q3, 3 * g - mp.mpf(1) / 6, q7, q8, q8, 3 * g - mp.mpf(1) / 6, q7]
for i in range(12):
    assert abs(Y[i] - d0["Y"][0][i]) < 1e-20 and abs(Q[i] - d0["YQ"][0][i]) < 1e-20, i
flat = (4 - 4 * r + 3 * r**2) / 144                                   # mu = 0 flat band, centroid-bound starts
# rows for other starts via the symmetry of the pickled matrices: use engine rows for the diagonal start 6
ch = lambda row, s: mp.fsum(((1 if j == s else 0) - y)**2 for j, y in enumerate(row)) / 4
tot0 = ch(Y, 0) + ch(Q, 0)
print("corner, lattice-arc start:", mp.nstr(tot0, 25), " engine", mp.nstr(d0["flat"][0][3], 25))
print("   flat mu=0 channel (centroid start) engine", mp.nstr(d0["flat"][6][2]["0"], 25), " closed", mp.nstr(flat, 25))
gc = 27 / mp.mpf(200) + mp.mpf(1) / 24 + mp.mpf(3) / 8 * (1 - 6 * g)**2
print("centroid start:", mp.nstr(gc, 25), " engine", mp.nstr(d1["flat"][0][3], 25))
print("g =", mp.nstr(g, 30), " g' =", mp.nstr(gp, 30))
# start 6: corner -> up-centroid arc
t = mp.mpf(1)
Y6 = [t / 10, (3 * r - 1) / 15, t / 10, (3 * r - 1) / 15, (5 - 6 * r) / 30, (5 - 6 * r) / 30,
      2 * t / 5, t / 30, t / 30, (4 - 3 * r) / 45, (4 - 3 * r) / 45, (6 * r - 2) / 45]
a = 33 * g + 378 * gp - 4 * t / 9
Q6 = [3 * g - t / 6, q7, 3 * g - t / 6, q7, q8, q8, 6 * g, 33 * g - 13 * t / 6, 33 * g - 13 * t / 6, a, a, 6 * g - 756 * gp - 34 * t / 9]
for i in range(12):
    assert abs(Y6[i] - d0["Y"][6][i]) < 1e-20 and abs(Q6[i] - d0["YQ"][6][i]) < 1e-20, i
tot6 = ch(Y6, 6) + ch(Q6, 6) + flat
print("corner, centroid-arc start:", mp.nstr(tot6, 25), " engine", mp.nstr(d0["flat"][6][3], 25))
import sympy as sp
R, Gs, Gp = sp.symbols("r g gp")
def S(x): return x
Ys = [sp.Rational(1, 5), (6 * R - 3) / 5, sp.Rational(1, 10), (2 - 3 * R) / 5, sp.Rational(1, 10), (2 - 3 * R) / 5,
      sp.Rational(1, 10), (3 * R - 1) / 15, (5 - 6 * R) / 30, (5 - 6 * R) / 30, sp.Rational(1, 10), (3 * R - 1) / 15]
Q1 = -13 - 264 * Gs - 6804 * Gp; Q3 = 231 * Gs + 3402 * Gp; Q7 = 5 + 3 * Gs + 1134 * Gp; Q8 = (1 - 468 * Gs - 6804 * Gp) / 6
Qs = [1 - 12 * Gs, Q1, (3 - 42 * Gs) / 2, Q3, (3 - 42 * Gs) / 2, Q3, 3 * Gs - sp.Rational(1, 6), Q7, Q8, Q8, 3 * Gs - sp.Rational(1, 6), Q7]
A6 = 33 * Gs + 378 * Gp - sp.Rational(4, 9)
Ys6 = [sp.Rational(1, 10), (3 * R - 1) / 15] * 2 + [(5 - 6 * R) / 30] * 2 + [sp.Rational(2, 5), sp.Rational(1, 30), sp.Rational(1, 30), (4 - 3 * R) / 45, (4 - 3 * R) / 45, (6 * R - 2) / 45]
Qs6 = [3 * Gs - sp.Rational(1, 6), Q7] * 2 + [Q8, Q8, 6 * Gs, 33 * Gs - sp.Rational(13, 6), 33 * Gs - sp.Rational(13, 6), A6, A6, 6 * Gs - 756 * Gp - sp.Rational(34, 9)]
chs = lambda row, s: sp.expand(sum(((1 if j == s else 0) - y)**2 for j, y in enumerate(row)) / 4)
print("lattice-arc start  +1:", chs(Ys, 0), "\n                   -1:", chs(Qs, 0))
print("centroid-arc start +1:", chs(Ys6, 6), "\n                   -1:", chs(Qs6, 6), "\n                 flat:", sp.expand((4 - 4 * R + 3 * R**2) / 144))
