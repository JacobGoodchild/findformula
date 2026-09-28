"""Flip-flop Grover walk on the king lattice: closed-form transfer currents vs the engine."""
import pickle, mpmath as mp, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 28
pkl = sys.argv[2] if len(sys.argv) > 2 else "king_d30.pkl"
pi = mp.pi; s2 = mp.sqrt(2); m = mp.mpf(7) / 8
L = mp.log(2 + mp.sqrt(3)) / (pi * mp.sqrt(3))
ka = s2 * mp.ellipk(m) / pi; ep = s2 * mp.ellipe(m) / pi; rh = s2 * mp.ellippi(mp.mpf(3) / 4, m) / pi
# arcs: 0 +x, 1 -x, 2 +y, 3 -y, 4 (+1,+1), 5 (-1,-1), 6 (+1,-1), 7 (-1,+1)
a, b, c, e = L, 3 * L - 2 / pi, (6 * L - 1) / 4, (1 - 2 * L) / 4
f = (1 + 2 / pi - 6 * L) / 2
Yax = [a, b, c, c, e, f, e, f]
Ydg = [e, f, e, f, (1 - 2 * L) / 2, 6 * L + 2 / pi - 2, (1 - 4 / pi + 2 * L) / 2, (1 - 4 / pi + 2 * L) / 2]
q0 = (rh - ka) / 6; q1 = (ka - 24 * ep + 5 * rh) / 18; q2 = (3 - 4 * ka + rh) / 12; q5 = (13 * ka + 12 * ep - 7 * rh) / 18
Qax = [q0, q1, q2, q2, q2, q5, q2, q5]
Qdg = [q2, q5, q2, q5, (3 + ka - rh) / 6, (-18 - 31 * ka + 24 * ep + 16 * rh) / 18, (9 + 7 * ka - 24 * ep - rh) / 18, (9 + 7 * ka - 24 * ep - rh) / 18]
d = pickle.load(open(pkl, "rb"))
err = 0
for r, Y, Q in [(0, Yax, Qax), (4, Ydg, Qdg)]:
    for j in range(8):
        err = max(err, abs(Y[j] - d["Y"][r][j]), abs(Q[j] - d["YQ"][r][j]))
ch = lambda row, s: mp.fsum(((1 if j == s else 0) - y)**2 for j, y in enumerate(row)) / 4
print("max entry error", mp.nstr(err, 3), " row sums", mp.nstr(mp.fsum(Yax), 5), mp.nstr(mp.fsum(Qax), 5), mp.nstr(mp.fsum(Ydg), 5), mp.nstr(mp.fsum(Qdg), 5))
print("axis start:     pbar =", mp.nstr(ch(Yax, 0) + ch(Qax, 0), mp.mp.dps - 3), " engine", mp.nstr(d["flat"][0][3], 22))
print("diagonal start: pbar =", mp.nstr(ch(Ydg, 4) + ch(Qdg, 4), mp.mp.dps - 3), " engine", mp.nstr(d["flat"][4][3], 22))
print("kappa, eps, rho =", mp.nstr(ka, 20), mp.nstr(ep, 20), mp.nstr(rh, 20), " L =", mp.nstr(L, 20))
