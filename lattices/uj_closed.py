"""Union Jack (tetrakis square) flip-flop Grover trapping: closed forms vs the engine."""
import pickle, mpmath as mp, sys
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 25
pi = mp.pi
c = mp.sqrt(2) * mp.acos(mp.mpf(1) / 3) / pi          # Laplacian constant
K = mp.ellipk(mp.mpf(1) / 2); k = K / pi; l = 1 / K    # lemniscatic K(1/sqrt2)
d0 = pickle.load(open(sys.argv[2] if len(sys.argv) > 2 else "uj0_d25.pkl", "rb"))
d1 = pickle.load(open("uj1_d25.pkl", "rb"))
# PSLQ of Laplacian row 4 (NE centre arc from a corner)
for i in range(8):
    print("Y[4][%d]" % i, mp.nstr(d0["Y"][4][i], 15), mp.pslq([d0["Y"][4][i], 1, 1 / pi, c], maxcoeff=10**4, maxsteps=10**6))

# ---- closed forms ----
# corner origin, arcs: 0 +x, 1 -x, 2 +y, 3 -y, 4 NE, 5 NW, 6 SE, 7 SW (centres)
Yx = [c / 2, 2 / pi - c, 2 * c - 1, 2 * c - 1, c / 4, mp.mpf(3) / 2 - 1 / pi - 2 * c, c / 4, mp.mpf(3) / 2 - 1 / pi - 2 * c]
Qx = [l / 2, (4 * k - 3 * l) / 6, k + l - 1, k + l - 1, (2 - 2 * k - l) / 4, (12 - 10 * k - 9 * l) / 12, (2 - 2 * k - l) / 4, (12 - 10 * k - 9 * l) / 12]
a, b = c / 4, mp.mpf(3) / 2 - 1 / pi - 2 * c
Yd = [a, b, a, b, (2 - c) / 4, 1 - 1 / (2 * pi) - 11 * c / 8, 1 - 1 / (2 * pi) - 11 * c / 8, (13 * c - 9 + 6 / pi) / 2]
a, b = (2 - 2 * k - l) / 4, (12 - 10 * k - 9 * l) / 12
Qd = [a, b, a, b, (2 - l) / 4, k / 6, k / 6, (28 * k + 27 * l - 30) / 12]
# centre origin (degree 4), arcs 0..3 to the four corners (0 and 3 opposite)
Yc = [(2 - c) / 4, (1 - c) / 2, (1 - c) / 2, (5 * c - 2) / 4]
Qc = [(2 - l) / 4, (1 - k) / 2, (1 - k) / 2, (4 * k + l - 2) / 4]
ch = lambda row, s: mp.fsum(((1 if j == s else 0) - y)**2 for j, y in enumerate(row)) / 4
for name, Y, Q, s, dd, r in [("corner, axis start", Yx, Qx, 0, d0, 0), ("corner, diagonal start", Yd, Qd, 4, d0, 4), ("centre start", Yc, Qc, 0, d1, 0)]:
    err = max(max(abs(Y[j] - dd["Y"][r][j]) for j in range(len(Y))), max(abs(Q[j] - dd["YQ"][r][j]) for j in range(len(Y))))
    tot = ch(Y, s) + ch(Q, s)
    print(f"{name}: p = {mp.nstr(tot, 30)}  engine {mp.nstr(dd['flat'][s][3], 22)}  max entry err {mp.nstr(err, 3)}  rowsums {mp.nstr(mp.fsum(Y), 5)}")
