import mpmath as mp, pickle
from identify import find_lin
mp.mp.dps = 35; pi = mp.pi; s2 = mp.sqrt(2); s3 = mp.sqrt(3)
(hp, _, _), Yh, (kp, kpp, kpm), Yk, YQk = pickle.load(open("hk.pkl", "rb"))
m = mp.mpf(5)/32; K = mp.ellipk(m); E = mp.ellipe(m)
for j in range(4):
    print("Y0%d" % j, Yk[0][j], find_lin(Yk[0][j], ["1", "s3/pi"], [1, s3/pi]))
for j in range(4):
    print("YQ0%d" % j, YQk[0][j], find_lin(YQk[0][j], ["1", "s2K/pi", "s2E/pi", "s2pi/K"], [1, s2*K/pi, s2*E/pi, s2*pi/K]))
flat = mp.mpf('0.010473571335890761167381588113434982')
print("flat channel", find_lin(flat, ["1", "s3/pi", "1/pi^2"], [1, s3/pi, 1/pi**2]))
print("+1 channel", kpp, find_lin(kpp, ["1", "s3/pi", "1/pi^2"], [1, s3/pi, 1/pi**2]))
