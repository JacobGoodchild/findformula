import mpmath as mp
from fcc import lgf
mp.mp.dps = 50
pi = mp.pi; s2 = mp.sqrt(2); s3 = mp.sqrt(3)
Wsc = mp.sqrt(6)/(32*pi**3)*mp.gamma(mp.mpf(1)/24)*mp.gamma(mp.mpf(5)/24)*mp.gamma(mp.mpf(7)/24)*mp.gamma(mp.mpf(11)/24)
Y = Wsc * s2
v200 = lgf((2,0,0), -1); v220 = lgf((2,2,0), -1)
print(v200, v220)
names = ["v", "1", "Y", "1/(pi^2 Y)", "1/pi", "s2", "s3", "s2*Y", "s3*Y", "s2/(pi^2Y)", "s3/(pi^2 Y)", "s3/pi", "s6/pi"]
for v in (v200, v220):
    basis = [v, 1, Y, 1/(pi**2*Y), 1/pi, s2, s3, s2*Y, s3*Y, s2/(pi**2*Y), s3/(pi**2*Y), s3/pi, mp.sqrt(6)/pi]
    r = mp.pslq(basis, maxcoeff=2000, maxsteps=10**6)
    print([(n, c) for n, c in zip(names, r) if c] if r else None)
