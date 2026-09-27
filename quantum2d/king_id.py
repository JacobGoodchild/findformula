import mpmath as mp
from lgf2d import a
mp.mp.dps = 40; pi = mp.pi; s = mp.sqrt
A = lambda k: 4 - mp.cos(k); C = lambda k: 1 + 2*mp.cos(k)
v = a((1,0), A, C)
names = ["v","1","1/pi","s2/pi","s3/pi","s5/pi","s7/pi","s15/pi","atan(s2)/pi","atan(s7)/pi","atan(1/s3)/pi","atan(s15)/pi","s2 atan(s2)/pi", "s7 atan(s7)/pi","s3 atan(s3/5)/pi"]
vals = [1,1/pi,s(2)/pi,s(3)/pi,s(5)/pi,s(7)/pi,s(15)/pi,mp.atan(s(2))/pi,mp.atan(s(7))/pi,mp.atan(1/s(3))/pi,mp.atan(s(15))/pi,s(2)*mp.atan(s(2))/pi,s(7)*mp.atan(s(7))/pi,s(3)*mp.atan(s(3)/5)/pi]
r = mp.pslq([v]+vals, maxcoeff=5000, maxsteps=10**6)
print(v, [(n,c) for n,c in zip(names,r) if c] if r else None)
# elliptic: king LGF may involve K at modulus from 4 - c1 - c2 - 2c1c2; try K and E at various k^2
for m in [mp.mpf(1)/2, mp.mpf(1)/4, mp.mpf(1)/9, mp.mpf(8)/9, mp.mpf(3)/4, mp.mpf(1)/16, mp.mpf(15)/16, (2-s(3))/4, (s(2)-1)**2, mp.mpf(1)/3, mp.mpf(2)/3, mp.mpf(1)/5, mp.mpf(4)/5, mp.mpf(1)/8, mp.mpf(7)/8]:
    K, E = mp.ellipk(m), mp.ellipe(m)
    r = mp.pslq([v, 1, K/pi, E/pi, pi/K, 1/pi, 1/(pi*K)], maxcoeff=2000, maxsteps=10**5)
    if r: print("m=", m, r)
