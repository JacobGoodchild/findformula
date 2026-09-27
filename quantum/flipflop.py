"""3D Grover walk with FLIP-FLOP shift: trapped (flat-band) part.
M = I/2 - N, with N built from simple-cubic lattice Green function values:
  N_aa = 1/(2d), N_{+i,-i} = alpha, N_{(+-)i,(+-)j} = beta (i != j),
  alpha = (1/4)[G(0) - 2G(e1) + G(2e1)],  beta = (1/4)[G(0) - 2G(e1) + G(e1+e2)],
  G(x) = E[cos(k.x) / sum_l (1 - cos k_l)] = int_0^inf e^{-3t} prod_i I_{x_i}(t) dt."""
import mpmath as mp
mp.mp.dps = 70
I = lambda n, t: mp.besseli(n, t)
f_alpha = lambda t: mp.exp(-3 * t) * I(0, t)**2 * (I(0, t) - 2 * I(1, t) + I(2, t)) / 4
f_beta = lambda t: mp.exp(-3 * t) * I(0, t) * (I(0, t) - I(1, t))**2 / 4
pts = [0, 1, 4, 16, 64, 256, 1024, mp.inf]
alpha = mp.quad(f_alpha, pts); beta = mp.quad(f_beta, pts)
W = mp.sqrt(6) / (32 * mp.pi**3) * mp.gamma(mp.mpf(1)/24) * mp.gamma(mp.mpf(5)/24) * mp.gamma(mp.mpf(7)/24) * mp.gamma(mp.mpf(11)/24)
print("alpha", alpha); print("beta", beta); print("Watson W = G(0)", W)
for name, v in [("alpha", alpha), ("beta", beta)]:
    print(name, mp.pslq([v, 1, W, 1 / (mp.pi**2 * W)], maxcoeff=10**6, maxsteps=10**6))
