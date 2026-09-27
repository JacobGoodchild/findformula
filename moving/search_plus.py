import mpmath as mp, pickle, itertools
mp.mp.dps = 32; pi = mp.pi; s3 = mp.sqrt(3); s2 = mp.sqrt(2)
V = pickle.load(open('tri_moving_plus32.pkl', 'rb'))
C = {'s3/pi': s3/pi, '1/pi': 1/pi, '1/pi^2': 1/pi**2, 's3/pi^2': s3/pi**2, 'log2/pi': mp.log(2)/pi, 'log3/pi': mp.log(3)/pi,
     's3log2/pi': s3*mp.log(2)/pi, 's3log3/pi': s3*mp.log(3)/pi, 'log2/pi^2': mp.log(2)/pi**2, 'log3/pi^2': mp.log(3)/pi**2,
     'Cl/pi^2': mp.clsin(2, pi/3)/pi**2, 's3Cl/pi^2': s3*mp.clsin(2, pi/3)/pi**2, 'G/pi^2': mp.catalan/pi**2, 's2/pi': s2/pi,
     's6/pi': mp.sqrt(6)/pi, 'acos13/pi': mp.acos(mp.mpf(1)/3)/pi, 'log(2+s3)/pi': mp.log(2+s3)/pi, 's3log(2+s3)/pi': s3*mp.log(2+s3)/pi,
     'zeta3/pi^2': mp.zeta(3)/pi**2, 'log2^2/pi^2': mp.log(2)**2/pi**2}
names = list(C)
v = V[(1, 1)]
hits = 0
for r_ in (2, 3, 4):
    for combo in itertools.combinations(names, r_):
        r = mp.pslq([v, 1] + [C[c] for c in combo], maxcoeff=600, maxsteps=20000)
        if r and r[0] != 0:
            print(r_, [(n, c) for n, c in zip(['v', '1'] + list(combo), r) if c], flush=True); hits += 1
            if hits > 5: raise SystemExit
print("done, hits", hits)
