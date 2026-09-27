import mpmath as mp
mp.mp.dps = 50
for d in (2, 3, 4):
    A = mp.quad(lambda s: mp.erfc(mp.sqrt(s))**d, [0, 1, 5, 20, mp.inf])
    print(d, A)
    basis = [A, 1, 1/mp.pi, mp.sqrt(2)/mp.pi, mp.sqrt(3)/mp.pi, mp.atan(mp.sqrt(2))/mp.pi, mp.atan(1/mp.sqrt(8))/mp.pi, mp.sqrt(2), mp.sqrt(3), mp.sqrt(6)/mp.pi]
    print("  pslq:", mp.pslq(basis, maxcoeff=10**6, maxsteps=10**6))
