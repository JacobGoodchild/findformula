import mpmath as mp
pi = mp.pi
def lin(v, basis, maxc=10**5):
    names = list(basis); vals = [basis[n] for n in names]
    r = mp.pslq([v] + vals, maxcoeff=maxc, maxsteps=10**6)
    return [(n, c) for n, c in zip(["v"] + names, r) if c] if r else None
def B_sq(): return {"1": 1, "1/pi": 1/pi, "1/pi^2": 1/pi**2}
def B_tri(): s3 = mp.sqrt(3); return {"1": 1, "s3/pi": s3/pi, "1/pi^2": 1/pi**2}
def B_mix(): s2, s3 = mp.sqrt(2), mp.sqrt(3); return {"1": 1, "1/pi": 1/pi, "s2/pi": s2/pi, "s3/pi": s3/pi, "1/pi^2": 1/pi**2, "s2/pi^2": s2/pi**2, "s3/pi^2": s3/pi**2, "s6/pi^2": mp.sqrt(6)/pi**2}
