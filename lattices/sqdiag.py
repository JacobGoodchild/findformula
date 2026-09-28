"""One-site lattice with bonds e1 (w1), e2 (w2), e1+e2 (wp), e1-e2 (wm): resistances by exact k2-average + quad."""
import mpmath as mp
def R(w1, w2, wp, wm, e, dps=30):
    mp.mp.dps = dps
    w1, w2, wp, wm = map(mp.mpf, (w1, w2, wp, wm))
    def num(k1, k2):
        return {"1": 4*mp.sin((k1)/2)**2, "2": 4*mp.sin((k2)/2)**2, "p": 4*mp.sin((k1 + k2)/2)**2, "m": 4*mp.sin((k1 - k2)/2)**2}[e]
    def D(k1, k2):
        return w1*(4*mp.sin((k1)/2)**2) + w2*(4*mp.sin((k2)/2)**2) + wp*(4*mp.sin((k1 + k2)/2)**2) + wm*(4*mp.sin((k1 - k2)/2)**2)
    f = lambda k1, k2: num(k1, k2) / D(k1, k2) if (k1 or k2) else 0
    return mp.quad(f, [-mp.pi, 0, mp.pi], [-mp.pi, 0, mp.pi]) / (4 * mp.pi**2)
def R1_formula(w1, w2, wp, wm):
    s = wp + wm; Dl = 4*wp*wm - w1*w1
    x = Dl / ((w1 + s) * (w1 + w2))
    if Dl > 0: return 2 / (mp.pi * mp.sqrt(Dl)) * mp.asinh(mp.sqrt(x))
    if Dl < 0: return 2 / (mp.pi * mp.sqrt(-Dl)) * mp.asin(mp.sqrt(-x))
    return 2 / (mp.pi * mp.sqrt((w1 + s) * (w1 + w2)))
