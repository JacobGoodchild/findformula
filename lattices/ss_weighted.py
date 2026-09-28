"""Shastry-Sutherland / snub-square lattice with square-edge conductance 1 and dimer conductance w:
tree entropy E_k ln det L(k),  det L = -4w(w+2)s^2 + 4(w+1)^2 t^2 - 16(w+2)^2 s + 16(w+2)(3w+4),
s = cos k1 + cos k2, t = cos k1 - cos k2."""
import mpmath as mp, sys
def Elogdet(w, dps=50):
    mp.mp.dps = dps
    w = mp.mpf(w)
    def inner(x):
        # det as quadratic in y = cos k2:  s = x + y, t = x - y
        a = -4*w*(w+2) + 4*(w+1)**2
        b = -8*w*(w+2)*x - 8*(w+1)**2*x - 16*(w+2)**2
        c = -4*w*(w+2)*x*x + 4*(w+1)**2*x*x - 16*(w+2)**2*x + 16*(w+2)*(3*w+4)
        D = mp.sqrt(b*b - 4*a*c + 0j)
        r1 = (-b + D)/(2*a); r2 = (-b - D)/(2*a)
        L = lambda r: mp.log((r + mp.sqrt(r - 1)*mp.sqrt(r + 1))/2)    # E_y ln|y - r|, r outside [-1,1]
        return mp.log(abs(a)) + L(r1) + L(r2)
    return mp.re(mp.quad(lambda k: inner(mp.cos(k)), [0, mp.pi/4, mp.pi/2, 3*mp.pi/4, mp.pi]))/mp.pi
if __name__ == "__main__":
    for w in sys.argv[1:]:
        print(w, Elogdet(mp.mpf(w)))
