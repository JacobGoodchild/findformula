"""King lattice (square + both diagonals, diagonal conductance w): resistances and tree entropy.
D(k) = 2(1-c1) + 2(1-c2) + 4w(1 - c1 c2)."""
import mpmath as mp
def Rdiag(w, dps=40):
    """Resistance across one diagonal bond = E[2(1 - c1 c2)/D] (both diagonals are equivalent)."""
    mp.mp.dps = dps; w = mp.mpf(w)
    def inner(k):
        x = mp.cos(k); om = 2*mp.sin(k/2)**2          # 1 - x without cancellation
        A = 2 + 2*om + 4*w; B = 2 + 4*w*x
        if abs(B) < mp.mpf('1e-6'):                   # removable point B = 0: integrate over k2 directly
            return mp.quad(lambda q: 2*(1 - x*mp.cos(q))/(A - B*mp.cos(q)), [0, mp.pi])/mp.pi
        AmB = 2*om*(1 + 2*w)
        J = 1/mp.sqrt(AmB*(A + B))
        return (2*x + 4*om*om*J)/B                    # = E_y[2(1 - x y)/(A - B y)]
    return mp.quad(inner, [0, mp.pi/4, mp.pi/2, 3*mp.pi/4, mp.pi])/mp.pi
def Elogdet(w, dps=40):
    mp.mp.dps = dps; w = mp.mpf(w)
    def inner(x):
        A = 4 - 2*x + 4*w; B = 2 + 4*w*x
        return mp.log((A + mp.sqrt(A*A - B*B))/2)
    return mp.quad(lambda k: inner(mp.cos(k)), [0, mp.pi/2, mp.pi])/mp.pi
if __name__ == "__main__":
    import sys
    for w in sys.argv[1:]:
        print(w, Rdiag(w), Elogdet(w))
