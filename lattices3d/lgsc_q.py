import mpmath as mp, sympy as sp, pickle, sys
sys.argv = ["x", "3", "45"]
src = open("lg_general.py").read().split("deg = 4 * P - 2")[0]
exec(src)
gs, res = pickle.load(open("sc_t7.pkl", "rb"))
g0, g1, g2 = sp.symbols('g0 g1 g2')
sym = {}
for x, (v, r) in res.items():
    sym[tuple(sorted(x))] = -(r[1] + r[2]*g0 + r[3]*g1 + r[4]*g2) / sp.Integer(r[0])
def Gs(x):
    key = tuple(sorted(abs(v) for v in x))
    return sym[key]
def BdBs(i, Ri, j, Rj):
    x = add(Rj, neg(Ri))
    return sum(Gs(add(add(x, a_), b_)) for a_, b_ in [(Z, Z), (neg(E[i]), Z), (Z, E[j]), (neg(E[i]), E[j])])
def Qs(i, Ri, j, Rj): return ((1 if (i, Ri) == (j, Rj) else 0) - BdBs(i, Ri, j, Rj)) / sp.Integer(4 * P - 4)
num = {g0: sp.Float(str(gs[0]), 45), g1: sp.Float(str(gs[1]), 45), g2: sp.Float(str(gs[2]), 45)}
for a in [nbrs[0], nbrs[2]]:
    Qr = [Qs(*o, *o) + Qs(*o, *b) + Qs(*a, *o) + Qs(*a, *b) for b in nbrs]
    pm = sp.expand(sum(((1 if b == a else 0) - q)**2 for b, q in zip(nbrs, Qr)) / 4)
    print("start", a, "\n  -1 =", pm, "\n     =", sp.N(pm.subs(num), 30))
