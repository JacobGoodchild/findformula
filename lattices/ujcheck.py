import numpy as np, pickle, mpmath as mp
from builders import LATTICES
from gen import out_arcs
nv, E = LATTICES["union_jack"]()
d = pickle.load(open("full_union_jack_0.pkl", "rb"))
arcs = d["arcs"]
N = 128
ks = (np.arange(N) + 0.5) * 2 * np.pi / N
YQ = np.zeros((len(arcs), len(arcs)))
for k1 in ks:
    for k2 in ks:
        A = np.zeros((nv, nv), complex); deg = np.zeros(nv)
        for (u, v, o) in E:
            ph = np.exp(1j * (k1 * o[0] + k2 * o[1]))
            A[u, v] += ph; A[v, u] += np.conj(ph); deg[u] += 1; deg[v] += 1
        Qi = np.linalg.inv(np.diag(deg) + A)
        def G(p, q): return Qi[p[0], q[0]] * np.exp(1j * (k1 * (p[1][0] - q[1][0]) + k2 * (p[1][1] - q[1][1])))
        O = (0, (0, 0))
        for i, hi in enumerate(arcs):
            for j, hj in enumerate(arcs):
                YQ[i, j] += np.real(G(O, O) + G(O, hj) + G(hi, O) + G(hi, hj))
YQ /= N * N
eng = np.array([[float(x) for x in r] for r in d["YQ"]])
print("max |grid - engine| over YQ entries:", np.abs(YQ - eng).max())
print(np.round(YQ - eng, 9))
