"""Generic flip-flop Grover walk simulator on a 2D periodic graph (unit cell with nv vertices).
State: amplitude on arcs. Arc (u -> v, offset o) at vertex u in cell R points to v in cell R+o."""
import numpy as np
def simulate(nv, edges, origin, T, start_arc=0):
    arcs = []
    for (u, v, o) in edges:
        arcs.append((u, v, o)); arcs.append((v, u, (-o[0], -o[1])))
    out = {x: [i for i, a in enumerate(arcs) if a[0] == x] for x in range(nv)}
    rev = [arcs.index((a[1], a[0], (-a[2][0], -a[2][1]))) for a in arcs]
    L = 2 * T + 5; c = L // 2
    psi = np.zeros((len(arcs), L, L), complex)
    first = out[origin][start_arc]
    psi[first, c, c] = 1
    rec = []
    for t in range(T):
        new = np.zeros_like(psi)
        for x in range(nv):          # Grover coin at each vertex over its out-arcs
            idx = out[x]; d = len(idx)
            s = psi[idx].sum(axis=0)
            for i in idx:
                new[i] = 2 * s / d - psi[i]
        psi2 = np.zeros_like(psi)
        for i, (u, v, o) in enumerate(arcs):   # flip-flop: move along arc, land on reverse arc
            psi2[rev[i]] = np.roll(new[i], o, axis=(0, 1))
        psi = psi2
        rec.append(sum(np.abs(psi[i, c, c])**2 for i in out[origin]))
    return np.array(rec), arcs[first]
if __name__ == "__main__":
    honey = [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))]
    r, a = simulate(2, honey, 0, 400)
    print("honeycomb sim avg:", r[150:].mean(), "(pred 1/12 = 0.083333)", a)
    kag = [(0, 1, (0, 0)), (0, 2, (0, 0)), (1, 2, (0, 0)), (1, 0, (1, 0)), (2, 0, (0, 1)), (1, 2, (1, -1))]
    r, a = simulate(3, kag, 0, 400)
    print("kagome sim avg:", r[150:].mean(), "(pred 0.1705152873)", a)
