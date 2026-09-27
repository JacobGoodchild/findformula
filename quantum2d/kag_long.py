from simgraph import simulate
kag = [(0, 1, (0, 0)), (0, 2, (0, 0)), (1, 2, (0, 0)), (1, 0, (1, 0)), (2, 0, (0, 1)), (1, 2, (1, -1))]
r, a = simulate(3, kag, 0, 800)
print("kagome T=800: avg t=401..800:", r[400:].mean(), " avg t=601..800:", r[600:].mean())
