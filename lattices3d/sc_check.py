import mpmath as mp
from gen3d import *
mp.mp.dps = 12
E = [(0, 0, (1, 0, 0)), (0, 0, (0, 1, 0)), (0, 0, (0, 0, 1))]
Y, arcs, deg = transfer_row(1, E, 0, 0, False)
pp = mp.fsum(((1 if j == 0 else 0) - Y[j])**2 for j in range(6)) / 4
print("SC +1 channel:", pp, " x2 =", 2 * pp, " (expect 0.26672720352085)")
