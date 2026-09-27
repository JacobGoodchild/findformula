import sys, mpmath as mp, pickle
sys.argv = ["x", "32"]
exec(open("tri_moving.py").read().split('if __name__ == "__main__":')[0])
a = arcs[0]
vals = {}
for b in [(-1, 0), (1, 1), (-1, 1)]:
    vals[b] = mp.re(avg(entry(1, b, a))); print(b, vals[b], flush=True)
pickle.dump(vals, open("tri_moving_plus32.pkl", "wb"))
