import numpy as np, sys
rng=np.random.default_rng(2); N=int(sys.argv[1]); h=rng.integers(1,7,size=(N,N,N))
def relax(h):
    while True:
        t=(h>6)
        if not t.any(): return h
        t=t.astype(np.int64); h=h-6*t
        for ax in range(3):
            s=np.zeros_like(t); sl=[slice(None)]*3; sd=[slice(None)]*3
            sl[ax]=slice(1,None); sd[ax]=slice(None,-1); s[tuple(sl)]+=t[tuple(sd)]
            s2=np.zeros_like(t); s2[tuple(sd)]+=t[tuple(sl)]
            h=h+s+s2
h=relax(h+6)
for _ in range(20000):
    i,j,k=rng.integers(0,N,3); h[i,j,k]+=1
    if h[i,j,k]>6: h=relax(h)
c=N//4; acc=[]
for it in range(150000):
    i,j,k=rng.integers(0,N,3); h[i,j,k]+=1
    if h[i,j,k]>6: h=relax(h)
    if it%100==0: acc.append((h[c:-c,c:-c,c:-c]==1).mean())
a=np.array(acc); print(N,a.mean(),'+-',a.std()/np.sqrt(len(a)/10))
