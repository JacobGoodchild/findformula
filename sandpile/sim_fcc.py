import numpy as np, sys
rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 3); N=int(sys.argv[1]); z=12
OFF=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,-1,0),(-1,1,0),(0,1,-1),(0,-1,1),(1,0,-1),(-1,0,1)]
def shift(t,o):
    s=np.zeros_like(t); src=[];dst=[]
    for d in o:
        src.append(slice(0,N-d) if d>0 else slice(-d,N)); dst.append(slice(d,N) if d>0 else slice(0,N+d))
    s[tuple(dst)]=t[tuple(src)]; return s
h=rng.integers(1,z+1,size=(N,N,N))
def relax(h):
    while True:
        t=(h>z)
        if not t.any(): return h
        t=t.astype(np.int64); h=h-z*t
        for o in OFF: h=h+shift(t,o)
h=relax(h+z)
for _ in range(20000):
    i,j,k=rng.integers(0,N,3); h[i,j,k]+=1
    if h[i,j,k]>z: h=relax(h)
c=N//4; acc=[]
for it in range(int(sys.argv[3]) if len(sys.argv)>3 else 150000):
    i,j,k=rng.integers(0,N,3); h[i,j,k]+=1
    if h[i,j,k]>z: h=relax(h)
    if it%100==0: acc.append((h[c:-c,c:-c,c:-c]==1).mean())
a=np.array(acc); print(N,a.mean(),'+-',a.std()/np.sqrt(len(a)/10))
