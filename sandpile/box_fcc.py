"""Exact (no sampling) P1 at the centre of a finite FCC box with dissipative boundary (Dirichlet toppling
matrix): P1 = det(I - M) with M from the finite-box Green function. Converges to the infinite-lattice value."""
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl, sys, itertools
OFF=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),(1,-1,0),(-1,1,0),(0,1,-1),(0,-1,1),(1,0,-1),(-1,0,1)]
def P1(N, OFF=OFF):
    z=len(OFF); idx=lambda x,y,w: x+N*(y+N*w); n=N**3; R=[];C=[]
    for x,y,w in itertools.product(range(N),repeat=3):
        for o in OFF:
            X,Y,Z=x+o[0],y+o[1],w+o[2]
            if 0<=X<N and 0<=Y<N and 0<=Z<N: R.append(idx(x,y,w)); C.append(idx(X,Y,Z))
    A=sps.csr_matrix((np.ones(len(R)),(R,C)),shape=(n,n))
    Delta=(z*sps.identity(n)-A).tocsc(); lu=spl.splu(Delta)
    c=N//2; o0=idx(c,c,c); nb=[idx(c+o[0],c+o[1],c+o[2]) for o in OFF]; cut=nb[1:]
    B=np.zeros((n,len(cut)))
    for j,v in enumerate(cut): B[o0,j]=1; B[v,j]=-1
    GB=np.column_stack([lu.solve(B[:,j]) for j in range(len(cut))])
    M=B.T@GB
    return np.linalg.det(np.eye(len(cut))-M)
if __name__=="__main__":
    for N in map(int,sys.argv[1:]): print(N, P1(N))
