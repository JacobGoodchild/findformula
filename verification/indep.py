import numpy as np, itertools
from scipy import integrate
from mpmath import mp, mpf, gamma, pi as mpi, sqrt as msqrt, ellipk, log as mlog, catalan
opts=dict(epsabs=1e-11,epsrel=1e-11,limit=400)
def E2(f):  # average over [0,pi]^2 of even integrand == average over full torus
    v,_=integrate.nquad(f,[[0,np.pi],[0,np.pi]],opts=[opts,opts]); return v/np.pi**2
def E2full(f):
    v,_=integrate.nquad(f,[[-np.pi,np.pi],[-np.pi,np.pi]],opts=[opts,opts]); return v/(4*np.pi**2)
c=np.cos
# ---- SC: R(x)=E[(1-cos k.x)/(3-sum c)], k3 done exactly
def sc_R(n):  # n = (a,b) in-plane displacement with zero k3 component
    return E2(lambda x,y: (1-c(n[0]*x)*c(n[1]*y))/np.sqrt((3-c(x)-c(y))**2-1))
R100,R110,R200=sc_R((1,0)),sc_R((1,1)),sc_R((2,0))
W=float(msqrt(6)/(32*mpi**3)*gamma(mpf(1)/24)*gamma(mpf(5)/24)*gamma(mpf(7)/24)*gamma(mpf(11)/24))
al=(7*W+54/(np.pi**2*W))/36
print("SC R100",R100,"(1/3)"); print("SC R110",R110,"alpha",al); print("SC R200",R200,"2-4a",2-4*al)
def P1(nbrs,Rf):
    z=len(nbrs); Y=np.array([[ (Rf(a)+Rf(b)-Rf(tuple(np.subtract(a,b))))/2 for b in nbrs] for a in nbrs])
    Yp=Y[1:,1:]; return np.linalg.det(np.eye(z-1)-Yp), Y
def sym(R):  # resistance lookup by sorted |coords|
    return lambda v: 0.0 if not any(v) else R[tuple(sorted(map(abs,v),reverse=True))]
nb=[s*np.eye(3,dtype=int)[i] for i in range(3) for s in (1,-1)]; nb=[tuple(v) for v in nb]
p,_=P1(nb,sym({(1,0,0):R100,(1,1,0):R110,(2,0,0):R200}))
print("SC P1 from integrals",p,"paper",4/3*al**3*(2-3*al)**2)
# ---- BCC: L=8(1-c1c2c3), R(x)=E[(1-cos k.x)/(4(1-c1c2c3))]
def bcc_R(n): return E2(lambda x,y:(1-c(n[0]*x)*c(n[1]*y))/(4*np.sqrt(1-(c(x)*c(y))**2)))
Wb=float(gamma(0.25)**4/(4*mpi**3)); q=1/(np.pi**2*Wb)
B200,B220=bcc_R((2,0)),bcc_R((2,2)); B111=0.25; B222=(8-3*Wb-36*q)/4  # R(2,2,2) from paper (Foster-consistent)
print("BCC R200",B200,"paper",(Wb-4*q)/4,"| R220",B220,"paper",4*q)
nbb=[tuple(v) for v in itertools.product((1,-1),repeat=3)]
p,_=P1(nbb,sym({(1,1,1):B111,(2,0,0):B200,(2,2,0):B220,(2,2,2):B222}))
print("BCC P1 (R222 from paper)",p,"paper",3*q/128*(Wb+4*q)**3*(4-Wb-12*q)**3)
# ---- FCC: L=12-4s, R(x)=E[(1-cos k.x)/(6-2s)], s=c1c2+c3(c1+c2)
def fcc_AB(x,y): return 6-2*c(x)*c(y), 2*(c(x)+c(y))
def fcc_R_plane(n):
    return E2(lambda x,y:(1-c(n[0]*x)*c(n[1]*y))/np.sqrt(fcc_AB(x,y)[0]**2-fcc_AB(x,y)[1]**2))
def fcc_R211():
    def f(x,y):
        A,B=fcc_AB(x,y); s=np.sqrt(A*A-B*B); Ec=(A/s-1)/B if abs(B)>1e-12 else 0.0
        return 1/s - c(2*x)*c(y)*Ec
    return E2(f)
Wf=float(3*gamma(mpf(1)/3)**6/(2**(mpf(14)/3)*mpi**4)); qf=1/(np.pi**2*Wf)
F200,F220,F211=fcc_R_plane((2,0)),fcc_R_plane((2,2)),fcc_R211()
F110=fcc_R_plane((1,1))
print("FCC R110",F110,"(1/6) R200",F200,"paper",(4*Wf-3*qf)/6)
print("FCC R220",F220,"paper",(8-12*Wf-9*qf)/3,"| R211",F211,"paper",(2*Wf+3*qf-1)/3)
nbf=sorted(set(tuple(p) for v in [(1,1,0),(1,-1,0),(-1,1,0),(-1,-1,0)] for p in itertools.permutations(v)))
def fR(v):
    if not any(v): return 0.0
    k=tuple(sorted(map(abs,v),reverse=True))
    return {(1,1,0):F110,(2,0,0):F200,(2,2,0):F220,(2,1,1):F211}[k]
p,_=P1(nbf,fR)
print("FCC P1 from integrals",p,"paper",(5-4*Wf-6*qf)**2*(7-8*Wf-3*qf)**3*(8*Wf+3*qf+1)**3*(16*Wf+15*qf-5)**3/6**10)
