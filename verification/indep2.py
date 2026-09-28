import numpy as np, itertools
from scipy import integrate
from mpmath import mp, mpf, ellipk, ellipe, gamma, catalan, quad, erfc, sqrt as ms, inf, exp, asin, pi as mpi, log as mlog, acos
opts=dict(epsabs=1e-12,epsrel=1e-12,limit=400)
def E2(f):
    v,_=integrate.nquad(f,[[0,np.pi],[0,np.pi]],opts=[opts,opts]); return v/np.pi**2
c=np.cos; s3=np.sqrt(3); pi=np.pi
# triangular S = c1+c2+cos(k1-k2); over [0,pi]^2 not symmetric -> use full torus via [0,2pi]
def Et(f):
    v,_=integrate.nquad(f,[[0,2*pi],[0,2*pi]],opts=[opts,opts]); return v/(4*pi**2)
S=lambda x,y:c(x)+c(y)+c(x-y)
print("kagome id",Et(lambda x,y:1/(11-S(x,y))), float(ellipk(mpf(5)/32)/(4*ms(2)*mpi)))
print("tri far edge X",Et(lambda x,y:1/(3+S(x,y))), float(2**(mpf(2)/3)*gamma(mpf(1)/3)**3/(8*mpi**2)))
ms_=mpf(1)/2-23*ms(15)/180
print("tri t=7",Et(lambda x,y:1/(7+S(x,y))), float(ellipk(ms_)/(mpf(135)**0.25*mpi)))
t=5.0; m=0.5-(t*t-3)/(2*(t-1)**1.5*np.sqrt(t+3))
print("tri real-mod t=5",Et(lambda x,y:1/(t+S(x,y))), float(2*ellipk(m)/(mpi*((t-1)**3*(t+3))**0.25)))
t=11.0; e=np.sqrt(2*t+3); m=16*e/((e-1)**3*(e+3))
print("Horiguchi t=11",Et(lambda x,y:1/(t-S(x,y))), float(4*ellipk(m)/(mpi*np.sqrt((e-1)**3*(e+3)))))
# king, diag weight w=1: lam = 2(1-c1)+2(1-c2)+2w(1-cos(k1+k2))+2w(1-cos(k1-k2)) = 4-2c1-2c2+4w(1-c1c2)
w=1.0; lam=lambda x,y:4-2*c(x)-2*c(y)+4*w*(1-c(x)*c(y))
R=lambda a,b: 0.0 if (a,b)==(0,0) else E2(lambda x,y:2*(1-c(a*x)*c(b*y))/lam(x,y))
L=np.log(2+s3)/(pi*s3); p=1/pi
Rv={}
for v in [(1,0),(1,1),(2,0),(2,1),(2,2)]: Rv[v]=R(*v)
print("king R10",Rv[(1,0)],L," R11",Rv[(1,1)],0.5-L," R20",Rv[(2,0)],4/pi-4*L)
print("king R21",Rv[(2,1)],6*L-0.5-2/pi," R22",Rv[(2,2)],5-4/pi-14*L)
def Rk(v):
    a,b=sorted(map(abs,v),reverse=True); return 0.0 if (a,b)==(0,0) else Rv[(a,b)]
nb=[(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]
Y=np.array([[(Rk(a)+Rk(b)-Rk((a[0]-b[0],a[1]-b[1])))/2 for b in nb] for a in nb])
print("king P1 indep",np.linalg.det(np.eye(7)-Y[1:,1:]),"paper",-(3+4*p-14*L)*(1+4*p-2*L)*(6*L+12*p-7)*(13-32*p+48*p**2-52*L-12*L**2)**2/2048)
G=float(catalan)
print("king entropy",E2(lambda x,y:np.log(lam(x,y))), 20*G/(3*pi))
print("snub E ln det",E2(lambda x,y:np.log(2*(170-72*(c(x)+c(y))-28*c(x)*c(y)+c(2*x)+c(2*y)))),40*G/(3*pi)+4/3*np.log(2+s3))
# honeycomb+NNN w=1: R_NN = ((3w+1)/(3w^2)) g_tri(3+(6w+1)/(2w^2))
tw=3+7/2; print("hc+nnn R_NN",4/3*Et(lambda x,y:1/(tw-S(x,y))), float(16*ellipk(mpf(64)/189)/(9*ms(21)*mpi)))
# square with 4 conductances, random: direct R(e1)=E[2(1-c1)/lam]
w1,w2,wp,wm=1,0.43,0.4,0.45
lam4=lambda x,y:2*w1*(1-c(x))+2*w2*(1-c(y))+2*wp*(1-c(x+y))+2*wm*(1-c(x-y))
d,_=integrate.nquad(lambda x,y:2*(1-c(x))/lam4(x,y),[[-pi,pi],[-pi,pi]],opts=[opts,opts]); d/=4*pi**2
sg=wp+wm; D=4*wp*wm-w1**2
print("Formula30",d,(2/(pi*np.sqrt(D)))*np.arcsinh(np.sqrt(D/((w1+sg)*(w1+w2)))) if D>0 else (2/(pi*np.sqrt(-D)))*np.arcsin(np.sqrt(-D/((w1+sg)*(w1+w2)))))
mp.dps=25
print("A4",quad(lambda s:erfc(ms(s))**4,[0,1,inf]), mpf(1)/2-(18-8*ms(3))/(3*mpi))
b=mpf(1)/3; print("lemma", quad(lambda x:exp(-b*x*x)*erfc(x)**2,[0,inf]), ms(mpi/b)/2*(1+2/mpi*(asin(1/(b+1))-2*asin(1/ms(b+1)))))
mp.dps=20; A=mpf(1)/2-(3-ms(3))/mpi; B=(1-3*A)/3; print("3D pbar",3*A*A+B*B)
