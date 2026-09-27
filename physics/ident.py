import mpmath as mp, itertools, sys
mp.mp.dps = 55
def basis():
    pi=mp.pi; L2=mp.log(2); z3=mp.zeta(3); G=mp.catalan
    return {"1":mp.mpf(1),"pi":pi,"pi^2":pi**2,"1/pi":1/pi,"1/pi^2":1/pi**2,"log2":L2,"z3":z3,"G":G,
            "z3/pi^2":z3/pi**2,"G/pi":G/pi,"log2/pi":L2/pi,"sqrt2":mp.sqrt(2),"sqrt3":mp.sqrt(3),"e":mp.e,
            "gamma":mp.euler,"log3":mp.log(3),"log(pi)":mp.log(pi),"z3/pi":z3/pi,"pi*log2":pi*L2}
def run(x):
    B=basis(); names=list(B)
    hits=[]
    for r in range(1,4):
        for combo in itertools.combinations(names[1:],r):
            vec=[x,B["1"]]+[B[c] for c in combo]
            rel=mp.pslq(vec,maxcoeff=10**5,maxsteps=10**5)
            if rel and rel[0]!=0:
                hits.append((combo,rel))
    return hits
if __name__=="__main__":
    x=mp.mpf(open(sys.argv[1]).read())
    for h in run(x): print(h)
    print("done")
