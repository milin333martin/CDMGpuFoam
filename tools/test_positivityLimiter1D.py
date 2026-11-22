import numpy as np
g=1.4
def prim(U):
    r=U[0]; u=U[1]/r; p=(g-1)*(U[2]-0.5*r*u*u); return r,u,p
def cons(r,u,p): return np.array([r,r*u,p/(g-1)+0.5*r*u*u])
def pad(a): return np.concatenate([[a[0]],a,[a[-1]]])
def vl(a):
    dm=a[1:-1]-a[:-2]; dp=a[2:]-a[1:-1]; s=np.where(dm*dp>0,2*dm*dp/(dm+dp+1e-300),0)
    return np.concatenate([[a[0]],a[1:-1]+0.5*s]), np.concatenate([a[1:-1]-0.5*s,[a[-1]]])
def euler(r,u,p):
    E=p/(g-1)+0.5*r*u*u; return np.array([r*u,r*u*u+p,(E+p)*u])
def hllc(rL,uL,pL,rR,uR,pR):
    cL=np.sqrt(g*pL/rL); cR=np.sqrt(g*pR/rR)
    sL_=np.sqrt(rL); sR_=np.sqrt(rR); k=1/(sL_+sR_); uh=(uL*sL_+uR*sR_)*k
    ch=np.sqrt((cL**2*sL_+cR**2*sR_)*k+0.5*sL_*sR_*k*k*(uR-uL)**2)
    SL=np.minimum(uL-cL,uh-ch); SR=np.maximum(uR+cR,uh+ch)
    aL=rL*(SL-uL); aR=rR*(SR-uR); Ss=(pR-pL+aL*uL-aR*uR)/(aL-aR)
    FL=euler(rL,uL,pL); FR=euler(rR,uR,pR); UL=cons(rL,uL,pL); UR=cons(rR,uR,pR)
    def star(r,u,p,S,a):
        k=a/(S-Ss); E=p/(g-1)/r+0.5*u*u
        return k*np.array([np.ones_like(r),Ss,E+(Ss-u)*(Ss+p/a)])
    FsL=FL+SL*(star(rL,uL,pL,SL,aL)-UL); FsR=FR+SR*(star(rR,uR,pR,SR,aR)-UR)
    return np.where(SL>=0,FL,np.where(SR<=0,FR,np.where(Ss>=0,FsL,FsR)))
def rhs(U,dx,dt,limit):
    r,u,p=prim(U); c=np.sqrt(g*np.maximum(p,1e-300)/r)
    rL,rR=vl(pad(r)); uL,uR=vl(pad(u)); pL,pR=vl(pad(p))
    F=hllc(rL,uL,pL,rR,uR,pR)            # N+1 faces
    nlim=0
    if limit:
        Up=np.concatenate([U[:,:1],U,U[:,-1:]],axis=1); rp,up,pp=prim(Up); cp=np.sqrt(g*pp/rp)
        Uo,Un=Up[:,:-1],Up[:,1:]
        a=np.maximum(abs(up[:-1])+cp[:-1],abs(up[1:])+cp[1:])
        FLF=0.5*(euler(rp[:-1],up[:-1],pp[:-1])+euler(rp[1:],up[1:],pp[1:]))+0.5*a*(Uo-Un)
        W=a[:-1]+a[1:]                    # per cell
        th=np.ones(len(a))
        for which in (0,1):
            Fc=FLF+th*(F-FLF); t=np.ones(len(a))
            for side,sgn in ((0,1.0),(1,-1.0)):   # owner (cell left of face), neighbour
                cells=np.arange(len(a))-1+side   # cell index for owner (i-1) / neighbour (i)
                ok=(cells>=0)&(cells<U.shape[1])
                ci=np.clip(cells,0,U.shape[1]-1)
                k=dt*W[ci]/(dx*a)
                Uho=U[:,ci]-k*sgn*Fc; Ulf=U[:,ci]-k*sgn*FLF
                if which==0: gho=Uho[0]; glf=Ulf[0]; eps=1e-3*1.0
                else:
                    gho=prim(Uho)[2]; glf=prim(Ulf)[2]; eps=1e-3*0.4
                tt=np.where(gho>=eps,1.0,np.where(glf<=eps,0.0,np.clip((glf-eps)/(glf-gho+1e-300),0,1)))
                t=np.minimum(t,np.where(ok,tt,1.0))
            th=th*t
        nlim=np.sum(th<1); F=FLF+th*(F-FLF)
    return -(F[:,1:]-F[:,:-1])/dx, nlim
def run(limit, u0=4.0, N=400):
    dx=1.0/N; x=(np.arange(N)+0.5)*dx
    U=cons(np.ones(N),np.where(x<0.5,-u0,u0),0.4*np.ones(N))
    t=0; T=0.1; nl=0; rmin=1
    while t<T:
        r,u,p=prim(U); c=np.sqrt(g*p/r); dt=min(0.4*dx/np.max(abs(u)+c),T-t)
        L1,n1=rhs(U,dx,dt,limit); U1=U+dt*L1
        if np.any(prim(U1)[0]<=0) or np.any(prim(U1)[2]<=0) or not np.all(np.isfinite(U1)): return "FAILED at t=%.4f"%t, nl
        L2,n2=rhs(U1,dx,dt,limit); U=0.5*U+0.5*(U1+dt*L2); nl+=n1+n2
        r,u,p=prim(U)
        if np.any(r<=0) or np.any(p<=0) or not np.all(np.isfinite(U)): return "FAILED at t=%.4f"%t, nl
        rmin=min(rmin,r.min()); t+=dt
    return "OK  min rho = %.3e, min p = %.3e"%(rmin, prim(U)[2].min()), nl
for u0 in (2.0, 4.0, 6.0):
    print(f"u0 = {u0}: no limiter: {run(False,u0)[0]:45s} | limiter: {run(True,u0)[0]}  (limited faces {run(True,u0)[1]})")
