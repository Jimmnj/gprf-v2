"""PREREG_RHOC_RULES_20261030. Reuses marg_fit_20261030 harness and cdr_test_20261030 recipe unchanged."""
import sys,time,numpy as np,pandas as pd
from multiprocessing import Pool
from scipy.optimize import minimize_scalar,brentq
import marg_fit_20261030 as M, cdr_test_20261030 as C
S=M.S; AGRID=np.arange(6.0,8.501,0.05)
BETA,LS0,LSC=0.578,np.log10(1678),np.log10(1263)
def cdr_target(Sstar): return 10**C.eq10([LS0,LSC,BETA],np.log10(Sstar))
def xb0_of(s):
    r,vb2,vobs,ev,Mb=s.curves(M.LUD0,M.LUB0,s.D0,s.i0); vb2c=np.clip(vb2,0,None); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2c]); xs=np.linspace(1e-4,rl,4001)
    phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2c[-1]; s2=np.sqrt(M.a0k*M.G*Mb)/2; return phib0/s2
def score(s,lrc,nuis0):
    m=s.model([lrc]+nuis0,'A')
    if m is None: return None
    r,vb2,vobs,ev,Mb,vm2=m; sel=s.sel; vm2=np.clip(vm2,1,None)
    dlg=np.log10(vobs**2/vm2)[sel]; gb=(np.clip(vb2,0,None)/r)[sel]; low=gb<0.1*M.a0k
    chi=np.sum(((np.sqrt(vm2)-vobs)/ev)[sel]**2)
    return dict(chi2=chi,N=int(sel.sum()),sq=np.sum(dlg**2),sqlow=np.sum(dlg[low]**2),nlow=int(low.sum()),Sdyn=C.toomre(r,vm2))
def one(g):
    s=M.Gal(g); nuis0=[M.LUD0]+([M.LUB0] if s.bul else [])+[0.0,0.0]; xb0=xb0_of(s); Sstar=C.sb0(g); out=[]
    base=dict(gal=g,inc=s.i0,Sstar=Sstar,xb0=xb0,Sb0=S(xb0))
    l0=minimize_scalar(lambda l:s.obj([l]+nuis0,'A'),bounds=(6,10),method='bounded',options={'xatol':0.02}).x
    for arm,l in (('A0',l0),('U0',7.45)):
        sc=score(s,l,nuis0); out.append(dict(base,arm=arm,lrc=l,**(sc or {})))
    # R0
    m=s.model(nuis0,'R')
    if m is not None:
        r,vb2,vobs,ev,Mb,vm2=m; Dp=1; Sh=Mb/(2*np.pi*(1.678*s.t.Rd*1e3)**2) if s.t.Rd>0 else np.nan; fg=1.33*s.t.MHI*1e9/Mb
        lR=11.063-0.365*np.log10(Mb)+0.268*np.log10(Sh)-0.929*fg if np.isfinite(Sh) else 7.45
        out.append(dict(base,arm='R0',lrc=lR,**score(s,lR,nuis0)))
    # BG grid
    for lA in AGRID:
        l=lA+S(xb0)/np.log(10); sc=score(s,l,nuis0)
        if sc: out.append(dict(base,arm='BG%.2f'%lA,lrc=l,**sc))
    # CDR rule
    tgt=cdr_target(Sstar); f=lambda l:(np.log10(score(s,l,nuis0)['Sdyn']/tgt) if score(s,l,nuis0) else 1.0)
    try:
        f6,f10=f(6.0),f(10.0)
        if f6*f10<0: lC=brentq(f,6.0,10.0,xtol=0.01); flag=0
        else: lC=6.0 if abs(f6)<abs(f10) else 10.0; flag=1
    except Exception: lC,flag=7.45,2
    sc=score(s,lC,nuis0); out.append(dict(base,arm='CDR',lrc=lC,clip=flag,**(sc or {})))
    return out
if __name__=='__main__':
    T,P,gals=M.load(); M.init(T,P); t0=time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv)>1 else 6,initializer=M.init,initargs=(T,P)) as pool: R=sum(pool.map(one,gals),[])
    D=pd.DataFrame(R); D.to_csv('rhoc_rules_points_20261030.csv',index=False)
    # pick BG constant by summed chi2
    bg=D[D.arm.str.startswith('BG')].groupby('arm').agg(chi2=('chi2','sum'),n=('gal','nunique'),sq=('sq','sum'),N=('N','sum'))
    full=bg[bg.n==bg.n.max()]; bestA=full.chi2.idxmin(); bestA_rms=full.sq.idxmin()
    D.loc[D.arm==bestA,'arm']='BG'; 
    lines=['PREREG_RHOC_RULES run %s  %.0f s  BG: log A = %s (chi2-optimal), %s (rms-optimal)'%(time.strftime('%Y-%m-%d %H:%M'),time.time()-t0,bestA[2:],bestA_rms[2:])]
    a0=D[D.arm=='A0'].set_index('gal'); res={}
    for arm in ['A0','U0','R0','BG','CDR']:
        a=D[(D.arm==arm)&np.isfinite(D.chi2)].set_index('gal'); g=a.index
        rms=np.sqrt(a.sq.sum()/a.N.sum()); rlow=np.sqrt(a.sqlow.sum()/a.nlow.sum()); chiN=np.median(a.chi2/a.N)
        c=a[a.inc<=85]; b,sb,rr=C.fit(np.log10(c.Sstar),np.log10(c.Sdyn)); C2=abs(b[2]-0.578)<=2*np.hypot(sb[2],0.048)
        d=a.lrc-a0.loc[g].lrc; rho=pd.Series(a.lrc.values).corr(pd.Series(a0.loc[g].lrc.values),method='spearman')
        res[arm]=dict(rms=rms,rlow=rlow,beta=b[2],eb=sb[2],C2=C2,chiN=chiN,rho=rho,drms=np.sqrt(np.mean(d**2)),sd=a.lrc.std(),n=len(g),clip=int(a.get('clip',pd.Series(0,index=g)).fillna(0).gt(0).sum()))
        lines.append('%-4s n=%3d  S1 rms %.4f (low-g %.4f)  S2 beta %.3f+-%.3f C2 %s  S3 med chi2/N %.2f  S4 rho_s %.2f rms dlrc %.2f  sd lrc %.2f  clipped %d'%(arm,len(g),rms,rlow,b[2],sb[2],C2,chiN,rho,np.sqrt(np.mean(d**2)),a.lrc.std(),res[arm]['clip']))
    for arm in ['BG','CDR']:
        r=res[arm]; D2=(r['rms']<=res['A0']['rms']+0.010) and (r['rms']<=res['U0']['rms']-0.005)
        lines.append('%s: D1 carries CDR %s | D2 holds RAR %s | D3 replaces fit %s'%(arm,r['C2'],D2,r['C2'] and D2))
    open('rhoc_rules_summary_20261030.txt','w').write('\n'.join(lines)); print('\n'.join(lines))
