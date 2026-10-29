# n16_closure_20261028.py -- N16 run under PREREG_N16_20261028.md (STAMPED JM 2026-10-28 00:24 "accepted - go"; FROZEN). [UNRATIFIED result]
# Closure: rho_c r_c = mu0_M (Milgrom 2009 eq 11, RAR mu, Sigma_b = 0.75 Sigma_h) -> rho_c = 4 pi G mu0_M^2 / (9 sigma^2), sigma^2 = sqrt(G a0 M_b)/2.
# Parts (argv): S1 (C0 + OLS + diagnostics), S3 (dwarfs C4/C5 + diagnostics), S2 (SPARC C1-C3, cases univ/rel3/clos), S2k (diagnostic: best constant k).
import sys, io, contextlib, warnings, numpy as np, pandas as pd; warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from gate2c_robustness_20260925 import load_table
PART=sys.argv[1]
Gk=4.30091e-6; a0SI=1.2e-10; kpc_m=3.0857e19; a0k=a0SI*kpc_m/1e6
SigM=a0SI/(2*np.pi*6.674e-11)*(3.0857e16**2)/1.989e30
nuRAR=lambda y: 1/(1-np.exp(-np.sqrt(y)))
nuSIM=lambda y: 0.5*(1+np.sqrt(1+4/y))
def make_sigc(nu):
    yg=np.logspace(-8,4,200001); Xg=yg*nu(yg); mug=1/nu(yg)
    Lg=np.gradient(np.log(mug),np.log(Xg))
    cumL=np.concatenate([[Xg[0]*Lg[0]],Xg[0]*Lg[0]+np.cumsum(0.5*(Lg[1:]+Lg[:-1])*np.diff(Xg))])
    def sigc(Sb):
        y=(2/3)*np.asarray(Sb,float)/SigM; X0=y*nu(y); m0=1/nu(y)
        return SigM*((6/np.pi)*X0*(1-m0)+np.interp(X0,Xg,cumL))
    return sigc,cumL[-1]
sigRAR,lamR=make_sigc(nuRAR); sigSIM,lamS=make_sigc(nuSIM)
def lrc_clos(Mb,Sh,fb=0.75,sig=sigRAR,lk=0.0):
    s2=0.5*np.sqrt(Gk*a0k*np.asarray(Mb,float)); mu=sig(fb*np.asarray(Sh,float))*1e6*10**lk   # Msun/kpc^2
    return np.log10(4*np.pi*Gk*mu**2/(9*s2))
print('Sigma_M %.1f (log %.3f); lambda RAR %.3f (gamma %.3f); lambda simple %.3f (gamma %.3f)'%(SigM,np.log10(SigM),lamR,2*lamR/np.pi,lamS,2*lamS/np.pi))
T=load_table()

if PART=='S1':
    F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal').join(T[['Rd','L','MHI']]); F=F[F.Rd>0]
    lM=np.log10(F.Mb.values); Sh=F.Mb.values/(2*np.pi*(1.678*F.Rd.values*1e3)**2); lS=np.log10(Sh); fg=1.33*F.MHI.values*1e9/F.Mb.values; y=F.l_tanh2.values
    X=np.column_stack([np.ones(len(y)),lM,lS,fg]); beta=np.linalg.lstsq(X,y,rcond=None)[0]; r3=X@beta
    s2chk=0.5*np.sqrt(Gk*a0k*F.Mb.values); print('S1 N=%d; eq.(2) sigma^2 vs TEST1 s2b: max |log ratio| %.2e'%(len(y),np.max(np.abs(np.log10(s2chk/F.s2b.values)))))
    print('relation (3) refit: %.3f %+.3f lMb %+.3f lSh %+.3f fg; rms vs fitted %.3f (sd %.3f)'%(*beta,np.sqrt(np.mean((y-r3)**2)),np.std(y-r3)))
    lc=lrc_clos(F.Mb.values,Sh); d=lc-y
    print('\nC0  closure vs fitted: rms %.3f  mean %+.3f  sd %.3f  median %+.3f   [PASS <=0.64, FAIL >0.80]'%(np.sqrt(np.mean(d**2)),d.mean(),d.std(),np.median(d)))
    print('    closure log rho_c: median %.2f, range %.2f..%.2f ; fitted median %.2f range %.2f..%.2f'%(np.median(lc),lc.min(),lc.max(),np.median(y),y.min(),y.max()))
    print('    Spearman(closure-fitted, log Sigma_h) %+.2f  (log M_b) %+.2f  (f_gas) %+.2f'%(spearmanr(d,lS)[0],spearmanr(d,lM)[0],spearmanr(d,fg)[0]))
    bc=np.linalg.lstsq(X,lc,rcond=None)[0]
    print('    OLS closure: %.3f %+.3f lMb %+.3f lSh %+.3f fg   (pre-stated expectation -0.50 / +0.80 / 0; rel3 %+.3f / %+.3f / %+.3f)'%(*bc,*beta[1:]))
    q=pd.qcut(lS,4,labels=['lowS','q2','q3','highS']); print('    median (closure-fitted) by Sigma_h quartile:',pd.Series(d).groupby(q).median().round(2).to_dict())
    print('\nDIAGNOSTICS (no verdict):')
    for lab,kw in (('Sigma_b=0.50 Sigma_h',dict(fb=0.5)),('Sigma_b=1.00 Sigma_h',dict(fb=1.0)),('simple mu',dict(sig=sigSIM))):
        dd=lrc_clos(F.Mb.values,Sh,**kw)-y; print('  %-22s C0 rms %.3f mean %+.3f sd %.3f'%(lab,np.sqrt(np.mean(dd**2)),dd.mean(),dd.std()))
    lk=-d.mean()/2; dd=lrc_clos(F.Mb.values,Sh,lk=lk)-y
    print('  best constant k (rho_c r_c = k mu0_M), fitted on S1: log k = %+.3f (k = %.2f); C0 rms %.3f'%(lk,10**lk,np.sqrt(np.mean(dd**2))))
    pd.DataFrame(dict(gal=F.index,lMb=lM,lSh=lS,fgas=fg,lrc_fit=y,lrc_rel3=r3,lrc_clos=lc)).to_csv('n16_S1_20261028.csv',index=False)

if PART=='S3':
    src=open('dwarf_prereg_20261020.py').read(); head=src.split('rows=[]')[0]
    with contextlib.redirect_stdout(io.StringIO()): exec(head)
    def run(fb=0.75,sig=sigRAR,lk=0.0,verbose=False):
        rows=[]
        for name,LV,rh,so,ep,em,vc,Dh,cls in dw:
            Mb=1.6*LV; Shd=Mb/(2*np.pi*(rh*1e3)**2); lrc=float(lrc_clos(Mb,Shd,fb,sig,lk))
            a=rh/1.305; s2=0.5*np.sqrt(G*a0k*Mb); amb=host_rho(vc,Dh)+rho_bg
            rmax=max(8*rh,3*np.sqrt(s2/(2*np.pi*G*amb))); s=solve(Mb,a,s2,amb,lrc,rmax); sg=sigma_los(s,a,rh,rmax)
            err=np.sqrt(0.10**2+(0.5*(ep+em)/so/np.log(10))**2); dl=np.log10(sg/so); ok=abs(dl)<=err
            rows.append(dict(name=name,cls=cls,Mb=Mb,rh=rh,lSig=np.log10(Shd),lrc=lrc,s_obs=so,s_clos=sg,dlog=dl,bar=err,PASS=ok))
            if verbose: print('%-10s %-3s Mb %.1e r_h %.2f lSig %.2f | rho_c %.2f | obs %5.1f | clos %5.1f (%+.2f dex, bar %.2f) %s'%(name,cls,Mb,rh,np.log10(Shd),lrc,so,sg,dl,err,'PASS' if ok else 'FAIL'))
        return pd.DataFrame(rows)
    R=run(verbose=True); R.to_csv('n16_S3_dwarfs_20261028.csv',index=False)
    emb=R[R.cls!='ISO'].dlog.mean(); iso=R[R.cls=='ISO'].dlog.mean(); kill=(emb>0.15) and not (iso>0.15)
    print('\nC4  %d/16 within bar   [PASS >=11]; rms log(pred/obs) %.3f mean %+.3f (rel3 0.148 / +0.028)'%(R.PASS.sum(),np.sqrt((R.dlog**2).mean()),R.dlog.mean()))
    print('C5  embedded mean log(pred/obs) %+.3f (N=%d), isolated %+.3f (N=%d); kill rule triggered: %s'%(emb,(R.cls!='ISO').sum(),iso,(R.cls=='ISO').sum(),'YES' if kill else 'NO'))
    print('    closure rho_c range %.2f..%.2f (rel3 dwarfs 8.37..9.10)'%(R.lrc.min(),R.lrc.max()))
    print('\nDIAGNOSTICS (no verdict):')
    for lab,kw in (('Sigma_b=0.50 Sigma_h',dict(fb=0.5)),('Sigma_b=1.00 Sigma_h',dict(fb=1.0)),('simple mu',dict(sig=sigSIM))):
        r=run(**kw); print('  %-22s %d/16, rms %.3f mean %+.3f'%(lab,r.PASS.sum(),np.sqrt((r.dlog**2).mean()),r.dlog.mean()))
    if len(sys.argv)>2:
        lk=float(sys.argv[2]); r=run(lk=lk); print('  best-k (log k %+.3f from S1)    %d/16, rms %.3f mean %+.3f'%(lk,r.PASS.sum(),np.sqrt((r.dlog**2).mean()),r.dlog.mean()))

if PART in ('S2','S2k'):
    src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
    with contextlib.redirect_stdout(io.StringIO()): exec(src)
    beta=np.array([11.063,-0.365,0.268,-0.929])
    Lsel=(P.Q<=2)&(P.inc>=30)&(P.ev/P.vobs<0.1)&(P.vobs>0)&(P.vbar2>0); gals=sorted(P[Lsel].gal.unique())
    lk=float(sys.argv[2]) if PART=='S2k' else 0.0
    def gpar(g):
        t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9; Sh=Mb/(2*np.pi*(1.678*t.Rd*1e3)**2) if t.Rd>0 else np.nan; fg=1.33*t.MHI*1e9/Mb
        return Mb,Sh,fg
    cases=(('clos',None),) if PART=='S2k' else (('univ',None),('rel3',None),('clos',None))
    out=[]
    for g in gals:
        Mb,Sh,fg=gpar(g)
        for lab,_ in cases:
            if lab=='univ': l=7.45
            elif lab=='rel3': l=beta[0]+beta[1]*np.log10(Mb)+beta[2]*np.log10(Sh)+beta[3]*fg if np.isfinite(Sh) else 7.45
            else: l=float(lrc_clos(Mb,Sh,lk=lk)) if np.isfinite(Sh) else 7.45
            o=solve_gal(g,l)
            if o is None: continue
            o=o[(o.ev/o.vobs<0.1)&(o.vobs>0)&(o.vbar2>0)].copy(); o['case']=lab; o['lrc']=l; out.append(o)
    A=pd.concat(out); A['lgb']=np.log10(A.vbar2/A.r*1e6/kpc); A['lgo']=np.log10(A.vobs**2/A.r*1e6/kpc); A['lgm']=np.log10(A.vmod2/A.r*1e6/kpc)
    A['lrar']=np.log10(10**A.lgb/(1-np.exp(-np.sqrt(10**A.lgb/a0)))); A['res']=A.lgo-A.lgm; A['resR']=A.lgo-A.lrar
    A.to_csv('n16_closure_points_20261028.csv' if PART=='S2' else 'n16_closure_points_k_20261028.csv',index=False)
    la0=np.log10(a0); edges=np.arange(-2.25,1.01,0.25)+la0; labs=[c for c,_ in cases]
    print('binned median log(g_obs/g_model), %d galaxies; cases %s'%(A.gal.nunique(),labs))
    for lo,hi in zip(edges[:-1],edges[1:]):
        s=[A[(A.case==c)&(A.lgb>=lo)&(A.lgb<hi)] for c in labs]
        if len(s[0])<10: continue
        print('  %+5.2f..%+5.2f %4d   '%(lo-la0,hi-la0,len(s[0]))+'   '.join('%s %+.3f'%(c,x.res.median()) for c,x in zip(labs,s)))
    for lab in labs:
        d=A[A.case==lab]
        print('%-5s C1 per-point rms %.3f (N %d) | C2 Spearman res vs log Sigma_h %+.2f | vs log r %+.2f, f_gas %+.2f, log M_b %+.2f'%(lab,np.sqrt((d.res**2).mean()),len(d),
              spearmanr(np.log10(d.Sh),d.res)[0],spearmanr(np.log10(d.r),d.res)[0],spearmanr(d.fgas,d.res)[0],spearmanr(np.log10(d.Mb),d.res)[0]))
        G2=[]
        for g,dd in d.groupby('gal'):
            last=dd.sort_values('r').iloc[-1]; G2.append(dict(gal=g,Sh=last.Sh,ratio=np.sqrt(last.vmod2)/last.vobs))
        G2=pd.DataFrame(G2); q=pd.qcut(np.log10(G2.Sh),4,labels=['lowSigma','q2','q3','highSigma']); m=G2.groupby(q).ratio.median()
        print('      C3 V_mod/V_obs at R_last by Sigma_h quartile %s ; spread low-high %+.3f ; overall median %.3f'%(m.round(3).to_dict(),m['lowSigma']-m['highSigma'],G2.ratio.median()))
