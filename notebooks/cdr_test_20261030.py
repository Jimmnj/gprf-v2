"""PREREG_CDR_20261030: central density relation test. Reuses marg_fit_20261030 harness unchanged."""
import os,sys,glob,numpy as np,pandas as pd,time
from multiprocessing import Pool
from scipy.optimize import minimize_scalar
import marg_fit_20261030 as M
from scipy import odr
G=M.G; Q0=0.2
def toomre(r,v2):
    r=np.asarray(r,float); v2=np.clip(np.asarray(v2,float),0,None)
    edges=np.concatenate([[0],(r[1:]+r[:-1])/2,[r[-1]+(r[-1]-r[-2])/2 if len(r)>1 else 1.5*r[-1]]])
    dR=np.diff(edges); s=np.sum(v2/r**2*dR)+v2[-1]/r[-1]
    return (1+Q0)/(2*np.pi*G)*s*1e-6   # Msun/pc^2
def sb0(g):
    d=np.atleast_2d(np.loadtxt(os.path.join(M.SD,g+'_rotmod.dat'),comments='#')); d=d[np.argsort(d[:,0])]
    return 0.5*d[0,6]+0.7*d[0,7]
def one(g):
    s=M.Gal(g); nuis0=[M.LUD0]+([M.LUB0] if s.bul else [])+[0.0,0.0]; out=[]
    r,vb2,vobs,ev,Mb=s.curves(M.LUD0,M.LUB0,s.D0,s.i0); Sobs=toomre(r,vobs**2); Sstar=sb0(g)
    inc=s.i0; out.append(dict(gal=g,arm='OBS',Sstar=Sstar,Sdyn=Sobs,inc=inc,Mb=Mb,lrc=np.nan))
    l0=minimize_scalar(lambda l:s.obj([l]+nuis0,'A'),bounds=(6,10),method='bounded',options={'xatol':0.02}).x
    for arm,x in (('A0',[l0]+nuis0),('U0',nuis0),('R0',nuis0),('B0',nuis0)):
        m=s.model(x,arm[0])
        if m is None: out.append(dict(gal=g,arm=arm,Sstar=Sstar,Sdyn=np.nan,inc=inc,Mb=Mb,lrc=np.nan)); continue
        out.append(dict(gal=g,arm=arm,Sstar=Sstar,Sdyn=toomre(m[0],m[5]),inc=inc,Mb=Mb,lrc=(l0 if arm=='A0' else np.nan)))
    # secondary: MAP parameters
    for arm,fn in (('A','marg_fit_params_20261030.csv'),('B','marg_fit_params_20261030.csv'),('U','marg_supp_params_20261030.csv'),('R','marg_supp_params_20261030.csv')):
        P=pd.read_csv(fn); row=P[(P.gal==g)&(P.arm==arm)]
        if not len(row): continue
        row=row.iloc[0]; x=([row.lrc] if arm=='A' else [])+[row.lud]+([row.lub] if s.bul else [])+[row.zD,row.zi]
        m=s.model(x,arm)
        if m is None: continue
        r2,vb2,vobs2,ev2,Mb2,vm2=m
        out.append(dict(gal=g,arm='OBS_'+arm,Sstar=Sstar,Sdyn=toomre(r2,vobs2**2),inc=row.inc,Mb=Mb2,lrc=np.nan))
        out.append(dict(gal=g,arm=arm+'_map',Sstar=Sstar,Sdyn=toomre(r2,vm2),inc=row.inc,Mb=Mb2,lrc=row.lrc if arm=='A' else np.nan))
    return out
def eq10(p,ls): lS0,lSc,beta=p; x=10**ls/10**lSc; return lS0+(1-beta)*np.log10(1+x)+beta*np.log10(x)
def fit(ls,ld):
    d=odr.RealData(ls,ld,sx=0.14,sy=0.14); o=odr.ODR(d,odr.Model(eq10),beta0=[np.log10(1468),np.log10(1571),0.62]).run()
    res=ld-eq10(o.beta,ls); return o.beta,o.sd_beta,np.sqrt(np.mean(res**2))
if __name__=='__main__':
    T,P,gals=M.load(); M.init(T,P); t0=time.time()
    with Pool(int(sys.argv[1]) if len(sys.argv)>1 else 4,initializer=M.init,initargs=(T,P)) as pool: R=sum(pool.map(one,gals),[])
    D=pd.DataFrame(R); D=D[D.inc<=85]; D.to_csv('cdr_points_20261030.csv',index=False)
    obs=D[D.arm=='OBS'].set_index('gal'); rows=[]; lines=['PREREG_CDR_20261030 run %s  N(OBS, i<=85)=%d  %.0f s'%(time.strftime('%Y-%m-%d %H:%M'),len(obs),time.time()-t0)]
    bo,so,rmso=fit(np.log10(obs.Sstar),np.log10(obs.Sdyn))
    for arm in ['OBS','A0','U0','R0','B0','A_map','U_map','R_map','B_map']:
        a=D[(D.arm==arm)&np.isfinite(D.Sdyn)&(D.Sdyn>0)].set_index('gal')
        if not len(a): continue
        ref=obs if not arm.endswith('_map') else D[D.arm=='OBS_'+arm[0]].set_index('gal')
        g=a.index.intersection(ref.index); dl=np.log10(a.loc[g].Sdyn)-np.log10(ref.loc[g].Sdyn)
        b,sb,rms=fit(np.log10(a.loc[g].Sstar),np.log10(a.loc[g].Sdyn))
        C1=(np.sqrt(np.mean(dl**2))<=0.10)&(abs(np.median(dl))<=0.05); C2=abs(b[2]-bo[2])<=2*np.hypot(sb[2],so[2]); C3=rms<=1.5*rmso
        rows.append(dict(arm=arm,N=len(g),lS0=b[0],lScrit=b[1],beta=b[2],ebeta=sb[2],rms_resid=rms,rms_delta=np.sqrt(np.mean(dl**2)),med_delta=np.median(dl),C1=C1,C2=C2,C3=C3,D1=(C1&C2&C3)))
        lines.append('%-6s N=%3d  beta=%.3f+-%.3f  Sig0=%6.0f Sigcrit=%6.0f  resid rms %.3f | Delta vs OBS: rms %.3f med %+.3f | C1 %s C2 %s C3 %s -> D1 %s'%(arm,len(g),b[2],sb[2],10**b[0],10**b[1],rms,np.sqrt(np.mean(dl**2)),np.median(dl),C1,C2,C3,C1&C2&C3))
    F=pd.DataFrame(rows); F.to_csv('cdr_fits_20261030.csv',index=False)
    u=F[F.arm=='U0'].iloc[0]; a=F[F.arm=='A0'].iloc[0]
    lines.append('D2 (CDR requires varying rho_c): %s  [beta_U0 %.3f vs beta_OBS %.3f+-%.3f; A0 C2 %s]'%('YES' if (not u.C2 and u.beta<bo[2]-2*np.hypot(u.ebeta,so[2]) and a.C2) else 'NO / not discriminating',u.beta,bo[2],so[2],a.C2))
    lines.append('external: Lelli+2016 beta 0.62+-0.04 (alpha=1 fixed); MOND low-Sigma slope 0.5')
    open('cdr_report_20261030.txt','w').write('\n'.join(lines)); print('\n'.join(lines))
