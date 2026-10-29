# gate2c_robustness_20260925.py  [UNRATIFIED]
# Rebuild of Gate 2c (RAR residual vs r/Rd) + check 4 (beam smearing / inner points)
# + check 5 (nonparametric g_bar control) + galaxy-clustered inference.
# Setup follows RESULT_RAR_RESIDUALS_20260923 sec 1. Paths: set SPARC_DIR, MRT.
import numpy as np, pandas as pd, glob, os, warnings; warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
SPARC_DIR=os.path.join('..','sparc_galaxy_data','sparc_database'); MRT=os.path.join('..','sparc_galaxy_data','SPARC_Lelli2016c.mrt')  # portable (Windows/WSL) 9/30  # unzip sparc_database.zip first; run from notebooks\
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6
rng=np.random.default_rng(1)
def load_table():
    rows=[]
    for line in open(MRT):
        p=line.split()
        if len(p)>=18 and p[1].isdigit():
            try: rows.append(dict(gal=p[0],D=float(p[2]),inc=float(p[5]),L=float(p[7]),Reff=float(p[9]),
                    Rd=float(p[11]),MHI=float(p[13]),RHI=float(p[14]),Q=int(p[17])))
            except ValueError: pass
    return pd.DataFrame(rows).set_index('gal')
def load(Yd=0.5,Yb=0.7):
    t=load_table(); out=[]
    for f in glob.glob(os.path.join(SPARC_DIR,'*_rotmod.dat')):
        g=os.path.basename(f).replace('_rotmod.dat','')
        if g not in t.index: continue
        d=np.atleast_2d(np.loadtxt(f,comments='#')); r=d[:,0]
        df=pd.DataFrame(dict(gal=g,r=r,vobs=d[:,1],ev=d[:,2],vg=d[:,3],vd=d[:,4],vb=d[:,5]))
        df['idx']=np.arange(len(r)); df['npts']=len(r)
        df['dr']=np.median(np.diff(r)) if len(r)>1 else r[0]; df['rmax']=r.max(); out.append(df)
    P=pd.concat(out).join(t,on='gal')
    P['vbar2']=np.sign(P.vg)*P.vg**2+Yd*P.vd**2+Yb*P.vb**2
    P=P[(P.Q<=2)&(P.inc>30)&(P.vbar2>0)&(P.vobs>0)].reset_index(drop=True)
    P['gbar']=P.vbar2/P.r*1e6/kpc; P['lgb']=np.log10(P.gbar)
    P['lgo']=np.log10(P.vobs**2/P.r*1e6/kpc)
    P['lgmond']=np.log10(P.gbar/(1-np.exp(-np.sqrt(P.gbar/a0))))
    P['res']=P.lgo-P.lgmond
    P['lr']=np.log10(P.r); P['lRd']=np.log10(P.Rd); P['lrRd']=P.lr-P.lRd
    P['theta']=P.r/(P.D*1e3)*206265; P['rdr']=P.r/P.dr     # smearing proxies
    P['elg']=2*P.ev/P.vobs/np.log(10)
    Mb=0.5*P.L*1e9+1.33*P.MHI*1e9; vinf2=np.sqrt(a0k*G*Mb)
    vh2=vinf2*(1-(P.Rd/P.r)*np.arctan(P.r/P.Rd))
    P['pred']=np.log10((P.vbar2+vh2)/P.r*1e6/kpc)-P.lgmond      # GPRF-predicted MOND residual
    P['lrRmax']=np.log10(P.r/P.rmax); P['lrReff']=np.log10(P.r/P.Reff)
    P['lrRHI']=np.log10(P.r/P.RHI.where(P.RHI>0)); P['fidx']=(P.idx+0.5)/P.npts
    return P
def np_detrend(df,y='res',w=0.1):
    x=df.lgb.values; v=df[y].values; o=np.argsort(x); xs,vs=x[o],v[o]
    m=np.array([np.median(vs[np.searchsorted(xs,xi-w/2):np.searchsorted(xs,xi+w/2)]) for xi in x])
    return v-m
def basis(x,nk=12):
    k=np.quantile(x,np.linspace(0,1,nk+2)[1:-1])
    return np.column_stack([np.ones_like(x),x,x**2,x**3]+[np.clip(x-q,0,None)**3 for q in k])
def ols(df,cov,y='res'):
    X=np.column_stack([basis(df.lgb.values)]+[df[c].values for c in cov])
    return np.linalg.lstsq(X,df[y].values,rcond=None)[0][-len(cov):]
def boot(df,cov,y='res',nb=400):
    df=df.dropna(subset=cov).reset_index(drop=True)
    gals=df.gal.unique(); grp={g:df.index[df.gal==g] for g in gals}; est=ols(df,cov,y); bs=[]
    for _ in range(nb):
        idx=np.concatenate([grp[g] for g in rng.choice(gals,len(gals))])
        bs.append(ols(df.loc[idx].reset_index(drop=True),cov,y))
    bs=np.array(bs); return est,*np.percentile(bs,[2.5,97.5],axis=0),len(df),len(gals)
def line(name,d,cov=['lrRd']):
    d=d.reset_index(drop=True); e,lo,hi,n,g=boot(d,cov); p=ols(d,cov,'pred')
    print(f'{name:30s} N={n:4d} g={g:3d} '+' '.join(f'{c}={a:+.3f}[{b:+.3f},{c2:+.3f}] pred={q:+.3f}' for c,a,b,c2,q in zip(cov,e,lo,hi,p)))
if __name__=='__main__':
    P=load(); print('N',len(P),'gal',P.gal.nunique())
    print('\n[5] nonparametric g_bar control (running median, 0.1 dex):')
    rn=np_detrend(P); print({v:round(spearmanr(rn,P[v])[0],3) for v in ['lrRd','lr','lRd']})
    print('\n[scale horse race] rho vs r/scale (np-detrended, RHI>0 subset):')
    S=P[P.RHI>0].reset_index(drop=True); rn=np_detrend(S)
    print({v:round(spearmanr(rn,S[v])[0],3) for v in ['lrRd','lrReff','lrRHI','lrRmax','fidx','lr']})
    line('joint r/Rd + r/Rmax',S,['lrRd','lrRmax']); line('joint r/Rd + fidx',S,['lrRd','fidx'])
    print('\n[4] smearing / quality cuts, spline(g_bar) control, galaxy-bootstrap 95% CI:')
    Q1=P.Q==1; I=P.inc<=80
    for k,m in {'ALL':P.Q>0,'Q=1':Q1,'Q=2':P.Q==2,'theta>=60"':P.theta>=60,'drop 3 innermost':P.idx>=3,
                'r>=4 steps':P.rdr>=4,'Q=1 & r>=2 steps':Q1&(P.rdr>=2),'Q=1 & r>=4 steps':Q1&(P.rdr>=4),
                'Q=1 & inc<=80 & r>=4':Q1&I&(P.rdr>=4),'Q=1 & theta>=60"':Q1&(P.theta>=60),
                'Q=1 innermost r<2 steps':Q1&(P.rdr<2),'inner half':P.idx<P.npts/2,'outer half':P.idx>=P.npts/2}.items():
        line(k,P[m])
