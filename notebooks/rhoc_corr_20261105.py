# rhoc_corr_20261105.py [UNRATIFIED] -- PREREG_RHOC_CORR_20261105
import numpy as np, pandas as pd, glob, os, io, contextlib, warnings; warnings.filterwarnings('ignore')
from scipy.stats import spearmanr, mannwhitneyu
rng=np.random.default_rng(5)
G=4.30091e-6; a0k=1.2e-10*3.0857e19/1e6
T=[]
for line in open('../sparc_galaxy_data/SPARC_Lelli2016c.mrt'):
    p=line.split()
    if len(p)>=18 and p[1].isdigit():
        try: T.append(dict(gal=p[0],T=int(p[1]),D=float(p[2]),fD=int(p[4]),inc=float(p[5]),L=float(p[7]),Reff=float(p[8]),SBeff=float(p[9]),Rd=float(p[11]),SBdisk=float(p[12]),MHI=float(p[13]),RHI=float(p[14]),Vf=float(p[15]),Q=int(p[17])))
        except ValueError: pass
T=pd.DataFrame(T).set_index('gal')
rows=[]
for f in glob.glob('../sparc_galaxy_data/sparc_database/*_rotmod.dat'):
    g=os.path.basename(f).replace('_rotmod.dat','')
    if g not in T.index: continue
    d=np.atleast_2d(np.loadtxt(f,comments='#')); r,vo,ev,vg,vd,vb=d[:,0],d[:,1],d[:,2],d[:,3],d[:,4],d[:,5]
    t=T.loc[g]; i=np.argmin(abs(r-t.Rd)) if t.Rd>0 else 0
    tot=lambda k: np.sign(vg[k])*vg[k]**2+0.5*vd[k]**2+0.7*vb[k]**2
    rows.append(dict(gal=g,BT=0.7*vb[-1]**2/max(tot(-1),1e-9),fgas_in=np.sign(vg[i])*vg[i]**2/max(tot(i),1e-9),rlast_Rd=r[-1]/t.Rd if t.Rd>0 else np.nan,npts=len(r),errfloor=np.median(ev/np.maximum(vo,1))))
R=pd.DataFrame(rows).set_index('gal')
fit=pd.read_csv('pass2_table2_part.csv').set_index('gal')[['lrc']]
mg=pd.read_csv('marg_fit_params_20261030.csv'); mg=mg[mg.arm=='A'].set_index('gal')[['lrc']].rename(columns={'lrc':'lrc_marg'})
D=fit.join(T).join(R).join(mg,how='left')
D['Mb']=0.5*D.L*1e9+1.33*D.MHI*1e9; D['lMb']=np.log10(D.Mb); D['fgas']=1.33*D.MHI*1e9/D.Mb; D['lSig']=np.log10(D.Mb/(2*np.pi*(1.678*D.Rd*1e3)**2))
D['btfr_res']=np.where(D.Vf>0,np.log10(D.Vf.clip(1e-3))-0.25*np.log10(G*a0k*D.Mb),np.nan)
D['RHI_Rd']=D.RHI/D.Rd; D['lSBdisk']=np.log10(D.SBdisk.clip(1e-3)); D['lD']=np.log10(D.D)
V=['btfr_res','BT','T','RHI_Rd','fgas_in','lSBdisk','rlast_Rd','npts','inc','Q','lD','errfloor']
def resid(y,X):
    X=np.column_stack([np.ones(len(y))]+[X[c].values for c in X.columns]); return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def run(ycol,label):
    print('\n=== %s (N=%d) ==='%(label,D[ycol].notna().sum())); out=[]
    for v in V:
        d=D[[ycol,v,'lMb','lSig','fgas']].dropna(); n=len(d)
        r,p=spearmanr(d[ycol],d[v]); b=[spearmanr(*(lambda s:(s[ycol],s[v]))(d.iloc[rng.integers(0,n,n)]))[0] for _ in range(600)]
        ry=resid(d[ycol].values,d[['lMb','lSig','fgas']]); rv=resid(d[v].values.astype(float),d[['lMb','lSig','fgas']]); rp,pp=spearmanr(ry,rv)
        out.append(dict(var=v,N=n,rho=r,lo=np.percentile(b,2.5),hi=np.percentile(b,97.5),p=p,rho_partial=rp,p_partial=pp))
    d=D[[ycol,'fD']].dropna(); hf=d[d.fD==1][ycol]; dr=d[d.fD>1][ycol]; u,p=mannwhitneyu(hf,dr)
    out.append(dict(var='fD (HubbleFlow vs direct): median diff',N=len(d),rho=hf.median()-dr.median(),lo=np.nan,hi=np.nan,p=p,rho_partial=np.nan,p_partial=np.nan))
    o=pd.DataFrame(out); m=len(o); o=o.sort_values('p'); o['holm_alpha']=0.05/(m-np.arange(m)); o['holm_sig']=False
    for k in range(m):
        if o.p.iloc[k]<o.holm_alpha.iloc[k]: o.iloc[k,o.columns.get_loc('holm_sig')]=True
        else: break
    pd.set_option('display.width',200); print(o.round(3).to_string(index=False)); return o
o1=run('lrc','PRIMARY: fixed-M/L rho_c'); o2=run('lrc_marg','S1: marginalised (MAP) rho_c')
o1.to_csv('rhoc_corr_20261105_primary.csv',index=False); o2.to_csv('rhoc_corr_20261105_marg.csv',index=False); D.to_csv('rhoc_corr_20261105_table.csv')
