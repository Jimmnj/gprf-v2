# feature_transmission_20261109.py -- PREREG_FEATURE_TRANSMISSION_20261109 (STAMPED JM 16:27). Fixed-M/L arm; spherical pass2 solve_gal. [UNRATIFIED result]
import numpy as np, pandas as pd, io, contextlib, os, warnings, sys; warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
assert abs(T.loc['NGC2841'].Rd-3.64)<1e-6
fix=pd.read_csv('clean_conc_part_20261109.csv').set_index('gal').lrc_fix0.dropna()
def refit(g):
    d=P[P.gal==g].sort_values('r'); vo=d.vobs.values; w=1/(0.1*vo)**2; m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values
    if m.sum()<3: return np.nan
    chi=lambda l: (lambda o: 1e12 if o is None else float(np.sum((w*(np.sqrt(o.vmod2.values)-vo)**2)[m])))(solve_gal(g,l))
    return minimize_scalar(chi,bounds=(4.5,10.5),method='bounded',options={'xatol':0.02}).x
def detr(r,v):
    lr=np.log10(r); lv=np.log10(v); c=np.polyfit(lr,lv,2); return lv-np.polyval(c,lr)
rows=[]; named=['NGC6946','NGC2403','NGC3198','NGC7331','NGC2841','IC2574','DDO154','UGC06787','UGC02953']
for g in gals:
    d=P[P.gal==g].sort_values('r')
    if len(d)<8: continue
    r=d.r.values; vo=d.vobs.values; ev=d.ev.values; vb2=np.clip(d.vbar2.values,0,None)
    if (vb2<=0).any() or (vo<=0).any(): 
        m=(vb2>0)&(vo>0); 
        if m.sum()<8: continue
        d=d[m]; r,vo,ev,vb2=r[m],vo[m],ev[m],vb2[m]
    o1=solve_gal(g,7.45)
    if o1 is None: continue
    lrf=fix.get(g,np.nan)
    if not np.isfinite(lrf): lrf=refit(g)
    o1f=solve_gal(g,lrf) if np.isfinite(lrf) else None
    v1=np.sqrt(o1.vmod2.values[np.isin(o1.r.values,r)]); v1f=np.sqrt(o1f.vmod2.values[np.isin(o1f.r.values,r)]) if o1f is not None else v1*np.nan
    gb=vb2/r; vR=np.sqrt(gb/(1-np.exp(-np.sqrt(np.clip(gb/a0k,1e-12,None))))*r)
    vb=np.sqrt(vb2); fb,fo,f1,f1f,fR=[detr(r,v) for v in (vb,vo,v1,v1f,vR)]
    # window from baryons: largest excursion with separation 2..8 points
    best=None
    for i in range(len(r)):
        for j in range(len(r)):
            if 2<=abs(i-j)<=8:
                a=fb[i]-fb[j]
                if best is None or a>best[0]: best=(a,i,j)
    Ab,ip,it=best; lo,hi=min(ip,it),max(ip,it)
    n=np.median(2*ev[lo:hi+1]/vo[lo:hi+1]/np.log(10)); gw=np.mean(gb[lo:hi+1])/a0k
    Tx=lambda f: (f[ip]-f[it])/Ab
    rows.append(dict(gal=g,Q=int(d.Q.iloc[0]),npts=len(r),A_bar=Ab,noise=n,sig=Ab>=2*n and hi-lo+1>=3,r_pk=r[ip],r_tr=r[it],gw=gw,
                     gclass='HIGH' if gw>1 else ('MID' if gw>=0.3 else 'LOW'),T_obs=Tx(fo),T_1=Tx(f1),T_1f=Tx(f1f),T_R=Tx(fR),lrc_fit=lrf,named=g in named))
R=pd.DataFrame(rows); R.to_csv('feature_transmission_20261109.csv',index=False)
S_=R[R.sig]
print('galaxies %d | significant baryon feature %d | by class: %s'%(len(R),len(S_),S_.gclass.value_counts().to_dict()))
rng=np.random.default_rng(3)
def bs(x,f=np.median):
    x=np.asarray(x); x=x[np.isfinite(x)]; b=[f(x[rng.integers(0,len(x),len(x))]) for _ in range(1000)]; return f(x),np.percentile(b,16),np.percentile(b,84),len(x)
def report(D,lab):
    print('\n'+lab+' (n=%d)'%len(D))
    for c in ('T_obs','T_1','T_1f','T_R'):
        m,l,h,n=bs(D[c]); print('  %-5s median %+.2f [%+.2f,%+.2f]'%(c,m,l,h))
    for c,lab2 in (('T_1','D_1 = T_obs-T_1'),('T_1f','D_1f = T_obs-T_1f'),('T_R','D_R = T_obs-T_R')):
        m,l,h,n=bs(D.T_obs-D[c]); print('  %-18s median %+.2f [%+.2f,%+.2f]  %s'%(lab2,m,l,h,'excludes 0' if l>0 or h<0 else 'covers 0'))
report(S_,'ALL significant'); report(S_[S_.gclass!='HIGH'],'SCORED: LOW+MID significant'); report(S_[S_.gclass=='LOW'],'LOW only'); report(S_[S_.gclass=='HIGH'],'HIGH only')
print('\nNamed galaxies:'); print(R[R.named][['gal','Q','A_bar','noise','sig','r_pk','r_tr','gw','gclass','T_obs','T_1','T_1f','T_R']].round(3).to_string(index=False))
