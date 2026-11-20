# errprop2_20261112.py -- PREREG_ERRPROP_20261112 sec 7 NEXT RUN (declared 22:01): perturbations with a GLOBAL fitter. [UNRATIFIED]
import sys, os, time, numpy as np, pandas as pd
os.environ['NOFIT']='1'
exec(open('errprop_globalfit_20261112.py').read())
GR=np.round(np.arange(4.5,10.5001,0.10),2)
def gfit2(f):
    c=np.array([f(l) for l in GR]); i=int(np.argmin(c)); lo,hi=max(4.5,GR[i]-0.15),min(10.5,GR[i]+0.15)
    r=minimize_scalar(f,bounds=(lo,hi),method='bounded',options={'xatol':0.005}); return r.x if r.fun<=c[i] else GR[i]
base=pd.read_csv('errprop_globalfit_20261112_part.csv').drop_duplicates('gal').set_index('gal').lrc_g
PART='errprop2_20261112_part.csv'; done=set(pd.read_csv(PART).gal) if os.path.exists(PART) else set(); t0=time.time()
for g in gals:
    if g in done: continue
    if time.time()-t0>150: break
    t=T0.loc[g]; eD=max(t.eD,1e-3*t.D); bul=bool((P0[P0.gal==g].vb>0).any()); i_alt=t.inc+t.einc if t.inc+t.einc<=90 else t.inc-t.einc
    row=dict(gal=g,lrc0=base[g],lrcDp=gfit2(chi_fn(g,Dp=1+eD/t.D)),lrcDm=gfit2(chi_fn(g,Dp=max(1-eD/t.D,0.2))),lrcI=gfit2(chi_fn(g,inc=i_alt)),
             lrcU=gfit2(chi_fn(g,ud=0.5*10**0.1)),lrcUb=gfit2(chi_fn(g,ub=0.7*10**0.1)) if bul else np.nan,D=t.D,eD=eD,inc=t.inc,einc=t.einc,bul=int(bul))
    pd.DataFrame([row]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False)
print('progress %d / %d (%.0f s)'%(len(pd.read_csv(PART)),len(gals),time.time()-t0))
