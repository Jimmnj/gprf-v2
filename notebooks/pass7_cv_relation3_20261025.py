# pass7_cv_relation3_20261025.py [UNRATIFIED] 2026-10-25 00:45. Sim-review item C1 (simulated Desmond; test ours): relation (3) scored out of sample.
# Per-galaxy fitted log rho_c (full 147-galaxy sample, pass2_table2_part.csv) regressed on [1, log M_b, log Sigma_h (Msun/pc^2), f_gas].
# (a) in-sample refit on all 147; (b) K-fold by galaxy: fit on train, predict test rho_c, solve the law, score per-point rms of log g.
# Usage: python3 pass7_cv_relation3_20261025.py <seed> <K>   (appends to pass7_cv_results.csv)
import numpy as np, pandas as pd, io, contextlib, sys, time, os, warnings; warnings.filterwarnings('ignore')
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
seed=int(sys.argv[1]); K=int(sys.argv[2])
F=pd.read_csv('pass2_table2_part.csv')[['gal','lrc']].set_index('gal')
def feats(g):
    t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9; Sh=Mb/(2*np.pi*(1.678*t.Rd*1e3)**2) if t.Rd>0 else np.nan
    return [1.0,np.log10(Mb),np.log10(Sh),1.33*t.MHI*1e9/Mb]
GL=[g for g in F.index if np.isfinite(feats(g)[2])]
X=np.array([feats(g) for g in GL]); y=F.loc[GL,'lrc'].values
def score(pred):   # pred: dict gal -> lrc ; returns per-point residual array
    res=[]
    for g,l in pred.items():
        o=solve_gal(g,l)
        if o is None: continue
        o=o[(o.ev/o.vobs<0.1)&(o.vobs>0)&(o.vbar2>0)]
        res.append(np.log10(o.vobs**2/o.r)-np.log10(o.vmod2/o.r))
    return np.concatenate(res)
t0=time.time(); rows=[]
if seed==0:
    b=np.linalg.lstsq(X,y,rcond=None)[0]
    r=score({g:float(X[i]@b) for i,g in enumerate(GL)})
    rows.append(dict(kind='insample_refit',seed=0,K=0,rms=np.sqrt(np.mean(r**2)),N=len(r),b0=b[0],b1=b[1],b2=b[2],b3=b[3]))
else:
    rng=np.random.default_rng(seed); idx=rng.permutation(len(GL)); folds=np.array_split(idx,K); pred={}
    for f in folds:
        tr=np.setdiff1d(idx,f); b=np.linalg.lstsq(X[tr],y[tr],rcond=None)[0]
        for i in f: pred[GL[i]]=float(X[i]@b)
    r=score(pred); rows.append(dict(kind='cv',seed=seed,K=K,rms=np.sqrt(np.mean(r**2)),N=len(r)))
out='pass7_cv_results.csv'; pd.DataFrame(rows).to_csv(out,mode='a',header=not os.path.exists(out),index=False)
print(rows, '%.0fs'%(time.time()-t0))
