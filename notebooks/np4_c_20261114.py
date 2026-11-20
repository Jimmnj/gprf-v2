# np4_c_20261114.py [UNRATIFIED] Item 4, part C: sec 3.3 second caution (disc M/L 0.7 and mass-dependent M/L) under the V3 first-minimum fitter.
# Rule 3a: test1d_yd_robust_20261018.run() is imported and used UNCHANGED; scipy.optimize.minimize_scalar is wrapped before the import (same wrapper logic
# as ep3_wrap_20261114.py: scan log rho_c 4.5-13.0 at 0.1, first interior local minimum, refine +-0.1; the wrapper returns V3 and records the Brent value).
# Reproduction: the recorded Brent value per galaxy must equal test1d_yd_<var>_20261018.csv l_rhoc. Resumable: python3 np4_c_20261114.py <var> <seconds>
import sys, os, time, numpy as np, pandas as pd, scipy.optimize as so
from multiprocessing import Pool
_orig=so.minimize_scalar; GR3=np.round(np.arange(4.5,13.0001,0.10),2); LAST={}
def _wrap(f,*a,**k):
    res=_orig(f,*a,**k); b=k.get('bounds')
    if b is None or not (b[0]>=4 and b[1]<=14 and b[1]-b[0]>=2): return res
    c=np.array([f(l) for l in GR3]); ks=[i for i in range(1,len(c)-1) if c[i]<c[i-1] and c[i]<=c[i+1]]
    if ks:
        i=ks[0]; r=_orig(f,bounds=(GR3[i]-0.1,GR3[i]+0.1),method='bounded',options={'xatol':0.005}); x,fx=(r.x,r.fun) if r.fun<=c[i] else (GR3[i],c[i])
    else: i=int(np.argmin(c)); x,fx=GR3[i],c[i]
    LAST.update(brent=res.x,chi_brent=res.fun,nmin=len(ks),edge=int(not ks),v3=x,chi_v3=fx); res.x=x; res.fun=fx; return res
so.minimize_scalar=_wrap
import test1d_yd_robust_20261018 as T1
def one(args):
    try: o=T1.run(args)
    except ValueError as e:   # run()'s Burkert fit to the model halo starts at l_rhoc and is bounded [3,11]: an edge fit above 11 cannot start it
        o=dict(gal=args[0],var=args[1],l_rhoc=LAST.get('v3'),chi=LAST.get('chi_v3'),err=str(e)[:40])
    o.update(LAST); return o
if __name__=='__main__':
    var=sys.argv[1]; budget=float(sys.argv[2]); out='np4_test1d_v3_%s_20261114_part.csv'%var; t0=time.time()
    done=set(pd.read_csv(out).gal) if os.path.exists(out) else set()
    todo=[g for g in T1.N14.index if g in set(T1.P0.gal) and g not in done]; print(var,'todo',len(todo))
    with Pool(2) as Pl:
        for o in Pl.imap_unordered(one,[(g,var) for g in todo]):
            cols=['gal','var','Yd','Mb','Q','xb0','l_rhoc','chi','lr0_m','lrho0_m','lr0_d','lrho0_d','chi_burk','rmax','brent','chi_brent','nmin','edge','v3','chi_v3','err']; pd.DataFrame([o]).reindex(columns=cols).to_csv(out,mode='a',header=not os.path.exists(out),index=False)
            if time.time()-t0>budget: Pl.terminate(); break
    print('rows',len(pd.read_csv(out)),'of',len([g for g in T1.N14.index if g in set(T1.P0.gal)]),'%.0fs'%(time.time()-t0))
