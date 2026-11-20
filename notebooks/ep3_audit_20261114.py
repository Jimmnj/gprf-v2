# ep3_audit_20261114.py -- EP-3 FIT-1 audit [UNRATIFIED]. JM 2026-11-14 01:39 "go ep3".
# Every bounded-Brent fit of log rho_c that feeds a number in the note, re-fitted with the V3 first-minimum fitter (BD-2 (b)).
# Rule 3a: canon sec 6 lines read at boot 01:34-01:40; solver = pass2_sparc_full_20261022.py solve_gal, exec'd verbatim; ffit verbatim from
# errprop3_first_20261113.py. Validation = reproduce pass2_table2_part.csv (V3, md5 82dce260) for all 147.
# EXPECTATIONS (written before any run of this script, between JM 01:39 "go ep3" and the first run; the stamp "01:58" first put here was a guess, not a clock reading -- corrected 02:04):
#  E0 validation: max |ffit - table| < 0.01 dex over 147.
#  E1 clean_conc lrc_fix0 differs from V3 by > 0.05 only for galaxies among NGC4013, NGC5907, UGC01281 that are in its set.
#  E2 clean_conc lrc_fix2 (inner two points dropped): <= 3 galaxies differ by > 0.05.
#  E3 sigmafree lrc0: same three at most; h2_refit: 0 of 18 differ.
# Usage: python3 ep3_audit_20261114.py <stage>   stage in val, cc, sf, h2. Resumable; outputs ep3_<stage>_20261114_part.csv.
import sys, os, time, io, contextlib, json, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar
STAGE=sys.argv[1]
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
GR3=np.round(np.arange(4.5,13.0001,0.10),2)
def ffit(f):
    c=np.array([f(l) for l in GR3]); ks=[k for k in range(1,len(c)-1) if c[k]<c[k-1] and c[k]<=c[k+1]]
    if not ks: i=int(np.argmin(c)); return GR3[i],1,c[i],0
    k=ks[0]; r=minimize_scalar(f,bounds=(GR3[k]-0.1,GR3[k]+0.1),method='bounded',options={'xatol':0.005})
    return (r.x if r.fun<=c[k] else GR3[k]),0,min(r.fun,c[k]),len(ks)
def chi_drop(g,drop):   # verbatim from clean_conc_20261109.refit
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vo=d.vobs.values; w=1/(0.1*vo)**2
    m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values.copy(); m[:drop]=False
    if m.sum()<3: return None
    def chi(l):
        o=solve_gal(g,l); return 1e12 if o is None else float(np.sum((w*(np.sqrt(o.vmod2.values)-vo)**2)[m]))
    return chi
PART='ep3_%s_20261114_part.csv'%STAGE; done=set(pd.read_csv(PART).key) if os.path.exists(PART) else set(); t0=time.time()
def put(row):
    pd.DataFrame([row]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False)
if STAGE=='val':
    Tb=pd.read_csv('pass2_table2_part.csv').set_index('gal')
    for g in Tb.index:
        if g in done: continue
        if time.time()-t0>140: break
        f=chi_drop(g,0); x,e,c,n=ffit(f) if f is not None else (np.nan,0,np.nan,0); put(dict(key=g,lrc_v3=x,edge=e,chi=c,nmin=n,lrc_tab=Tb.loc[g,'lrc'],chi_tab=Tb.loc[g,'chi_law']))
if STAGE=='cc':
    C=pd.read_csv('clean_conc_part_20261109.csv').drop_duplicates('gal').set_index('gal')
    for g in C.index:
        if g in done: continue
        if time.time()-t0>140: break
        row=dict(key=g,fix0_old=C.loc[g,'lrc_fix0'],fix2_old=C.loc[g,'lrc_fix2'])
        for d in (0,2):
            f=chi_drop(g,d)
            if f is None: row.update({'fix%d_v3'%d:np.nan,'edge%d'%d:0,'nmin%d'%d:0}); continue
            x,e,c,n=ffit(f); row.update({'fix%d_v3'%d:x,'edge%d'%d:e,'nmin%d'%d:n})
        put(row)
R=pd.read_csv(PART); print('stage %s progress %d (%.0f s)'%(STAGE,len(R),time.time()-t0))
