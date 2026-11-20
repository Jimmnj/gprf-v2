# errprop3_first_20261113.py -- BD-3 (JM 2026-11-13 19:42 BD-2 = (b)): redo ERRPROP perturbed fits for the 12 multi-minimum
# galaxies under the V3 FIRST-MINIMUM definition of rho_c. [UNRATIFIED]
# Solver + chi^2 verbatim via errprop_globalfit_20261112.chi_fn (rule 3a; reproduced 147/147 in run 1).
# Fitter V3: scan log rho_c 4.5-13.0 at 0.1; take the FIRST (lowest-rho_c) interior local minimum; refine +-0.1 (bounded Brent).
#   If no interior local minimum exists, take the grid argmin and set edge=1.
# Perturbations identical to errprop2_20261112.py (D +-eD, inc +-einc, disc M/L x10^0.1, bulge M/L x10^0.1).
# FIX 19:55: first run wrote 12- and 13-column rows (lrcUb only for bulge galaxies) -> CSV unreadable; lrcUb now always written (NaN if no bulge). Bad file moved to _to_delete/.
# EXPECTATIONS (written 2026-11-13 ~19:50, before any run):
#  X1 baseline reproduces bound_check_first_20261112.csv lrc_first for all 12 to < 0.01 dex.
#  X2 NGC 0289 sigma_tot falls from 5.3 dex to < 1 dex (no branch jumps under V3).
#  X3 tau (unflagged, all) moves by <= 0.05 from 0.48; median unchanged to 0.02; direct/Hubble split ordering unchanged (0.53 vs 0.37).
#  X4 no galaxy outside the 12 needs V3 (no perturbed fit jumps > 2.5 dex from its baseline in errprop2).
import sys, os, time, numpy as np, pandas as pd
EXTRA=sys.argv[1].split(',') if len(sys.argv)>1 else []   # captured before exec (errprop_20261112 overwrites sys.argv)
os.environ['NOFIT']='1'
exec(open('errprop_globalfit_20261112.py').read())
GR3=np.round(np.arange(4.5,13.0001,0.10),2)
def ffit(f):
    c=np.array([f(l) for l in GR3]); ks=[k for k in range(1,len(c)-1) if c[k]<c[k-1] and c[k]<=c[k+1]]
    if not ks: i=int(np.argmin(c)); return GR3[i],1
    k=ks[0]; r=minimize_scalar(f,bounds=(GR3[k]-0.1,GR3[k]+0.1),method='bounded',options={'xatol':0.005})
    return (r.x if r.fun<=c[k] else GR3[k]),0
SEL=list(pd.read_csv('bound_check_first_20261112.csv').gal)+EXTRA  # X4 extension 20:00: extra galaxies via argv
PART='errprop3_first_20261113_part.csv'; done=set(pd.read_csv(PART).gal) if os.path.exists(PART) else set(); t0=time.time()
for g in SEL:
    if g in done: continue
    if time.time()-t0>140: break
    t=T0.loc[g]; eD=max(t.eD,1e-3*t.D); bul=bool((P0[P0.gal==g].vb>0).any()); i_alt=t.inc+t.einc if t.inc+t.einc<=90 else t.inc-t.einc
    fits={'lrc0':ffit(chi_fn(g)),'lrcDp':ffit(chi_fn(g,Dp=1+eD/t.D)),'lrcDm':ffit(chi_fn(g,Dp=max(1-eD/t.D,0.2))),
          'lrcI':ffit(chi_fn(g,inc=i_alt)),'lrcU':ffit(chi_fn(g,ud=0.5*10**0.1))}
    fits['lrcUb']=ffit(chi_fn(g,ub=0.7*10**0.1)) if bul else (np.nan,0)
    row=dict(gal=g,**{k:v[0] for k,v in fits.items()},edge=int(any(v[1] for v in fits.values())),D=t.D,eD=eD,inc=t.inc,einc=t.einc,bul=int(bul))
    pd.DataFrame([row]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False)
print('progress %d / %d (%.0f s)'%(len(pd.read_csv(PART)),len(SEL),time.time()-t0))
