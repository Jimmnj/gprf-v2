# errprop_20261112.py -- PREREG_ERRPROP_20261112.md (written 21:51 before any perturbed fit). [UNRATIFIED]
# Solver: pass2_sparc_full_20261022.py solve_gal VERBATIM (exec'd); fit: pass2_table2_full_20261022.py law_chi VERBATIM; nuisances = input rescalings only.
# Usage: python3 errprop_20261112.py [--validate N] [--budget SECONDS]; resumable via errprop_20261112_part.csv.
import numpy as np, pandas as pd, os, sys, io, contextlib, time, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar
ns={}
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src,ns)
P0=ns['P'].copy(); T0=ns['T'].copy(); solve_gal=ns['solve_gal']
Lsel=(P0.Q<=2)&(P0.inc>=30)&(P0.ev/P0.vobs<0.1)&(P0.vobs>0)&(P0.vbar2>0)
gals=sorted(P0[Lsel].gal.unique())
def fit(g,Dp=1.0,inc=None,ud=0.5,ub=0.7):
    d0=P0[P0.gal==g].sort_values('r'); t0=T0.loc[g]; i0=t0.inc; inc=i0 if inc is None else inc
    m=((d0.ev/d0.vobs<0.1)&(d0.vobs>0)&(d0.vbar2>0)).values          # Lelli mask at nominal values (as marg_fit)
    d=d0.copy(); d['r']=d0.r*Dp
    d['vbar2']=(np.sign(d0.vg)*d0.vg**2+ud*d0.vd**2+ub*d0.vb**2)*Dp
    d['vobs']=d0.vobs*np.sin(np.radians(i0))/np.sin(np.radians(inc))
    t=t0.copy(); t['L']=t0.L*(ud/0.5)*Dp**2; t['MHI']=t0.MHI*Dp**2; t['Rd']=t0.Rd*Dp
    ns['P']=d; ns['T']=pd.DataFrame([t],index=[g])
    vo=d.vobs.values; w=1/(0.1*vo)**2
    def law_chi(l):
        o=solve_gal(g,l)
        if o is None: return 1e12
        return float(np.sum((w*(np.sqrt(o.vmod2.values)-vo)**2)[m]))
    return minimize_scalar(law_chi,bounds=(4.5,10.5),method='bounded',options={'xatol':0.02}).x
args=sys.argv[1:]; budget=float(args[args.index('--budget')+1]) if '--budget' in args else 160
if '--validate' in args:
    n=int(args[args.index('--validate')+1]); ref=pd.read_csv('pass2_table2_part.csv').set_index('gal').lrc
    t0=time.time()
    for g in gals[:n]:
        a=time.time(); x=fit(g); print('%-10s lrc %.3f  pass2 %.3f  diff %+.3f  (%.1f s)'%(g,x,ref[g],x-ref[g],time.time()-a),flush=True)
    print('validate: %.1f s for %d galaxies'%(time.time()-t0,n)); sys.exit()
PART='errprop_20261112_part.csv'; done=set(pd.read_csv(PART).gal) if os.path.exists(PART) else set(); tstart=time.time()
for g in gals:
    if g in done: continue
    if time.time()-tstart>budget: break
    t=T0.loc[g]; eD=max(t.eD,1e-3*t.D); bul=bool((P0[P0.gal==g].vb>0).any())
    i_alt=t.inc+t.einc if t.inc+t.einc<=90 else t.inc-t.einc
    row=dict(gal=g,lrc0=fit(g),lrcDp=fit(g,Dp=1+eD/t.D),lrcDm=fit(g,Dp=max(1-eD/t.D,0.2)),lrcI=fit(g,inc=i_alt),lrcU=fit(g,ud=0.5*10**0.1),
             lrcUb=fit(g,ub=0.7*10**0.1) if bul else np.nan,D=t.D,eD=eD,inc=t.inc,einc=t.einc,bul=int(bul))
    pd.DataFrame([row]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False)
n=len(pd.read_csv(PART)) if os.path.exists(PART) else 0; print('progress %d / %d (%.0f s this call)'%(n,len(gals),time.time()-tstart))
