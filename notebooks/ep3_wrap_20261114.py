# ep3_wrap_20261114.py -- EP-3 runner [UNRATIFIED]: executes an existing script UNCHANGED with scipy.optimize.minimize_scalar wrapped.
# For every bounded fit whose bounds lie within [4,14] (the log rho_c fits), the wrapper also scans log rho_c 4.5-13.0 at 0.1 with the SAME
# objective, finds the first (lowest-rho_c) interior local minimum, refines +-0.1 (V3, BD-2 (b)), and prints: Brent result vs V3.
# With EP3_USE_V3=1 the wrapper RETURNS the V3 fit, so the script's downstream numbers are the V3 numbers.
# Run only in a scratch copy of notebooks (the scripts write their usual output files). Usage: python3 ep3_wrap_20261114.py <script> [args]
import sys, os, numpy as np, scipy.optimize as so
_orig=so.minimize_scalar; GR3=np.round(np.arange(4.5,13.0001,0.10),2); USE=os.environ.get('EP3_USE_V3')=='1'
def _wrap(f,*a,**k):
    res=_orig(f,*a,**k); b=k.get('bounds')
    if b is None or not (b[0]>=4 and b[1]<=14 and b[1]-b[0]>=2): return res
    c=np.array([f(l) for l in GR3]); ks=[i for i in range(1,len(c)-1) if c[i]<c[i-1] and c[i]<=c[i+1]]
    if ks:
        i=ks[0]; r=_orig(f,bounds=(GR3[i]-0.1,GR3[i]+0.1),method='bounded',options={'xatol':0.005}); x,fx=(r.x,r.fun) if r.fun<=c[i] else (GR3[i],c[i])
    else: i=int(np.argmin(c)); x,fx=GR3[i],c[i]
    print('EP3FIT brent %.3f chi %.4g | v3 %.3f chi %.4g | nmin %d edge %d | diff %+.3f'%(res.x,res.fun,x,fx,len(ks),int(not ks),x-res.x),flush=True)
    if USE: res.x=x; res.fun=fx
    return res
so.minimize_scalar=_wrap
script=sys.argv[1]; sys.argv=sys.argv[1:]
exec(compile(open(script).read(),script,'exec'),globals())   # in the real __main__ so multiprocessing can pickle the script's functions
