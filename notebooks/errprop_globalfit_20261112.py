# errprop_globalfit_20261112.py -- PREREG_ERRPROP_20261112 sec 7 (declared 21:55): GLOBAL fixed-arm fit for all 147 galaxies. [UNRATIFIED]
# Solver + chi^2 verbatim (errprop_20261112.py internals = pass2 solve_gal + pass2_table2 law_chi). Global fit = 0.05-dex grid on [4.5,10.5] then bounded refine +-0.1.
import sys, os, time, numpy as np, pandas as pd
sys.argv=['x','--validate','0']
exec(open('errprop_20261112.py').read().split("args=sys.argv[1:]")[0])
def chi_fn(g,Dp=1.0,inc=None,ud=0.5,ub=0.7):
    d0=P0[P0.gal==g].sort_values('r'); t0=T0.loc[g]; i0=t0.inc; inc=i0 if inc is None else inc
    m=((d0.ev/d0.vobs<0.1)&(d0.vobs>0)&(d0.vbar2>0)).values
    d=d0.copy(); d['r']=d0.r*Dp; d['vbar2']=(np.sign(d0.vg)*d0.vg**2+ud*d0.vd**2+ub*d0.vb**2)*Dp
    d['vobs']=d0.vobs*np.sin(np.radians(i0))/np.sin(np.radians(inc))
    t=t0.copy(); t['L']=t0.L*(ud/0.5)*Dp**2; t['MHI']=t0.MHI*Dp**2; t['Rd']=t0.Rd*Dp
    vo=d.vobs.values; w=1/(0.1*vo)**2
    def f(l):
        ns['P']=d; ns['T']=pd.DataFrame([t],index=[g]); o=solve_gal(g,l)
        return 1e12 if o is None else float(np.sum((w*(np.sqrt(o.vmod2.values)-vo)**2)[m]))
    return f
GRID=np.round(np.arange(4.5,10.5001,0.05),2)
def gfit(f):
    c=np.array([f(l) for l in GRID]); i=int(np.argmin(c)); lo,hi=max(4.5,GRID[i]-0.1),min(10.5,GRID[i]+0.1)
    r=minimize_scalar(f,bounds=(lo,hi),method='bounded',options={'xatol':0.005}); x,cx=(r.x,r.fun) if r.fun<=c[i] else (GRID[i],c[i])
    nmin=sum(1 for k in range(1,len(c)-1) if c[k]<c[k-1] and c[k]<=c[k+1])
    return x,cx,nmin
if __name__=='__main__' and 'NOFIT' not in os.environ:
    PART='errprop_globalfit_20261112_part.csv'; done=set(pd.read_csv(PART).gal) if os.path.exists(PART) else set(); t0=time.time()
    for g in gals:
        if g in done: continue
        if time.time()-t0>150: break
        x,cx,nm=gfit(chi_fn(g)); pd.DataFrame([dict(gal=g,lrc_g=x,chi_g=cx,n_localmin=nm)]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False)
    print('progress %d / %d (%.0f s)'%(len(pd.read_csv(PART)),len(gals),time.time()-t0))
