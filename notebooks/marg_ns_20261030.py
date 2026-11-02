# marg_ns_20261030.py -- nested-sampling (dynesty) version of the marginalised fit: posteriors AND evidences per galaxy, per arm, per intrinsic-scatter value.
# Pre-registration: PREREG_NS_20261030.md (status in its header). Same harness as marg_fit_20261030.py (imports it). [UNRATIFIED]
# Arms: U relation at universal rho_c 10^7.45 | B RAR function | A relation with log rho_c free (flat [6,10]).
# Nuisances (all arms): log Ud ~ N(log 0.5, 0.1), log Ub ~ N(log 0.7, 0.1) (bulge), D ~ N(D0, eD) truncated at 0.2 D0, i ~ N(i0, ei) truncated [20, 90].
# Likelihood: Gaussian in v with sigma_v^2 = ev^2 + (sig_int * v * ln10 / 2)^2, i.e. an intrinsic scatter sig_int in log g added in quadrature; sig_int fixed per run from SIG_INT list.
# Output: one JSON line per (galaxy, arm, sig_int) in marg_ns_20261030.jsonl (resumable; finished keys are skipped). Fields: logZ, logZerr, posterior medians/16/84 of each parameter, posterior-median residual rms, n_like, seconds.
# Usage: python3 marg_ns_20261030.py --arms UBA --sig 0,0.03,0.06 --procs 30 [--slice a:b] [--pilot] [--nlive 300]
import os, sys, json, time, argparse, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.special import erfinv
import dynesty
sys.argv_backup=sys.argv; sys.argv=['x']
exec(open('marg_fit_20261030.py').read().split("if __name__")[0])   # Gal, load, init, constants
sys.argv=sys.argv_backup
LN10=np.log(10)
def ptform_factory(s,arm):
    def ptform(u):
        k=0; x=[]
        if arm=='A': x.append(6+4*u[k]); k+=1
        x.append(LUD0+SIGU*np.sqrt(2)*erfinv(2*u[k]-1)); k+=1
        if s.bul: x.append(LUB0+SIGU*np.sqrt(2)*erfinv(2*u[k]-1)); k+=1
        zD=np.sqrt(2)*erfinv(2*u[k]-1); k+=1
        zi=np.sqrt(2)*erfinv(2*u[k]-1); k+=1
        return np.array(x+[zD,zi])
    return ptform
def loglike_factory(s,arm,sig_int):
    def loglike(x):
        lrc,lud,lub,D,inc=s.unpack(x,arm)
        if D<0.2*s.D0: return -1e100
        m=s.model(x,arm)
        if m is None: return -1e100
        r,vb2,vobs,ev,Mb,vm2=m; vm=np.sqrt(np.clip(vm2,1,None))
        var=ev**2+(sig_int*vobs*LN10/2)**2
        ll=-0.5*np.sum((((vm-vobs)**2/var)+np.log(2*np.pi*var))[s.sel])
        return ll
    return loglike
def run_one(args):
    g,arm,sig_int,nlive=args; t0=time.time(); s=Gal(g)
    ndim=(1 if arm=='A' else 0)+1+(1 if s.bul else 0)+2
    rng=np.random.default_rng(abs(hash((g,arm,sig_int)))%2**32)
    sam=dynesty.NestedSampler(loglike_factory(s,arm,sig_int),ptform_factory(s,arm),ndim,nlive=nlive,rstate=rng,bound='multi',sample='rwalk')
    sam.run_nested(dlogz=0.1,print_progress=False)
    res=sam.results; w=np.exp(res.logwt-res.logz[-1]); samp=res.samples
    q=lambda j:np.percentile(np.random.default_rng(0).choice(samp[:,j],5000,p=w/w.sum()),[16,50,84])
    names=(['lrc'] if arm=='A' else [])+['lud']+(['lub'] if s.bul else [])+['zD','zi']
    post={n:q(j).tolist() for j,n in enumerate(names)}
    xmed=np.array([post[n][1] for n in names]); m=s.model(xmed,arm)
    rms=float(np.sqrt(np.mean((np.log10(m[2]**2/np.clip(m[5],1,None))[s.sel])**2))) if m is not None else None
    return dict(gal=g,arm=arm,sig_int=sig_int,logZ=float(res.logz[-1]),logZerr=float(res.logzerr[-1]),logLmax=float(res.logl.max()),post=post,rms_med=rms,bul=s.bul,npts=int(s.sel.sum()),
                D0=float(s.D0),eD=float(s.eD),i0=float(s.i0),ei=float(s.ei),n_like=int(res.ncall.sum()),sec=time.time()-t0)
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--arms',default='UBA'); ap.add_argument('--sig',default='0,0.03,0.06'); ap.add_argument('--procs',type=int,default=2)
    ap.add_argument('--slice',default=None); ap.add_argument('--pilot',action='store_true'); ap.add_argument('--nlive',type=int,default=300); ap.add_argument('--out',default='marg_ns_20261030.jsonl'); a=ap.parse_args()
    T,P,gals=load()
    if a.pilot: gals=[g for g in gals if T.loc[g,'Q']==1][:2]
    if a.slice: lo,hi=map(int,a.slice.split(':')); gals=gals[lo:hi]
    sigs=[float(x) for x in a.sig.split(',')]
    done=set()
    if os.path.exists(a.out):
        for line in open(a.out):
            try: d=json.loads(line); done.add((d['gal'],d['arm'],d['sig_int']))
            except: pass
    jobs=[(g,arm,sg,a.nlive) for sg in sigs for arm in a.arms for g in gals if (g,arm,sg) not in done]
    # (2026-10-31 PREREG_CDM2: A-at-sigs[0] filter removed; A runs at every requested sig_int)
    print('galaxies %d, jobs %d (skipping %d done), procs %d, nlive %d'%(len(gals),len(jobs),len(done),a.procs,a.nlive),flush=True); t0=time.time()
    with Pool(a.procs,initializer=init,initargs=(T,P)) as pool, open(a.out,'a') as f:
        for k,d in enumerate(pool.imap_unordered(run_one,jobs),1):
            f.write(json.dumps(d)+'\n'); f.flush()
            if k%10==0 or k==len(jobs): print('  %d/%d done, %.0f min'%(k,len(jobs),(time.time()-t0)/60),flush=True)
    print('done; %.0f min'%((time.time()-t0)/60))
