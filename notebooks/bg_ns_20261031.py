# bg_ns_20261031.py -- PREREG_RHOC_BG_20261031: background-normalisation rule as a zero-freedom NS arm.
#   Arm G: rho_nb = A exp[S(x_b) - dphi_nb/sigma^2] - rho_amb with log A = LOGA fixed (default 6.75, RESULT_RHOC_RULES), i.e. rho_c = A e^{S(x_b0)} per galaxy, no per-galaxy parameter.
#   Nuisances, likelihood, sampler as marg_ns_20261030.py (arm U). Output bg_ns_20261031.jsonl. Rows carry the implied log rho_c at the posterior-median nuisances.
# Usage: python3 bg_ns_20261031.py --sig 0.03,0.06 --procs 30 [--logA 6.75] [--nlive 300] [--pilot] [--out bg_ns_20261031.jsonl]
import os, sys, json, time, argparse, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.special import erfinv
from scipy.integrate import solve_ivp
import dynesty
sys.argv_backup=sys.argv; sys.argv=['x']
src=open('marg_fit_20261030.py').read(); exec(src[:src.index("if __name__")])
sys.argv=sys.argv_backup
LN10=np.log(10); LOGA=6.75
class GalG(Gal):
    def unpack(s,x,arm):
        lud=x[0]; lub=x[1] if s.bul else LUB0; k=2 if s.bul else 1; zD=x[k]; zi=x[k+1]
        D=max(s.D0+zD*s.eD,0.1*s.D0); inc=float(np.clip(s.i0+zi*s.ei,20,90)); return None,lud,lub,D,inc
    def model(s,x,arm,return_lrc=False):
        _,lud,lub,D,inc=s.unpack(x,arm); r,vb2,vobs,ev,Mb=s.curves(lud,lub,D,inc); vb2c=np.clip(vb2,0,None)
        rl=r[-1]; rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2c])
        vbar=lambda q: np.interp(q,rr,vv) if q<=rl else vb2c[-1]*rl/q
        xs=np.linspace(1e-4,rl,4001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2c[-1]
        s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**LOGA*np.exp(Sb0)     # <-- the rule
        if return_lrc: return np.log10(rc),xb0
        fn=lambda q,y:[vbar(max(q,1e-4))/max(q,1e-4), G*y[2]/max(q,1e-4)**2, 4*np.pi*q*q*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
        so=solve_ivp(fn,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
        if not so.success or so.y.shape[1]!=len(r): return None
        return r,vb2,vobs,ev,Mb,vb2c+G*so.y[2]/r
def ptform_factory(s):
    def ptform(u):
        x=[LUD0+SIGU*np.sqrt(2)*erfinv(2*u[0]-1)]; k=1
        if s.bul: x.append(LUB0+SIGU*np.sqrt(2)*erfinv(2*u[k]-1)); k+=1
        x.append(np.sqrt(2)*erfinv(2*u[k]-1)); x.append(np.sqrt(2)*erfinv(2*u[k+1]-1)); return np.array(x)
    return ptform
def loglike_factory(s,sig_int):
    def loglike(x):
        _,lud,lub,D,inc=s.unpack(x,'G')
        if D<0.2*s.D0: return -1e100
        m=s.model(x,'G')
        if m is None: return -1e100
        r,vb2,vobs,ev,Mb,vm2=m; vm=np.sqrt(np.clip(vm2,1,None)); var=ev**2+(sig_int*vobs*LN10/2)**2
        return -0.5*np.sum((((vm-vobs)**2/var)+np.log(2*np.pi*var))[s.sel])
    return loglike
def run_one(args):
    g,arm,sig_int,nlive=args; t0=time.time(); s=GalG(g); ndim=1+(1 if s.bul else 0)+2
    rng=np.random.default_rng(abs(hash((g,arm,sig_int)))%2**32)
    sam=dynesty.NestedSampler(loglike_factory(s,sig_int),ptform_factory(s),ndim,nlive=nlive,rstate=rng,bound='multi',sample='rwalk')
    sam.run_nested(dlogz=0.1,print_progress=False)
    res=sam.results; w=np.exp(res.logwt-res.logz[-1]); samp=res.samples
    q=lambda j:np.percentile(np.random.default_rng(0).choice(samp[:,j],5000,p=w/w.sum()),[16,50,84])
    names=['lud']+(['lub'] if s.bul else [])+['zD','zi']; post={n:q(j).tolist() for j,n in enumerate(names)}
    xmed=np.array([post[n][1] for n in names]); m=s.model(xmed,'G'); lrc,xb0=s.model(xmed,'G',return_lrc=True)
    rms=float(np.sqrt(np.mean((np.log10(m[2]**2/np.clip(m[5],1,None))[s.sel])**2))) if m is not None else None
    return dict(gal=g,arm=arm,sig_int=sig_int,logA=LOGA,logZ=float(res.logz[-1]),logZerr=float(res.logzerr[-1]),logLmax=float(res.logl.max()),post=post,lrc_implied=float(lrc),xb0=float(xb0),rms_med=rms,bul=s.bul,npts=int(s.sel.sum()),
                D0=float(s.D0),eD=float(s.eD),i0=float(s.i0),ei=float(s.ei),n_like=int(res.ncall.sum()),sec=time.time()-t0)
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--sig',default='0.03,0.06'); ap.add_argument('--procs',type=int,default=2); ap.add_argument('--logA',type=float,default=6.75)
    ap.add_argument('--slice',default=None); ap.add_argument('--pilot',action='store_true'); ap.add_argument('--nlive',type=int,default=300); ap.add_argument('--out',default='bg_ns_20261031.jsonl'); a=ap.parse_args()
    LOGA=a.logA; arm='G' if abs(LOGA-6.75)<1e-6 else 'G%.2f'%LOGA
    T,P,gals=load()
    if a.pilot: gals=[g for g in gals if T.loc[g,'Q']==1][:2]
    if a.slice: lo,hi=map(int,a.slice.split(':')); gals=gals[lo:hi]
    sigs=[float(x) for x in a.sig.split(',')]; done=set()
    if os.path.exists(a.out):
        for line in open(a.out):
            try: d=json.loads(line); done.add((d['gal'],d['arm'],d['sig_int']))
            except: pass
    jobs=[(g,arm,sg,a.nlive) for sg in sigs for g in gals if (g,arm,sg) not in done]
    print('galaxies %d, jobs %d (skipping %d done), procs %d, nlive %d, logA %.2f'%(len(gals),len(jobs),len(done),a.procs,a.nlive,LOGA),flush=True); t0=time.time()
    with Pool(a.procs,initializer=init,initargs=(T,P)) as pool, open(a.out,'a') as f:
        for k,d in enumerate(pool.imap_unordered(run_one,jobs),1):
            f.write(json.dumps(d)+'\n'); f.flush()
            if k%10==0 or k==len(jobs): print('  %d/%d done, %.0f min'%(k,len(jobs),(time.time()-t0)/60),flush=True)
    print('done; %.0f min'%((time.time()-t0)/60))
