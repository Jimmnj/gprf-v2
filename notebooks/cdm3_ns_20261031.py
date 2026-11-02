# cdm3_ns_20261031.py -- PREREG_CDM3_20261031: cored-halo arms under the PREREG_NS evidence machinery.
#   D: DC14 feedback-cored halo (Di Cintio+2014 alpha,beta,gamma from X=log10(M*/M200), clamped to [-4.1,-1.3]);
#      priors log M200 ~ U[8.5,14]; log c200(NFW) = Dutton-Maccio(M200) + delta, delta ~ N(0,0.11); c_DC14 = c_NFW (1 + 1e-5 exp[3.4 (X+4.5)]);
#      c200 = r200 / r_-2, r_-2 = r_s ((2-gamma)/(beta-2))^(1/alpha). M* = Ud L (fitted Ud; the MRT has no bulge luminosity).
#   K: Burkert halo, flat priors log rho0 ~ U[5,10] Msun/kpc^3, log r_c ~ U[-1,2] kpc. Two free parameters, no cosmological prior.
# Same nuisances, likelihood, sampler settings as cdm_ns_20261031.py. Output cdm3_ns_20261031.jsonl.
# Usage: python3 cdm3_ns_20261031.py --arms DK --sig 0.03,0.06 --procs 30 [--nlive 300] [--pilot] [--out cdm3_ns_20261031.jsonl]
import os, sys, json, time, argparse, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.special import erfinv
import dynesty
sys.argv_backup=sys.argv; sys.argv=['x']
exec(open('marg_fit_20261030.py').read().split("if __name__")[0])
sys.argv=sys.argv_backup
LN10=np.log(10); H0=70.0/1000.0; RHOC=3*H0**2/(8*np.pi*G)
def cM(lM): return 0.905-0.101*(lM-np.log10(1e12/0.7))
def dc14_abg(X):
    X=float(np.clip(X,-4.1,-1.3)); t=10**(X+2.33); u=10**(X+2.56)
    al=2.94-np.log10(t**-1.08+t**2.29); be=4.23+1.34*X+0.26*X*X; ga=-0.06+np.log10(u**-0.68+u); return al,be,ga
XG=np.logspace(-3,3.5,400)
def vdc14_2(r,lM,lcN,lMstar):
    M=10**lM; X=lMstar-lM; al,be,ga=dc14_abg(X); c=10**lcN*(1+1e-5*np.exp(3.4*(np.clip(X,-4.1,-1.3)+4.5)))
    r200=(3*M/(4*np.pi*200*RHOC))**(1/3); rm2=r200/c; rs=rm2/((2-ga)/(be-2))**(1/al)
    rho=XG**-ga*(1+XG**al)**(-(be-ga)/al); m=np.concatenate([[0],np.cumsum(0.5*(rho[1:]*XG[1:]**2+rho[:-1]*XG[:-1]**2)*np.diff(XG))])
    m200=np.interp(r200/rs,XG,m); mr=np.interp(r/rs,XG,m)*M/m200
    return G*mr/r
def vburk2(r,lr0,lrc):
    rho0=10**lr0; rc=10**lrc; x=r/rc
    return G*2*np.pi*rho0*rc**3*(np.log(1+x)-np.arctan(x)+0.5*np.log(1+x*x))/r
class GalC(Gal):
    def unpack(s,x,arm):
        p1,p2=x[0],x[1]; lud=x[2]; lub=x[3] if s.bul else LUB0; k=4 if s.bul else 3; zD=x[k]; zi=x[k+1]
        D=max(s.D0+zD*s.eD,0.1*s.D0); inc=float(np.clip(s.i0+zi*s.ei,20,90)); return (p1,p2),lud,lub,D,inc
    def model(s,x,arm):
        (p1,p2),lud,lub,D,inc=s.unpack(x,arm); r,vb2,vobs,ev,Mb=s.curves(lud,lub,D,inc)
        if arm=='D': lMstar=np.log10(10**lud*s.t.L*1e9*(D/s.D0)**2); vh2=vdc14_2(r,p1,cM(p1)+p2,lMstar)
        else: vh2=vburk2(r,p1,p2)
        return r,vb2,vobs,ev,Mb,np.clip(vb2,0,None)+vh2
def ptform_factory(s,arm):
    def ptform(u):
        if arm=='D': x=[8.5+5.5*u[0], 0.11*np.sqrt(2)*erfinv(2*u[1]-1)]
        else: x=[5+5*u[0], -1+3*u[1]]
        k=2; x.append(LUD0+SIGU*np.sqrt(2)*erfinv(2*u[k]-1)); k+=1
        if s.bul: x.append(LUB0+SIGU*np.sqrt(2)*erfinv(2*u[k]-1)); k+=1
        x.append(np.sqrt(2)*erfinv(2*u[k]-1)); x.append(np.sqrt(2)*erfinv(2*u[k+1]-1)); return np.array(x)
    return ptform
def loglike_factory(s,arm,sig_int):
    def loglike(x):
        _,lud,lub,D,inc=s.unpack(x,arm)
        if D<0.2*s.D0: return -1e100
        r,vb2,vobs,ev,Mb,vm2=s.model(x,arm); vm=np.sqrt(np.clip(vm2,1,None))
        var=ev**2+(sig_int*vobs*LN10/2)**2
        return -0.5*np.sum((((vm-vobs)**2/var)+np.log(2*np.pi*var))[s.sel])
    return loglike
def run_one(args):
    g,arm,sig_int,nlive=args; t0=time.time(); s=GalC(g); ndim=2+1+(1 if s.bul else 0)+2
    rng=np.random.default_rng(abs(hash((g,arm,sig_int)))%2**32)
    sam=dynesty.NestedSampler(loglike_factory(s,arm,sig_int),ptform_factory(s,arm),ndim,nlive=nlive,rstate=rng,bound='multi',sample='rwalk')
    sam.run_nested(dlogz=0.1,print_progress=False)
    res=sam.results; w=np.exp(res.logwt-res.logz[-1]); samp=res.samples
    q=lambda j:np.percentile(np.random.default_rng(0).choice(samp[:,j],5000,p=w/w.sum()),[16,50,84])
    names=(['lM200','dlc'] if arm=='D' else ['lrho0','lrc'])+['lud']+(['lub'] if s.bul else [])+['zD','zi']
    post={n:q(j).tolist() for j,n in enumerate(names)}
    xmed=np.array([post[n][1] for n in names]); m=s.model(xmed,arm)
    rms=float(np.sqrt(np.mean((np.log10(m[2]**2/np.clip(m[5],1,None))[s.sel])**2)))
    extra={}
    if arm=='D':
        _,lud,_,D,_=s.unpack(xmed,arm); X=np.log10(10**lud*s.t.L*1e9*(D/s.D0)**2)-post['lM200'][1]; al,be,ga=dc14_abg(X); extra=dict(X=float(X),abg=[al,be,ga],lcN_med=float(cM(post['lM200'][1])+post['dlc'][1]))
    return dict(gal=g,arm=arm,sig_int=sig_int,logZ=float(res.logz[-1]),logZerr=float(res.logzerr[-1]),logLmax=float(res.logl.max()),post=post,rms_med=rms,bul=s.bul,npts=int(s.sel.sum()),
                D0=float(s.D0),eD=float(s.eD),i0=float(s.i0),ei=float(s.ei),n_like=int(res.ncall.sum()),sec=time.time()-t0,**extra)
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--arms',default='DK'); ap.add_argument('--sig',default='0.03,0.06'); ap.add_argument('--procs',type=int,default=2)
    ap.add_argument('--slice',default=None); ap.add_argument('--pilot',action='store_true'); ap.add_argument('--nlive',type=int,default=300); ap.add_argument('--out',default='cdm3_ns_20261031.jsonl'); a=ap.parse_args()
    T,P,gals=load()
    if a.pilot: gals=[g for g in gals if T.loc[g,'Q']==1][:2]
    if a.slice: lo,hi=map(int,a.slice.split(':')); gals=gals[lo:hi]
    sigs=[float(x) for x in a.sig.split(',')]; done=set()
    if os.path.exists(a.out):
        for line in open(a.out):
            try: d=json.loads(line); done.add((d['gal'],d['arm'],d['sig_int']))
            except: pass
    jobs=[(g,arm,sg,a.nlive) for sg in sigs for arm in a.arms for g in gals if (g,arm,sg) not in done]
    print('galaxies %d, jobs %d (skipping %d done), procs %d, nlive %d'%(len(gals),len(jobs),len(done),a.procs,a.nlive),flush=True); t0=time.time()
    with Pool(a.procs,initializer=init,initargs=(T,P)) as pool, open(a.out,'a') as f:
        for k,d in enumerate(pool.imap_unordered(run_one,jobs),1):
            f.write(json.dumps(d)+'\n'); f.flush()
            if k%10==0 or k==len(jobs): print('  %d/%d done, %.0f min'%(k,len(jobs),(time.time()-t0)/60),flush=True)
    print('done; %.0f min'%((time.time()-t0)/60))
