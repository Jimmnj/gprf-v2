# rhoc_sigmaD_1b_20261031.py -- PREREG_RHOC_SIGMAD Stage 1b (rule at A's posterior-median nuisances) + P5 (CDR slope under the Milgrom rule, RESULT_CDR recipe).
import sys,json,os,numpy as np,pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import spearmanr
import marg_fit_20261030 as M
from cdr_test_20261030 import toomre,fit,sb0 as sb0_nom
G,a0k,ramb,S,LUD0,LUB0=M.G,M.a0k,M.ramb,M.S,M.LUD0,M.LUB0
T,P,gals=M.load(); M.init(T,P)
SIGM=a0k/(2*np.pi*G)*1e-6
def sigc_star(sb,gam=1.0): y=sb/SIGM; return SIGM/np.sqrt(gam**-2+(2.4*np.sqrt(y))**-2)
def solve(s,lrc,lud,lub,D,inc,rmax_fac=50):
    r,vb2,vobs,ev,Mb=s.curves(lud,lub,D,inc); vb2c=np.clip(vb2,0,None)
    rl=r[-1]; rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2c])
    vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2c[-1]*rl/x
    xs=np.linspace(1e-4,rl,4001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2c[-1]
    s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
    rho=lambda x,y: max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)
    fn=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*rho(x,y)]
    ev_=lambda x,y: rho(x,y)+ramb-1e-6*rc; ev_.terminal=True
    g1=np.unique(np.concatenate([r,np.geomspace(1e-3,rl,200)])); g2=np.geomspace(rl,rmax_fac*rl,150)[1:]
    s1=solve_ivp(fn,[1e-4,rl],[0,0,0],t_eval=g1,rtol=1e-6,atol=1e-3,max_step=rl/200)
    s2_=solve_ivp(fn,[rl,rmax_fac*rl],s1.y[:,-1],t_eval=g2,rtol=1e-6,atol=1e-3,max_step=rl,events=ev_)
    rg=np.concatenate([s1.t,s2_.t]); Y=np.concatenate([s1.y,s2_.y],axis=1)
    rh=np.clip(np.array([rho(x,Y[:,i]) for i,x in enumerate(rg)]),0,None)
    col=2*(np.trapezoid(rh,rg)+rh[0]*rg[0]); Mnb=np.interp(r,rg,Y[2]); return col*1e-6, r, vb2c+G*Mnb/r, vobs, Mb
def root(s,sb,lud,lub,D,inc):
    f=lambda l:(2/np.pi)*solve(s,l,lud,lub,D,inc)[0]-sigc_star(sb)
    try:
        if f(4.0)*f(11.0)>0: return np.nan
        return brentq(f,4.0,11.0,xtol=0.01)
    except Exception: return np.nan
A={}
for l in open('marg_ns_20261030.jsonl'):
    d=json.loads(l)
    if d['arm']=='A' and d['sig_int']==0.03: A[d['gal']]=d['post']
rows=[]; cdr=[]
for g in gals:
    s=M.Gal(g); p=A[g]; lud=p['lud'][1]; lub=p['lub'][1] if s.bul else LUB0; D=max(s.D0+p['zD'][1]*s.eD,0.1*s.D0); inc=float(np.clip(s.i0+p['zi'][1]*s.ei,20,90)); lrcA=p['lrc'][1]
    dd=np.atleast_2d(np.loadtxt(os.path.join(M.SD,g+'_rotmod.dat'),comments='#')); dd=dd[np.argsort(dd[:,0])]; t=T.loc[g]
    sg0=1.33*t.MHI*1e9/(2*np.pi*(2*t.Rd*1e3)**2) if t.Rd>0 else 0.0
    sbA=10**lud*dd[0,6]+10**lub*dd[0,7]+sg0*(D/s.D0)**0      # SB distance-independent; gas term kept nominal
    sbN=0.5*dd[0,6]+0.7*dd[0,7]+sg0
    lrc_1b=root(s,sbA,lud,lub,D,inc); lrc_nom=root(s,sbN,LUD0,LUB0,s.D0,s.i0)
    rows.append(dict(gal=g,lrcA=lrcA,lrc_1b=lrc_1b,lrc_nom=lrc_nom,lud=lud,zD=p['zD'][1],zi=p['zi'][1]))
    # P5: CDR at nominal nuisances, arm M0 (rule rho_c) -- same recipe as cdr_test one()
    if s.i0<=85:
        r,vb2,vobs,ev,Mb=s.curves(LUD0,LUB0,s.D0,s.i0); Sstar=sb0_nom(g); Sobs=toomre(r,vobs**2)
        if np.isfinite(lrc_nom):
            _,r2,vm2,_,_=solve(s,lrc_nom,LUD0,LUB0,s.D0,s.i0); cdr.append(dict(gal=g,Sstar=Sstar,Sobs=Sobs,SM0=toomre(r2,vm2)))
    print(g,'A %.2f  rule@A-nuis %.2f  rule@nom %.2f'%(lrcA,lrc_1b,lrc_nom),flush=True)
R=pd.DataFrame(rows); R.to_csv('rhoc_sigmaD_1b_20261031.csv',index=False); C=pd.DataFrame(cdr); C.to_csv('rhoc_sigmaD_cdr_20261031.csv',index=False)
o=R.dropna(subset=['lrc_1b','lrcA']); d=o.lrcA-o.lrc_1b
print('\nSTAGE 1b (rule at A nuisances) N=%d: median %+.2f rms %.2f sd %.2f Spearman %.2f slope %.2f'%(len(o),d.median(),np.sqrt((d**2).mean()),d.std(),spearmanr(o.lrcA,o.lrc_1b)[0],np.polyfit(o.lrc_1b,o.lrcA,1)[0]))
o2=R.dropna(subset=['lrc_nom','lrcA']); d2=o2.lrcA-o2.lrc_nom
print('nominal (Stage 1 repeat) N=%d: median %+.2f rms %.2f Spearman %.2f'%(len(o2),d2.median(),np.sqrt((d2**2).mean()),spearmanr(o2.lrcA,o2.lrc_nom)[0]))
print('rule shift from nuisances (1b - nom): median %+.2f sd %.2f; corr(A - rule_nom, zD) %.2f, corr(A - rule_nom, lud) %.2f'%((o.lrc_1b-o.lrc_nom).median(),(o.lrc_1b-o.lrc_nom).std(),spearmanr(d2,o2.zD)[0],spearmanr(d2,o2.lud)[0]))
Cc=C[(C.SM0>0)&(C.Sobs>0)]
bo,so,rmso=fit(np.log10(Cc.Sstar),np.log10(Cc.Sobs)); bm,sm,rmsm=fit(np.log10(Cc.Sstar),np.log10(Cc.SM0)); dl=np.log10(Cc.SM0)-np.log10(Cc.Sobs)
print('P5 CDR N=%d: OBS beta %.3f+-%.3f | M0 (Milgrom rule) beta %.3f+-%.3f resid rms %.3f | Delta vs OBS rms %.3f med %+.3f | C2 (|dbeta|<=2sig): %s'%(len(Cc),bo[2],so[2],bm[2],sm[2],rmsm,np.sqrt(np.mean(dl**2)),np.median(dl),abs(bm[2]-bo[2])<=2*np.hypot(sm[2],so[2])))
F=pd.read_csv('cdr_fits_20261030.csv'); print(F[['arm','N','beta','ebeta','rms_delta','med_delta','C2']].to_string(index=False))
