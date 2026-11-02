# rhoc_sigmaD_20261031.py -- PREREG_RHOC_SIGMAD_20261031 Stage 1: rho_c from Milgrom's central-surface-density bound vs the fitted per-galaxy rho_c (arm A @0.03).
import sys,json,glob,os,numpy as np,pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import spearmanr
sys.argv_backup=sys.argv; sys.argv=['x']
src=open('marg_fit_20261030.py').read(); exec(src[:src.index("if __name__")])
sys.argv=sys.argv_backup
T,P,gals=load(); init(T,P)
SIGM=a0k/(2*np.pi*G)*1e-6   # Msun/pc^2  (a0k in (km/s)^2/kpc, G in kpc (km/s)^2/Msun -> Msun/kpc^2; *1e-6 -> pc^2)
print('Sigma_M = %.1f Msun/pc^2 (log %.2f)'%(SIGM,np.log10(SIGM)))
def sigc_star(sb0,gam): y=sb0/SIGM; return SIGM/np.sqrt(gam**-2+(2.4*np.sqrt(y))**-2)
def solve(s,lrc,rmax_fac=50):
    r=s.r0; vb2c=np.clip(s.curves(LUD0,LUB0,s.D0,s.i0)[1],0,None); Mb=s.curves(LUD0,LUB0,s.D0,s.i0)[4]
    rl=r[-1]; rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2c])
    vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2c[-1]*rl/x
    xs=np.linspace(1e-4,rl,4001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2c[-1]
    s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
    rho=lambda x,y: max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)
    fn=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*rho(x,y)]
    ev=lambda x,y: rho(x,y)+ramb-1e-6*rc; ev.terminal=True   # stop where rho_nb falls to ~1e-6 rho_c (ambient floor)
    g1=np.unique(np.concatenate([r,np.geomspace(1e-3,rl,200)])); g2=np.geomspace(rl,rmax_fac*rl,150)[1:]
    s1=solve_ivp(fn,[1e-4,rl],[0,0,0],t_eval=g1,rtol=1e-6,atol=1e-3,max_step=rl/200)
    s2_=solve_ivp(fn,[rl,rmax_fac*rl],s1.y[:,-1],t_eval=g2,rtol=1e-6,atol=1e-3,max_step=rl,events=ev)
    rg=np.concatenate([s1.t,s2_.t]); Y=np.concatenate([s1.y,s2_.y],axis=1)
    class so: pass
    so.t=rg; so.y=Y
    rh=np.array([rho(x,Y[:,i]) for i,x in enumerate(rg)]); rh=np.clip(rh,0,None)
    col=2*(np.trapezoid(rh,rg)+rh[0]*rg[0])            # Msun/kpc^2 through the centre
    Mnb=np.interp(r,rg,so.y[2]); vm2=vb2c+G*Mnb/r
    return col*1e-6, vm2, Mb, s2
OUT='rhoc_sigmaD_20261031.jsonl'; done={}
if os.path.exists(OUT):
    for l in open(OUT): d=json.loads(l); done[d['gal']]=d
rows=list(done.values())
A={}
for l in open('marg_ns_20261030.jsonl'):
    d=json.loads(l)
    if d['arm']=='A' and d['sig_int']==0.03: A[d['gal']]=d['post']['lrc'][1]
for g in gals:
    if g in done: continue
    s=Gal(g); d=P[P.gal==g].sort_values('r'); f=os.path.join(SD,g+'_rotmod.dat'); dd=np.atleast_2d(np.loadtxt(f,comments='#'))
    sbd,sbb=dd[0,6],dd[0,7]; t=T.loc[g]; sg0=1.33*t.MHI*1e9/(2*np.pi*(2*t.Rd*1e3)**2) if t.Rd>0 else 0.0
    sb0=0.5*sbd+0.7*sbb+sg0; sb0_nogas=0.5*sbd+0.7*sbb
    out=dict(gal=g,lMb=np.log10(s.curves(LUD0,LUB0,s.D0,s.i0)[4]),sb0=sb0,sb0_nogas=sb0_nogas,lrcA=A.get(g,np.nan))
    for tag,gam,sb in [('g1',1.0,sb0),('g2pi',2/np.pi,sb0),('g1_nogas',1.0,sb0_nogas)]:
        target=sigc_star(sb,gam)
        fcol=lambda lrc: (2/np.pi)*solve(s,lrc)[0]-target
        try:
            lo,hi=4.0,11.0; flo,fhi=fcol(lo),fcol(hi)
            if flo*fhi>0: lrc=np.nan
            else: lrc=brentq(fcol,lo,hi,xtol=0.01)
        except Exception as e: lrc=np.nan
        out['lrc_'+tag]=lrc; out['sigc_'+tag]=target
    # RAR per-point rms (dex in g) at rule, universal, A
    sel=s.sel; r=s.r0; vobs=s.vobs0
    for tag,lrc in [('rule',out['lrc_g1']),('univ',7.45),('A',out['lrcA'])]:
        if np.isfinite(lrc):
            vm2=solve(s,lrc)[1]; out['rms_'+tag]=float(np.sqrt(np.mean((np.log10(vobs[sel]**2/np.clip(vm2[sel],1,None)))**2)))
        else: out['rms_'+tag]=np.nan
    rows.append(out); open(OUT,'a').write(json.dumps({k:(None if (isinstance(v,float) and np.isnan(v)) else v) for k,v in out.items()})+'\n'); print(g,'lMb %.2f sb0 %6.1f sigc* %6.1f  lrc rule %.2f  A %.2f  rms rule %.3f univ %.3f A %.3f'%(out['lMb'],sb0,out['sigc_g1'],out['lrc_g1'],out['lrcA'],out['rms_rule'],out['rms_univ'],out['rms_A']),flush=True)
R=pd.DataFrame(rows); R.to_csv('rhoc_sigmaD_20261031.csv',index=False)
ok=R.dropna(subset=['lrc_g1','lrcA'])
print('\nN with solution:',len(ok),'of',len(R))
for tag in ['g1','g2pi','g1_nogas']:
    o=R.dropna(subset=['lrc_'+tag,'lrcA']); dlt=o['lrcA']-o['lrc_'+tag]
    print('P1 %-9s median %+.2f  rms %.2f  sd %.2f  | P2 Spearman %.2f slope %.2f'%(tag,dlt.median(),np.sqrt((dlt**2).mean()),dlt.std(),spearmanr(o['lrcA'],o['lrc_'+tag])[0],np.polyfit(o['lrc_'+tag],o['lrcA'],1)[0]))
print('P3 A vs lMb: Spearman %.2f slope %.2f ; sd of A about its mean %.2f'%(spearmanr(ok.lrcA,ok.lMb)[0],np.polyfit(ok.lMb,ok.lrcA,1)[0],ok.lrcA.std()))
print('rule vs lMb slope %.2f'%np.polyfit(ok.lMb,ok.lrc_g1,1)[0])
print('P4 pooled rms: rule %.4f  univ %.4f  A %.4f  (median per galaxy: %.3f %.3f %.3f)'%(np.sqrt((ok.rms_rule**2).mean()),np.sqrt((ok.rms_univ**2).mean()),np.sqrt((ok.rms_A**2).mean()),ok.rms_rule.median(),ok.rms_univ.median(),ok.rms_A.median()))
