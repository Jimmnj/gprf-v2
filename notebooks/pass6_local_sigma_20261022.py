# pass6_local_sigma_20261022.py [UNRATIFIED] 2026-10-22 22:05. Owed test from RESULT_PASS5 reading (iii): is the gas-fraction pattern a LOCAL sigma?
# Local-sigma proxy v2 (22:12): sigma^2 = sqrt(G a0 M_b(<R_last))/2, M_b(<R_last) = stars from SPARC SBdisk/SBbul integrated (M/L 0.5/0.7)
#   + gas M_g [1-(1+x)e^-x], x = R_last/h_g, h_g = R_HI/2 (pass 4; fallback 1.5 R_d). v1 used R V_bar^2/G: over-states disc enclosed mass (median 1.17 M_b) -> DISCARDED.
# (Motivation: Famaey, Khoury & Penco 2018, local v^4 ~ a0 G M_b(<r).) Compare with the BTFR (total M_b) sigma: (a) one global rho_c chosen by
# per-point rms of log g; last-point ratio vs Sigma_h, f_gas; (b) per-galaxy rho_c fit chi2 (catalogue errors, 3% floor) vs pass5 chi1.
import numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.argv=['x','none']; exec(open('pass5_sigma_free_20261022.py').read().split("if sys.argv[1]=='time'")[0])   # reuse data + constants
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr

RHI={}
for line in open(MRT):
    p_=line.split()
    if len(p_)>=18 and p_[1].isdigit():
        try: RHI[p_[0]]=float(p_[14])
        except ValueError: pass
def Mb_enclosed(g,rl):
    d=np.atleast_2d(np.loadtxt(os.path.join(SD,g+'_rotmod.dat'),comments='#')); d=d[d[:,0]<=rl+1e-9]
    R=np.concatenate([[0],d[:,0]]); sb=0.5*d[:,6]+0.7*d[:,7]; sb=np.concatenate([[sb[0]],sb])
    Mst=np.trapezoid(2*np.pi*R*sb*1e6,R)
    t=T.loc[g]; Mg=1.33*t.MHI*1e9; hg=RHI.get(g,0)/2 if RHI.get(g,0)>0 else 1.5*t.Rd
    x=rl/hg; return Mst+Mg*(1-(1+x)*np.exp(-x)), Mst
def model2(g,lrc,mode):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2])
    vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2[-1]
    t=T.loc[g]; Mtot=0.5*t.L*1e9+1.33*t.MHI*1e9; Mloc=Mb_enclosed(g,rl)[0]
    Ms=Mtot if mode=='btfr' else Mloc
    s2=np.sqrt(a0k*G*Ms)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-5,atol=1e-3,max_step=rl/100)
    if not s.success or s.y.shape[1]!=len(r): return None,d,Mtot,Mloc
    return np.sqrt(vb2+G*s.y[2]/r),d,Mtot,Mloc
info={}
for mode in ('btfr','local'):
    best=None
    for lrc in np.arange(7.0,8.01,0.1):
        res=[];last=[]
        for g in gals:
            v,d,Mt,Ml=model2(g,lrc,mode)
            if v is None: continue
            res.append(2*np.log10(v/d.vobs.values)); last.append((g,v[-1]/d.vobs.values[-1],Mt,Ml))
        rms=np.sqrt(np.mean(np.concatenate(res)**2))
        if best is None or rms<best[0]: best=(rms,lrc,last)
    rms,lrc,last=best; L=pd.DataFrame(last,columns=['gal','ratio','Mt','Ml']).set_index('gal')
    t=T.loc[L.index]; L['Sh']=L.Mt/(2*np.pi*(1.678*t.Rd*1e3)**2); L['fgas']=1.33*t.MHI*1e9/L.Mt; L=L[np.isfinite(L.Sh)]
    y=np.log10(L.ratio); q=pd.qcut(np.log10(L.Sh),4,labels=False); qf=pd.qcut(L.fgas,4,labels=False)
    print('%-5s global log rho_c %.1f  per-point rms(log g) %.3f | last-pt median %.3f | Spearman vs logSh %+.2f (p %.0e), vs fgas %+.2f (p %.0e)'%(mode,lrc,rms,L.ratio.median(),*spearmanr(y,np.log10(L.Sh)),*spearmanr(y,L.fgas)))
    print('      by Sigma_h quartile:',' '.join('%.3f'%L.ratio[q==k].median() for k in range(4)),' | by f_gas quartile (low->high):',' '.join('%.3f'%L.ratio[qf==k].median() for k in range(4)))
    print('      stars(<R_last)/(0.5 L) median %.2f; '%np.median([Mb_enclosed(g,P[P.gal==g].r.max())[1]/(0.5*T.loc[g].L*1e9) for g in L.index[:40]])+'M_b(<R_last)/M_b median %.2f (16-84 %.2f..%.2f)'%((L.Ml/L.Mt).median(),*(L.Ml/L.Mt).quantile([.16,.84])))
    info[mode]=L
# (b) per-galaxy rho_c fits, chi2
F5=pd.read_csv('pass5_sigma_free_fits.csv').set_index('gal')
ch=[]
for g in gals:
    def c(lr):
        v,d,_,_=model2(g,lr,'local')
        if v is None: return 1e9
        e=np.maximum(d.ev.values,0.03*d.vobs.values); return float(np.sum(((v-d.vobs.values)/e)**2))
    r=minimize_scalar(c,bounds=(5.5,9.5),method='bounded',options=dict(xatol=0.01)); ch.append((g,r.fun))
C=pd.DataFrame(ch,columns=['gal','chi_loc']).set_index('gal').join(F5[['chi1','chi2']])
print('(b) per-galaxy rho_c fit: median chi2 BTFR-sigma %.1f | local-sigma %.1f | sigma-free (2 par) %.1f; local better than BTFR in %d/%d'%(C.chi1.median(),C.chi_loc.median(),C.chi2.median(),(C.chi_loc<C.chi1).sum(),len(C)))
C.to_csv('pass6_local_sigma_chi2.csv')
