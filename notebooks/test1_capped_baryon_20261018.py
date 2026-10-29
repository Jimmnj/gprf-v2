# test1_capped_baryon_20261018.py  [UNRATIFIED]  TEST 1 of the FRM-4-FREEZE law (BALLOT 2026-10-18 15:30).
# LAW: rho_nb(r) = rho_c exp( S(x_b(r)) - S(x_b0) - dpsi_nb(r)/s2 ) ;  x_b = phi_b/s2 (phi_b bump-positive baryon potential, 0 at infinity);
#      s2 = sqrt(G a0 M_b)/2 (BTFR, FIXED); S saturating, S(x)~x for x<<xsat, S -> xsat. Cap on the BARYON exponent only. rho_c free (1 param).
#      Same ODE, weights, data cuts, Yd 0.5/Yb 0.7, rho_amb 35.5, galaxy list as N14 (n14_sigma_free_20261001.csv).
# LAWS: lin (S=x; VALIDATION: must reproduce N14 lrhoc_btfr/chi_btfr), tanh2 (S=2 tanh(x/2); MAIN), sensitivity tanh1, tanh4, rat2 (x/(1+x/2)),
#       plus diagnostic tanh2_free (s2 free, 2 params).
# PRE-DECLARED (before any run): P1 median chi(tanh2) <= median chi(lin) at BTFR s2.  P2 std(lrhoc_fit - lrhoc_pred) after one universal offset
#       <= 0.30 dex AND < the lin baseline; lrhoc_pred = log rho_amb + S(x_b0)/ln10. Both reported for all and Q=1.
# Run from notebooks/: python3 test1_capped_baryon_20261018.py [workers]   (resumable: appends to test1_part.csv); then: ... report
import numpy as np, pandas as pd, sys, os
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import solve_ivp
from multiprocessing import Pool
from gate2c_robustness_20260925 import load, G, a0k
rho_amb=35.5; PART='test1_part.csv'; OUT='test1_capped_baryon_20261018.csv'
LAWS={'lin':lambda x:x,'tanh2':lambda x:2*np.tanh(x/2),'tanh1':lambda x:np.tanh(x),'tanh4':lambda x:4*np.tanh(x/4),'rat2':lambda x:x/(1+x/2)}
N14=pd.read_csv('n14_sigma_free_20261001.csv',index_col=0)
P=load(); P=P[P.gal.isin(N14.index)]
def setup(g):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); vobs=d.vobs.values
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2]); rl=r[-1]
    vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2[-1]   # bump-positive central baryon potential
    Mb=float(0.5*d.L.iloc[0]*1e9+1.33*d.MHI.iloc[0]*1e9); s2b=np.sqrt(a0k*G*Mb)/2
    w=1/(0.1*vobs)**2
    def chi(l,s2,S):
        rc=10**l; xb0=phib0/s2; Sb0=S(xb0)
        def f(x,y):
            x_=max(x,1e-4); vb=vbar(x_)
            return [vb/x_, G*y[2]/x_**2, 4*np.pi*x*x*(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-rho_amb)]
        s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
        if not s.success or s.y.shape[1]!=len(r): return 1e12
        m=vb2+G*s.y[2]/r
        return float(np.sum(w*(np.sqrt(np.clip(m,1,None))-vobs)**2))
    return chi,s2b,phib0,Mb,int(d.Q.iloc[0])
def fit(g):
    chi,s2b,phib0,Mb,Q=setup(g); out=dict(gal=g,Mb=Mb,Q=Q,s2b=s2b,phib0=phib0,xb0=phib0/s2b)
    for k,S in LAWS.items():
        res=minimize_scalar(lambda l:chi(l,s2b,S),bounds=(4,11),method='bounded',options={'xatol':0.02})
        out['l_'+k]=res.x; out['chi_'+k]=res.fun
    row=N14.loc[g]; S=LAWS['tanh2']
    res=minimize(lambda p:chi(np.clip(p[1],4,11),10**p[0],S),[np.log10(row.s2_free),row.lrhoc_free],method='Nelder-Mead',
                 options={'xatol':0.01,'fatol':0.1,'maxfev':300})
    out['s2_tf']=10**res.x[0]; out['l_tf']=float(np.clip(res.x[1],4,11)); out['chi_tf']=res.fun; out['xb0_tf']=phib0/out['s2_tf']
    with open(PART,'a') as fh: fh.write(pd.DataFrame([out]).to_csv(header=not os.path.exists(PART) or os.path.getsize(PART)==0,index=False))
    return out
def report():
    C=pd.read_csv(PART).drop_duplicates('gal').set_index('gal').sort_index(); C.to_csv(OUT)
    C=C.join(N14[['lrhoc_btfr','chi_btfr','lrhoc_free','chi_free']])
    print('N =',len(C),' | validation lin vs N14 btfr: dlrhoc median %+.3f max|%.3f|, chi ratio median %.3f'%(
        (C.l_lin-C.lrhoc_btfr).median(),(C.l_lin-C.lrhoc_btfr).abs().max(),(C.chi_lin/C.chi_btfr).median()))
    print('x_b0 at BTFR s2: median %.2f 16-84 %.2f..%.2f ; frac > 2: %.2f'%(C.xb0.median(),*np.percentile(C.xb0,[16,84]),(C.xb0>2).mean()))
    for sub,lab in ((C,'ALL'),(C[C.Q==1],'Q=1')):
        print('--- %s N=%d ---  median chi: N14 free-s2 %.1f'%(lab,len(sub),sub.chi_free.median()))
        base=None
        for k,S in list(LAWS.items()):
            pred=np.log10(rho_amb)+S(sub.xb0.values)/np.log(10); d=sub['l_'+k].values-pred; sd=d.std()
            A=np.polyfit(pred,sub['l_'+k].values,1); rr=np.corrcoef(pred,sub['l_'+k].values)[0,1] if np.std(pred)>0 else np.nan
            if k=='lin': base=sd
            print('%-6s chi med %6.1f | lrhoc med %.2f std %.2f | fit-pred offset %+.2f resid %.2f dex | slope %+.2f r %.2f | self-term/s2 med %.1f'%(
                k,sub['chi_'+k].median(),sub['l_'+k].median(),sub['l_'+k].std(),np.median(d),sd,A[0],rr,np.median(d)*np.log(10)))
        S=LAWS['tanh2']; pred=np.log10(rho_amb)+S(sub.xb0_tf.values)/np.log(10); d=sub.l_tf.values-pred
        print('tanh2 s2-free: chi med %.1f | resid %.2f dex | log(s2/s2_btfr) med %+.2f'%(sub.chi_tf.median(),d.std(),np.log10(sub.s2_tf/sub.s2b).median()))
        p1=sub.chi_tanh2.median()<=sub.chi_lin.median(); dt=sub.l_tanh2.values-(np.log10(rho_amb)+S(sub.xb0.values)/np.log(10)); p2=(dt.std()<=0.30) and (dt.std()<base)
        print('P1 (capped chi <= lin chi): %s   P2 (resid <= 0.30 and < lin): %s'%('PASS' if p1 else 'FAIL','PASS' if p2 else 'FAIL'))
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='report': report(); sys.exit()
    nw=int(sys.argv[1]) if len(sys.argv)>1 else 2
    done=set(pd.read_csv(PART).gal) if os.path.exists(PART) and os.path.getsize(PART)>0 else set()
    gals=[g for g in N14.index if g not in done and g in set(P.gal)]
    print('todo',len(gals),flush=True)
    with Pool(nw) as Pl:
        for o in Pl.imap_unordered(fit,gals): print(o['gal'],'%.2f %.1f %.1f'%(o['xb0'],o['chi_lin'],o['chi_tanh2']),flush=True)
