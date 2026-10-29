# n12b_mond_freeML_20261020.py [UNRATIFIED] PRE-REGISTERED 2026-10-20 13:30. Fair one-knob comparison (framing (az'')): MOND RAR nu with Upsilon_disk FREE
# per galaxy (Upsilon_bulge 0.7 fixed; D, i fixed) vs GPRF frozen law with rho_c free (Upsilon_disk 0.5 fixed) -- same 120 SPARC galaxies, same points,
# same weighting (10% velocity errors). Reported: per-point rms log(g_obs/g_model) and per-galaxy rms. Variants: Upsilon unconstrained in [0.1, 2.0];
# with log-normal prior 0.5 +/- 0.1 dex (Li+2018-like) added to chi; GPRF reference = fitted rho_c (from N12 csv). DECLARED: if MOND(free Y) rms < GPRF(free rho_c)
# rms 0.107, the per-galaxy claim (az') is DROPPED.
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from gate2c_robustness_20260925 import load, kpc, a0
P=load(); N12=pd.read_csv('n12_emergent_nu_20261020.csv'); gals=sorted(set(N12.gal))
nu=lambda gb: 1/(1-np.exp(-np.sqrt(gb/a0)))
rows=[]
for g in gals:
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vobs=d.vobs.values; w=1/(0.1*vobs)**2
    base=np.sign(d.vg.values)*d.vg.values**2+0.7*d.vb.values**2; vd2=d.vd.values**2
    def model(Y):
        vb2=np.clip(base+Y*vd2,1e-3,None); gb=vb2/r*1e6/kpc; return np.sqrt(gb*nu(gb)*r*kpc/1e6)   # v_model
    def chi(Y,prior):
        c=np.sum(w*(model(Y)-vobs)**2)
        return c+ (np.log10(Y/0.5)/0.1)**2 if prior else c
    out=dict(gal=g,n=len(r))
    for lab,pr in (('free',False),('prior',True)):
        res=minimize_scalar(lambda Y:chi(Y,pr),bounds=(0.1,2.0),method='bounded'); Y=res.x
        lg=np.log10((vobs**2/r)/(model(Y)**2/r)); out['Y_'+lab]=Y; out['rms_'+lab]=np.sqrt(np.mean(lg**2)); out['chi_'+lab]=np.sum(w*(model(Y)-vobs)**2)
        out['res_'+lab]=lg
    lg=np.log10((vobs**2/r)/(model(0.5)**2/r)); out['rms_fixed']=np.sqrt(np.mean(lg**2)); out['res_fixed']=lg
    rows.append(out)
R=pd.DataFrame(rows).set_index('gal')
gp=N12[N12.case=='fit'].copy(); gp['res']=gp.lgo-gp.lgm; gpr=gp.groupby('gal').res.apply(lambda s:np.sqrt((s**2).mean()))
gu=N12[N12.case=='u7.45'].copy(); gu['res']=gu.lgo-gu.lgm
def allpts(col): return np.sqrt(np.mean(np.concatenate(R[col].values)**2))
print('N galaxies %d, points %d'%(len(R),R.n.sum()))
print('per-point rms log(g_obs/g_model):  MOND Y fixed 0.5 %.3f | MOND Y free %.3f | MOND Y prior(0.1 dex) %.3f | GPRF rho_c free (Y fixed) %.3f | GPRF rho_c universal %.3f'%(
    allpts('res_fixed'),allpts('res_free'),allpts('res_prior'),np.sqrt(np.mean(gp.res**2)),np.sqrt(np.mean(gu.res**2))))
print('per-galaxy rms median:              %.3f | %.3f | %.3f | %.3f | %.3f'%(R.rms_fixed.median(),R.rms_free.median(),R.rms_prior.median(),gpr.median(),gu.groupby('gal').res.apply(lambda s:np.sqrt((s**2).mean())).median()))
j=R.join(gpr.rename('rms_gprf'))
print('galaxy-by-galaxy: GPRF(rho_c) better than MOND(Y free) in %d/%d; than MOND(Y prior) in %d/%d'%((j.rms_gprf<j.rms_free).sum(),len(j),(j.rms_gprf<j.rms_prior).sum(),len(j)))
print('MOND free Y: median %.2f, 16-84%% %.2f-%.2f; at bounds (0.1 or 2.0): %d galaxies'%(R.Y_free.median(),*np.percentile(R.Y_free,[16,84]),((R.Y_free<0.11)|(R.Y_free>1.99)).sum()))
print('chi (TEST1 weighting) median: MOND fixed %.1f | free %.1f | prior %.1f   (TEST1 GPRF tanh2 rho_c-free median 13.8; free-sigma 5.6)'%(np.median([np.sum((1/(0.1*P[P.gal==g].vobs.values)**2)*0) for g in gals])*0+R.chi_free.median()*0+np.nan if False else np.nan, R.chi_free.median(), R.chi_prior.median()))
R.drop(columns=['res_free','res_prior','res_fixed']).to_csv('n12b_mond_freeML_20261020.csv')
