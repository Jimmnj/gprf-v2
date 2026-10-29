# n8d_chae_method_20261021.py [UNRATIFIED] PRE-REGISTERED 2026-10-21 02:05 (Chae-persona review, required change 1). Chae+2021's statistic applied (i) to the observed SPARC
# curves and (ii) to the FROZEN LAW's own curves (fitted rho_c; universal 10^7.45) at the same radii with the same errors.
# Model: g = nu_eN(y) g_N, nu_eN (Chae+2021 eq. nueN): 1/2 + sqrt(D^2/4 + D/y) - etilde*C/y, D = 1 + |eN|/y, C = sqrt(1 + |eN|/4), etilde = eN/sqrt|eN| (so eN = etilde|etilde|).
# Nuisance: Y_disk (lognormal 0.5 +/- 0.1 dex), Y_bulge (0.7 +/- 0.1 dex), distance factor (SPARC e_D), inclination shift (SPARC e_inc), gas 1.33 fixed. chi2 with SPARC velocity errors.
# Selection: Chae's x_{0,3} < -10.6 from his Table 2 (data/chae2021). Statistics: sample median etilde, fraction > 0, sign-test p; correlation with log e_N,env.
# DECLARED: if the law's own curves return median etilde ~ 0 (|median| < 0.01, fraction positive ~0.5) while the data return etilde > 0 as in Chae+2021, his detection is evidence
# AGAINST the law. If the law's curves also return etilde > 0 of comparable size, his statistic does not discriminate the law from the EFE.
import numpy as np, pandas as pd
from scipy.optimize import minimize
from gate2c_robustness_20260925 import load, kpc, a0
rows=[]
for line in open('../sparc_galaxy_data/SPARC_Lelli2016c.mrt'):
    p=line.split()
    if len(p)>=18 and p[1].isdigit():
        try: rows.append(dict(gal=p[0],D=float(p[2]),eD=float(p[3]),inc=float(p[5]),einc=float(p[6])))
        except ValueError: pass
E=pd.DataFrame(rows).set_index('gal'); P=load(); N12=pd.read_csv('n12_emergent_nu_20261020.csv')
fit=pd.read_csv('../data/chae2021/chae2021_table2_fits.csv').set_index('galaxy'); env=pd.read_csv('../data/chae2021/chae2021_table3_eNenv.csv').set_index('galaxy')
def nu_eN(y,et):
    eN=et*abs(et); D=1+abs(eN)/y; C=np.sqrt(1+abs(eN)/4); return 0.5+np.sqrt(D**2/4+D/y)-et*C/y
def fit_gal(g,vtarget):
    d=P[P.gal==g].sort_values('r'); r0=d.r.values; ev=np.maximum(d.ev.values,2.0); i0=np.radians(E.loc[g,'inc']); D0,eD,ei=E.loc[g,['D','eD','einc']]
    vg2=np.sign(d.vg.values)*d.vg.values**2; vd2=d.vd.values**2; vb2=d.vb.values**2; hasb=vb2.max()>0
    def model(th):
        Yd,Yb,fD,di,et=th; r=r0*fD; vbar2=np.clip((vg2+Yd*vd2+Yb*vb2)*fD,1e-3,None); gN=vbar2/r*1e6/kpc; y=gN/a0
        vm=np.sqrt(np.maximum(nu_eN(y,et)*gN,0)*r*kpc/1e6); vo=vtarget*np.sin(i0)/np.sin(i0+np.radians(di)); return vm,vo
    def chi(th):
        Yd,Yb,fD,di,et=th
        if Yd<=0.05 or Yb<=0.05 or fD<=0.3 or not(5<np.degrees(i0)+di<89) or abs(et)>0.6: return 1e12
        vm,vo=model(th); c=np.sum(((vm-vo)/ev)**2)+(np.log10(Yd/0.5)/0.1)**2+(np.log10(Yb/0.7)/0.1)**2+((fD-1)*D0/eD)**2+(di/ei)**2
        return c
    best=None
    for et0 in (-0.05,0.0,0.05,0.15):
        res=minimize(chi,[0.5,0.7,1.0,0.0,et0],method='Nelder-Mead',options={'xatol':1e-4,'fatol':1e-3,'maxiter':2000})
        if best is None or res.fun<best.fun: best=res
    return best.x[4], best.fun
out=[]
gals=sorted(set(N12.gal))
for g in gals:
    if g not in fit.index or g not in E.index: continue
    d=P[P.gal==g].sort_values('r'); n=len(d)
    rec=dict(gal=g,x03=fit.loc[g,'x03'],et_chae=fit.loc[g,'etilde'],leN=env.loc[g,'log_eN_env_maxclust'] if g in env.index else np.nan)
    rec['et_obs'],_=fit_gal(g,d.vobs.values)
    for case,lab in (('fit','et_law_fit'),('u7.45','et_law_uni')):
        m=N12[(N12.case==case)&(N12.gal==g)].sort_values('r'); vmod=np.sqrt(10**m.lgm.values*m.r.values*kpc/1e6)
        rec[lab],_=fit_gal(g,vmod) if len(vmod)==n else (np.nan,np.nan)
    out.append(rec)
R=pd.DataFrame(out).set_index('gal'); R.to_csv('n8d_chae_method_20261021.csv')
from scipy.stats import binomtest, spearmanr
sel=R[R.x03<-10.6]
print('N fitted %d; Chae selection x03 < -10.6: %d'%(len(R),len(sel)))
print('%-28s %8s %8s %8s %10s'%('curves','median','frac>0','N','sign-test p'))
for col,lab in (('et_chae','Chae+2021 own fits (table)'),('et_obs','observed, our refit'),('et_law_fit','LAW curves, fitted rho_c'),('et_law_uni','LAW curves, rho_c 10^7.45')):
    x=sel[col].dropna(); k=int((x>0).sum()); print('%-28s %+8.3f %8.2f %8d %10.1e'%(lab,x.median(),k/len(x),len(x),binomtest(k,len(x)).pvalue))
m=sel.dropna(subset=['leN','et_obs']); print('corr(etilde_obs, log eN_env): spearman %+.3f (N=%d); corr(Chae table etilde, log eN): %+.3f'%(spearmanr(m.et_obs,m.leN)[0],len(m),spearmanr(m.et_chae,m.leN)[0]))
print('our refit vs Chae table etilde: spearman %+.3f, median diff %+.3f'%(spearmanr(sel.et_obs,sel.et_chae)[0],(sel.et_obs-sel.et_chae).median()))
