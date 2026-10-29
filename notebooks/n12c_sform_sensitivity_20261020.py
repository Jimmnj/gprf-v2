# n12c_sform_sensitivity_20261020.py [UNRATIFIED] PRE-REGISTERED 2026-10-20 16:05 (cold review A1). Emergent nu at UNIVERSAL rho_c 10^7.45 for different
# saturating forms S(x) and cap values: tanh2 (baseline, cap 2), tanh1 (cap 1), tanh1.5, tanh3, tanh4 (cap 4), rat2 (x/(1+x/2), cap 2), lin (no cap; reference).
# DECLARED: the 'saturating potential-space law' claim SURVIVES if, across cap >= 1.5 forms, the binned median log nu_eff differs from tanh2 by < 0.05 dex in
# all transition bins (|log g_bar/a0| < 0.5) and < 0.10 dex over -2..+1. If the shape depends on the cap value by more than that, the claim is 'cap-tuned'.
import numpy as np, pandas as pd, sys
sys.argv=['x']; src=open('test1_capped_baryon_20261018.py').read().split('def fit(')[0]; exec(src)
from scipy.integrate import solve_ivp
from gate2c_robustness_20260925 import kpc, a0
FORMS={'tanh2':lambda x:2*np.tanh(x/2),'tanh1':lambda x:np.tanh(x),'tanh1.5':lambda x:1.5*np.tanh(x/1.5),'tanh3':lambda x:3*np.tanh(x/3),'tanh4':lambda x:4*np.tanh(x/4),'rat2':lambda x:x/(1+x/2),'lin':lambda x:x}
F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal'); l=7.45
def curve(g,S):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None)
    chi,s2b,phib0,Mb,Q=setup(g); rc=10**l; xb0=phib0/s2b; Sb0=S(xb0); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2]); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*(rc*np.exp(S(xb0-y[0]/s2b)-Sb0-y[1]/s2b)-rho_amb)]
    try: s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
    except Exception: return None
    if not s.success or s.y.shape[1]!=len(r): return None
    gb=vb2/r*1e6/kpc; gm=(vb2+G*s.y[2]/r)/r*1e6/kpc
    return pd.DataFrame(dict(gal=g,lgb=np.log10(gb),lnu=np.log10(gm/gb)))
out=[]
for name,S in FORMS.items():
    for g in F.index:
        c=curve(g,S)
        if c is not None: c['form']=name; out.append(c)
D=pd.concat(out); D.to_csv('n12c_sform_sensitivity_20261020.csv',index=False)
la0=np.log10(a0); edges=np.arange(-2.0,1.01,0.25)+la0
D['lnu_rar']=-np.log10(1-np.exp(-np.sqrt(10**D.lgb/a0)))
print('binned median log nu_eff - log nu_tanh2 (dex), universal rho_c 7.45; last col = tanh2 offset from RAR fn')
hdr='%12s'%'lgb/a0'+''.join('%9s'%f for f in FORMS if f!='tanh2')+'%10s'%'tanh2-RAR'; print(hdr)
rows=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    base=D[(D.form=='tanh2')&(D.lgb>=lo)&(D.lgb<hi)]; b=base.lnu.median(); row={'lo':lo-la0,'hi':hi-la0,'N':len(base)}
    line='%+5.2f..%+5.2f'%(lo-la0,hi-la0)
    for f in FORMS:
        if f=='tanh2': continue
        m=(D.form==f)&(D.lgb>=lo)&(D.lgb<hi); v=D.lnu[m].median()-b; row[f]=v; line+='%+9.3f'%v
    row['tanh2_dRAR']=(base.lnu-base.lnu_rar).median(); line+='%+10.3f'%row['tanh2_dRAR']; print(line); rows.append(row)
R=pd.DataFrame(rows); R.to_csv('n12c_sform_bins_20261020.csv',index=False)
t=R[(R.lo>=-0.5)&(R.hi<=0.5)]; w=R[(R.lo>=-2.0)]
for f in ('tanh1','tanh1.5','tanh3','tanh4','rat2','lin'):
    print('%-7s max|d| transition %.3f  full(-2..+1) %.3f  -> %s'%(f,t[f].abs().max(),w[f].abs().max(),'within bar' if (t[f].abs().max()<0.05 and w[f].abs().max()<0.10) else 'EXCEEDS'))
print('N galaxies solved per form:',D.groupby('form').gal.nunique().to_dict())
