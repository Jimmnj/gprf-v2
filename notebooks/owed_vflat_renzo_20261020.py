# owed_vflat_renzo_20261020.py [UNRATIFIED] 2026-10-20 17:15 -- two owed checks (cold review A9; debate table rows 3, 4) from n12_emergent_nu csv.
# (1) BTFR scatter: model V_flat (mean of outer-third model V) vs M_b, residual std in log V at fixed M_b, for fitted rho_c / universal / MOND RAR fn / data.
#     DECLARED: pass if fitted-rho_c model BTFR scatter <= data scatter + 0.02 dex.
# (2) Renzo per-feature: features = log V - quadratic-in-log-r smooth per galaxy (>= 8 points). Pooled Pearson r and slope of feature_model vs feature_obs for
#     fitted / universal / MOND RAR fn / baryons-only (V_bar, Y fixed). DECLARED: 'potential explanation' holds if GPRF-fitted r is within 0.05 of MOND's.
import numpy as np, pandas as pd
from gate2c_robustness_20260925 import load_table, kpc, a0
D=pd.read_csv('n12_emergent_nu_20261020.csv'); T=load_table()
D['lV_obs']=0.5*(D.lgo+np.log10(D.r*kpc/1e6)); D['lV_mod']=0.5*(D.lgm+np.log10(D.r*kpc/1e6))
D['lV_bar']=0.5*(D.lgb+np.log10(D.r*kpc/1e6))
nurar=lambda lgb: -np.log10(1-np.exp(-np.sqrt(10**lgb/a0)))
D['lV_mond']=0.5*(D.lgb+nurar(D.lgb)+np.log10(D.r*kpc/1e6))
# (1) BTFR
rows=[]
for (g,case),d in D.groupby(['gal','case']):
    d=d.sort_values('r'); n=len(d); o=d.iloc[-max(3,n//3):]
    Mb=0.5*T.loc[g,'L']*1e9+1.33*T.loc[g,'MHI']*1e9
    rows.append(dict(gal=g,case=case,lMb=np.log10(Mb),vf_obs=o.lV_obs.mean(),vf_mod=o.lV_mod.mean(),vf_mond=o.lV_mond.mean()))
B=pd.DataFrame(rows)
def btfr(x,y,lab):
    A=np.polyfit(x,y,1); r=y-np.polyval(A,x); print('  %-22s slope %.3f  scatter %.3f dex (N=%d)'%(lab,A[0],r.std(),len(x))); return r.std()
print('(1) BTFR in model V_flat (log V vs log M_b, residual std):')
b=B[B.case=='fit']; s_data=btfr(b.lMb,b.vf_obs,'data'); s_fit=btfr(b.lMb,b.vf_mod,'GPRF fitted rho_c'); btfr(b.lMb,b.vf_mond,'MOND RAR fn')
u=B[B.case=='u7.45']; btfr(u.lMb,u.vf_mod,'GPRF universal rho_c')
print('  -> fitted-rho_c scatter %.3f vs data %.3f: %s (bar data+0.02)'%(s_fit,s_data,'PASS' if s_fit<=s_data+0.02 else 'FAIL'))
# (2) Renzo features
def feat(d,col):
    x=np.log10(d.r.values); y=d[col].values
    if len(x)<8: return None
    return y-np.polyval(np.polyfit(x,y,2),x)
res={}
for case in ('fit','u7.45'):
    F_obs=[];F_mod=[];F_bar=[];F_mond=[]
    for g,d in D[D.case==case].groupby('gal'):
        d=d.sort_values('r'); fo=feat(d,'lV_obs')
        if fo is None: continue
        F_obs.append(fo); F_mod.append(feat(d,'lV_mod')); F_bar.append(feat(d,'lV_bar')); F_mond.append(feat(d,'lV_mond'))
    fo=np.concatenate(F_obs); fm=np.concatenate(F_mod); fb=np.concatenate(F_bar); fmo=np.concatenate(F_mond)
    res[case]=(fo,fm,fb,fmo)
print('(2) Renzo per-feature: pooled Pearson r and slope of feature_model vs feature_obs (features = log V minus per-galaxy quadratic in log r):')
fo,fm,fb,fmo=res['fit']; n=len(fo)
for lab,f in (('GPRF fitted rho_c',fm),('GPRF universal rho_c',res['u7.45'][1]),('MOND RAR fn',fmo),('baryons only (V_bar)',fb)):
    r=np.corrcoef(fo,f)[0,1]; sl=np.polyfit(f,fo,1)[0]; print('  %-22s r = %.3f  slope(obs on model) = %.2f  feature rms model %.3f (obs %.3f)  N=%d'%(lab,r,sl,f.std(),fo.std(),n))
print('  feature rms of the nb-only part of the model (log V_mod - log V_bar, detrended) vs MOND phantom part:')
for lab,(a_,b_) in (('GPRF fitted',(fm,fb)),('MOND',(fmo,fb))):
    print('     %-12s rms %.3f'%(lab,(a_-b_).std()))
