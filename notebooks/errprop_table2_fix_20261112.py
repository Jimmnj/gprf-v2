# errprop_table2_fix_20261112.py -- Table 2 "law, rho_c free" per-point rms with the 3 corrected global fits (NGC 4013, NGC 0289, UGC 1281). [UNRATIFIED]
# Expectation (written 22:11 before the run): rms falls from 0.081 by <= 0.002 (three galaxies of 147; NGC 4013 has 36 points).
import os, json, numpy as np, pandas as pd
os.environ['NOFIT']='1'; exec(open('errprop_globalfit_20261112.py').read())
R=pd.read_csv('pass2_table2_part.csv').drop_duplicates('gal').set_index('gal'); Gf=pd.read_csv('errprop_globalfit_20261112_part.csv').drop_duplicates('gal').set_index('gal').lrc_g
kpc=ns['kpc']
def law_res(g,l):
    d=P0[P0.gal==g].sort_values('r'); m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values; ns['P']=d; ns['T']=pd.DataFrame([T0.loc[g]],index=[g])
    o=solve_gal(g,l); r=d.r.values; vo=d.vobs.values; lgo=np.log10(vo**2/r*1e6/kpc); return list((lgo-np.log10(o.vmod2.values/r*1e6/kpc))[m])
rms=lambda col: np.sqrt(np.mean(np.concatenate([np.array(json.loads(x)) for x in col])**2))
old=rms(R.law); fixed=R.law.copy()
for g in ('NGC4013','NGC0289','UGC01281'):
    a=np.array(json.loads(R.law[g])); b=np.array(law_res(g,Gf[g])); print('%-9s per-galaxy rms %.3f -> %.3f (n %d)'%(g,np.sqrt(np.mean(a**2)),np.sqrt(np.mean(b**2)),len(b))); fixed[g]=json.dumps(list(b))
print('Table 2 law, rho_c free (147 galaxies, Lelli-style points): %.4f -> %.4f ; RAR M/L free %.4f ; RAR fixed %.4f'%(old,rms(fixed),rms(R.rfree),rms(R.rfix)))
