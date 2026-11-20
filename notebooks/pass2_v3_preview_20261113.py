# pass2_v3_preview_20261113.py -- BD-2 (b) preview: fixed-arm table stats if NGC 0289 and NGC 5907 take their V3 first-minimum values.
# READ-ONLY on pass2_table2_part.csv (no in-place change without a JM ruling). law_res copied verbatim from pass2_globalfix_20261112.py. [UNRATIFIED]
import os, json, numpy as np, pandas as pd
os.environ['NOFIT']='1'; exec(open('errprop_globalfit_20261112.py').read()); kpc=ns['kpc']
def law_res(g,l):
    d=P0[P0.gal==g].sort_values('r'); m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values; ns['P']=d; ns['T']=pd.DataFrame([T0.loc[g]],index=[g])
    o=solve_gal(g,l); r=d.r.values; vo=d.vobs.values; lgo=np.log10(vo**2/r*1e6/kpc); return [float(x) for x in (lgo-np.log10(o.vmod2.values/r*1e6/kpc))[m]]
R=pd.read_csv('pass2_table2_part.csv').drop_duplicates('gal').set_index('gal'); F=pd.read_csv('bound_check_first_20261112.csv').set_index('gal')
hw=lambda x:(np.percentile(x,84)-np.percentile(x,16))/2
def stats(R,lab):
    rms=np.sqrt(np.mean(np.concatenate([np.array(json.loads(x) if isinstance(x,str) else x) for x in R.law])**2))
    print('%-10s N %d | median %.3f | std %.3f | robust %.3f | >=10.45 %d | law rms %.4f'%(lab,len(R),R.lrc.median(),R.lrc.std(ddof=1),hw(R.lrc),(R.lrc>=10.45).sum(),rms))
stats(R,'current')
V=R.copy().astype({'law':object})
for g in ('NGC0289','NGC5907'):
    lr=F.loc[g,'lrc_first']; res=law_res(g,lr); V.at[g,'lrc']=lr; V.at[g,'chi_law']=F.loc[g,'chi_first']; V.at[g,'law']=json.dumps(res)
    print('  %s lrc %.3f -> %.3f  chi %.2f -> %.2f  galaxy rms %.4f -> %.4f'%(g,R.loc[g,'lrc'],lr,R.loc[g,'chi_law'],F.loc[g,'chi_first'],np.sqrt(np.mean(np.array(json.loads(R.loc[g,'law']))**2)),np.sqrt(np.mean(np.array(res)**2))))
stats(V,'V3')
