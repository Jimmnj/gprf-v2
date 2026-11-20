# pass2_globalfix_20261112.py -- EP-1 (JM 2026-11-12 22:26 "accepted"): correct pass2_table2_part.csv IN PLACE for the 3 galaxies whose bounded fit
# missed the global chi^2 minimum (RESULT_ERRPROP_20261112 sec 1). Only lrc, chi_law and law (residuals) of those 3 rows change; every other byte is kept.
# Expectation (written 22:26 before the run): afterwards median 7.42, std 0.63, robust half-width 0.50, law rms 0.0806; 3 lines differ from the backup.
import os, csv, io, json, numpy as np, pandas as pd
os.environ['NOFIT']='1'; exec(open('errprop_globalfit_20261112.py').read()); kpc=ns['kpc']
Gf=pd.read_csv('errprop_globalfit_20261112_part.csv').drop_duplicates('gal').set_index('gal')
FIX=('NGC4013','NGC0289','UGC01281')
def law_res(g,l):
    d=P0[P0.gal==g].sort_values('r'); m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values; ns['P']=d; ns['T']=pd.DataFrame([T0.loc[g]],index=[g])
    o=solve_gal(g,l); r=d.r.values; vo=d.vobs.values; lgo=np.log10(vo**2/r*1e6/kpc); return [float(x) for x in (lgo-np.log10(o.vmod2.values/r*1e6/kpc))[m]]
lines=open('pass2_table2_part.csv',newline='').read().splitlines(keepends=True); hdr=next(csv.reader([lines[0]]))
ix={k:i for i,k in enumerate(hdr)}; out=[]; n=0
for ln in lines:
    row=next(csv.reader([ln]))
    if row and row[0] in FIX:
        g=row[0]; row[ix['lrc']]=repr(float(Gf.lrc_g[g])); row[ix['chi_law']]=repr(float(Gf.chi_g[g])); row[ix['law']]=json.dumps(law_res(g,Gf.lrc_g[g]))
        b=io.StringIO(); csv.writer(b,lineterminator=ln[len(ln.rstrip('\r\n')):] or '\n').writerow(row); ln=b.getvalue(); n+=1
    out.append(ln)
open('pass2_table2_part.csv','w',newline='').write(''.join(out)); print('rows rewritten:',n)
R=pd.read_csv('pass2_table2_part.csv').drop_duplicates('gal'); hw=lambda x:(np.percentile(x,84)-np.percentile(x,16))/2
rms=np.sqrt(np.mean(np.concatenate([np.array(json.loads(x)) for x in R.law])**2))
print('N %d | median %.2f | std %.2f | robust %.2f | at bound %d | law rms %.4f'%(len(R),R.lrc.median(),R.lrc.std(ddof=1),hw(R.lrc),(R.lrc>=10.45).sum(),rms))
print(R.set_index('gal').loc[list(FIX),['lrc','chi_law']].round(3).to_string())
