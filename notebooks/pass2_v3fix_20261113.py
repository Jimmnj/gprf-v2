# pass2_v3fix_20261113.py -- BD-6 (JM 2026-11-13 22:04 "update the table in place"): apply BD-2 (b) to pass2_table2_part.csv IN PLACE for the
# two galaxies whose fixed-arm fit sat on a second (JLO) minimum. Only lrc, chi_law and law of those 2 rows change; every other byte is kept.
# Procedure and law_res copied verbatim from pass2_globalfix_20261112.py (EP-1). Values from bound_check_first_20261112.csv (V3).
# Expectation (written 22:05 before the run; = pass2_v3_preview_20261113 report): median 7.38, std 0.53, robust 0.50, at bound 0, law rms 0.0810; 2 lines differ from the backup.
import os, csv, io, json, numpy as np, pandas as pd
os.environ['NOFIT']='1'; exec(open('errprop_globalfit_20261112.py').read()); kpc=ns['kpc']
F=pd.read_csv('bound_check_first_20261112.csv').set_index('gal')
FIX=('NGC0289','NGC5907')
def law_res(g,l):
    d=P0[P0.gal==g].sort_values('r'); m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values; ns['P']=d; ns['T']=pd.DataFrame([T0.loc[g]],index=[g])
    o=solve_gal(g,l); r=d.r.values; vo=d.vobs.values; lgo=np.log10(vo**2/r*1e6/kpc); return [float(x) for x in (lgo-np.log10(o.vmod2.values/r*1e6/kpc))[m]]
lines=open('pass2_table2_part.csv',newline='').read().splitlines(keepends=True); hdr=next(csv.reader([lines[0]]))
ix={k:i for i,k in enumerate(hdr)}; out=[]; n=0
for ln in lines:
    row=next(csv.reader([ln]))
    if row and row[0] in FIX:
        g=row[0]; row[ix['lrc']]=repr(float(F.lrc_first[g])); row[ix['chi_law']]=repr(float(F.chi_first[g])); row[ix['law']]=json.dumps(law_res(g,F.lrc_first[g]))
        b=io.StringIO(); csv.writer(b,lineterminator=ln[len(ln.rstrip('\r\n')):] or '\n').writerow(row); ln=b.getvalue(); n+=1
    out.append(ln)
open('pass2_table2_part.csv','w',newline='').write(''.join(out)); print('rows rewritten:',n)
R=pd.read_csv('pass2_table2_part.csv').drop_duplicates('gal'); hw=lambda x:(np.percentile(x,84)-np.percentile(x,16))/2
rms=np.sqrt(np.mean(np.concatenate([np.array(json.loads(x)) for x in R.law])**2))
print('N %d | median %.2f | std %.2f | robust %.2f | at bound %d | law rms %.4f'%(len(R),R.lrc.median(),R.lrc.std(ddof=1),hw(R.lrc),(R.lrc>=10.45).sum(),rms))
print(R.set_index('gal').loc[list(FIX),['lrc','chi_law']].round(3).to_string())
