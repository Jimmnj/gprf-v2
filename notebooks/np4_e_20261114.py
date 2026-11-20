# np4_e_20261114.py [UNRATIFIED] Item 4, part E: V3 floor check on the table of record. For each galaxy: chi2 at the scan floor (log rho_c 4.5) against
# chi2 at the fitted first minimum (ep3_val_20261114_part.csv). A floor chi2 BELOW the first interior minimum would mean the first-minimum rule skipped a
# floor fit (the NGC 4085 case at M/L 0.7, np4_c). Objective: ep3_audit_20261114.chi_drop (verbatim from clean_conc_20261109.refit), exec'd from that file.
import sys, numpy as np, pandas as pd
src=open('ep3_audit_20261114.py').read().split("PART='ep3_")[0]; sys.argv=['x','none']; exec(src)
V=pd.read_csv('ep3_val_20261114_part.csv').set_index('key'); rows=[]
for g in V.index:
    f=chi_drop(g,0)
    if f is None: continue
    rows.append(dict(gal=g,lrc=V.lrc_v3[g],chi_min=V.chi[g],chi_floor=f(4.5),chi_5=f(5.0)))
R=pd.DataFrame(rows).set_index('gal'); R.to_csv('np4_e_floor_20261114.csv')
print('N %d | floor chi2 < chi2 at the first minimum: %d | min ratio chi_floor/chi_min %.2f (%s) | chi(5.0) < chi(4.5) in %d (falling from the floor)'%(len(R),(R.chi_floor<R.chi_min).sum(),(R.chi_floor/R.chi_min).min(),(R.chi_floor/R.chi_min).idxmin(),(R.chi_5<=R.chi_floor).sum()))
print(R[(R.chi_floor<R.chi_min*1.5)].round(3).to_string())
