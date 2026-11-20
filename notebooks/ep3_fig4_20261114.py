# ep3_fig4_20261114.py -- EP-3: sec 3.3 / Fig 4 mu0 statistics with the test1 (2026-10-18, bounds 4-11) fits vs the V3 table of record. [UNRATIFIED]
# Runs fig4_mu0_mond_20261027.py's statistics block verbatim (no figure, no csv written); only F.l_tanh2 is swapped.
import sys, io, contextlib
src=open('fig4_mu0_mond_20261027.py').read().split("fig,ax=plt.subplots")[0]
key="F=pd.read_csv('test1_capped_baryon_20261018.csv'); T=load_table()"; assert src.count(key)==1
patch=key+"\nV3=pd.read_csv('pass2_table2_part.csv').set_index('gal').lrc\nif MODE=='v3': F['l_tanh2']=F.gal.map(V3).fillna(F.l_tanh2)  # 2 of 120 (UGC05005, UGCA444) are not in the 147 table: test1 value kept\n"
for MODE in ('t1','v3'):
    print('=====',MODE); exec(src.replace(key,patch),{'MODE':MODE,'__name__':'x'})
