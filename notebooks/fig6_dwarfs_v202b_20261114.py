# fig6_dwarfs_v202b_20261114.py -- K5 (JM 2026-11-14 21:33 "opens accepted then batch 8"; BALLOT RULING 21:33): label offsets only.
# Runs make_figs_v2_20261020.py's header and Fig 6 block verbatim on the F1 rerun csv, as fig6_dwarfs_v202_20261114.py does, with one change:
# the annotate call takes a per-name offset and alignment from OFF; names not in OFF keep the old (3,3) offset and default alignment.
# With OFF empty the output is byte-identical to fig6_dwarfs_v202_20261114.py's (md5 594dafb5), checked 2026-11-14 before this run.
import sys; sys.argv=['x','fig6']
OFF={'Tucana':((-4,3),'right','bottom'), 'AndXIX':((-4,3),'right','bottom'),
     'UrsaMinor':((4,-3),'left','top'), 'LeoI':((-4,3),'right','bottom'), 'AndXXI':((3,-3),'left','top')}
src=open('make_figs_v2_20261020.py').read(); head=src.split('# Fig 1')[0]; f6=src.split('# Fig 6')[1]
key="pd.read_csv('dwarf_prereg_20261020.csv')"; assert f6.count(key)==1 and src.count('# Fig 6')==1
ann="for _,r in w.iterrows(): ax.annotate(r['name'],(r.s_obs,r.s_GPRF),fontsize=6,xytext=(3,3),textcoords='offset points')"
new="for _,r in w.iterrows(): o,ha,va=OFF.get(r['name'],((3,3),'left','baseline')); ax.annotate(r['name'],(r.s_obs,r.s_GPRF),fontsize=6,xytext=o,textcoords='offset points',ha=ha,va=va)"
assert f6.count(ann)==1
exec(head+'# Fig 6'+f6.replace(key,"pd.read_csv('np4_dwarf_prereg_C_20261114.csv')").replace(ann,new),{'__name__':'x','OFF':OFF})
