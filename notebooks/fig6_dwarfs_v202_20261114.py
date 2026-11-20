# fig6_dwarfs_v202_20261114.py -- runs make_figs_v2_20261020.py's header and Fig 6 block verbatim on the F1 rerun csv
# (np4_dwarf_prereg_C_20261114.csv: Plummer scale radius = projected R_e, ruling DF1-F1 2026-11-14 17:34). Only the csv name differs.
import sys; sys.argv=['x','fig6']
src=open('make_figs_v2_20261020.py').read(); head=src.split('# Fig 1')[0]; f6=src.split('# Fig 6')[1]
key="pd.read_csv('dwarf_prereg_20261020.csv')"; assert f6.count(key)==1 and src.count('# Fig 6')==1
exec(head+'# Fig 6'+f6.replace(key,"pd.read_csv('np4_dwarf_prereg_C_20261114.csv')"),{'__name__':'x'})
