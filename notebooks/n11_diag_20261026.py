# n11_diag_20261026.py -- EXPLORATORY diagnostics of N11 pass 1 (NOT pre-registered). Brackets the mean-field (Jensen) bias.
import numpy as np, warnings; warnings.filterwarnings('ignore')
import n11_field_20261026 as N
o=N.run(verbose=False)
n,Mex,lMs,s2=o['n'],o['Mex'],o['lMs'],o['s2']
def Mex_at(j,rl):
    p=o['excess'](j,rl); return p[2] if p else 0.0
# (1) convergence check for the two capped runs: residual at the last iterate
for kw in [dict(gam=1.95),dict(lumdep=True)]:
    a=N.run(verbose=False,**kw)
    print(kw,'iter',a['it'],'H',round(a['H'],2),'rho_f',round(a['rho_f'],2))
# (2) H by mass, fiducial mean field vs NO clustering truncation (rho_loc = rho_bar for all)
H_unc=np.array([n[j]*Mex_at(j,N.rho_bar) for j in range(len(n))])
H_mf=n*Mex
print('\nlog M*   share_unclustered  H_unc_cum   H_mf_cum')
for cut in [8,9,9.5,10,10.5,11,11.5,12.01]:
    m=lMs<cut; print(f'  <{cut:5.2f}  {H_unc[m].sum()/H_unc.sum():6.3f}   {H_unc[m].sum():7.1f}   {H_mf[m].sum():7.1f}')
print(f'TOTAL unclustered H = {H_unc.sum():.1f} ({H_unc.sum()/N.rho_bar:.2f} rho_bar)  vs mean field {H_mf.sum():.1f}')
# (3) bracket: mean field above a mass split, untruncated (rho_bar) below it -- i.e. 'typical dwarf sees no neighbour'
for split in [9.0,9.5,10.0,10.5]:
    m=lMs<split; H=H_unc[m].sum()+H_mf[~m].sum()
    print(f'  dwarfs below logM*={split} untruncated: H={H:6.1f} ({H/N.rho_bar:.2f} rho_bar) rho_f={N.rho_bar-H:+6.1f}  K1 {"PASS" if H<=N.rho_bar else "FAIL"}')
# (4) where does the conditional density at a dwarf's r_t come from: 1+xi at that radius
r0=5.5/0.7*N.Mpc
for j in [0,20,40,60,80]:
    rt=o['rt'][j]; print(f'  logM* {lMs[j]:.2f}: r_t {rt:7.1f} kpc, xi(r_t)={(rt/r0)**-1.85:.2e}')
