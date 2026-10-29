# body_rule_test_20261022.py [UNRATIFIED] 2026-10-22 12:30 (JM: a gas cloud is a body as a whole, with a diffuse-core halo; molecules are not).
# Two candidate criteria per system: (H) Hill/tidal radius r_J > own size R; (B) bound: random speed s < escape-like speed v_e = sqrt(G M/R)
# [current sec 2.3 rule, baryons only], and (B*) bound once the system's OWN nb (law, BTFR sigma, isothermal 2 sigma^2 R/G inside R) is included --
# the bootstrap fixed point. r_J: flat host r_J = (G m D^2/(2 v^2))^(1/3); for a dominant body, against the ambient: r = (3 m/(4 pi rho_amb))^(1/3).
import numpy as np
G=4.30091e-6; a0k=1.2e-10*3.0857e19/1e6; ramb=35.5; kpc_m=3.0857e19; Msun_kg=1.989e30
def rJ_flat(m,D,v): return (G*m*D**2/(2*v**2))**(1/3)
def rJ_amb(m): return (3*m/(4*np.pi*ramb))**(1/3)
def nb_in(m,R): s2=np.sqrt(G*a0k*m)/2; return 2*s2*R/G, np.sqrt(s2)
cases=[ # name, m (Msun), R (kpc), random speed (km/s), host: ('flat',D,v) or ('amb',)
 ('H2 molecule, MW at 8 kpc',2*1.67e-27/Msun_kg,1e-10/kpc_m,0.3,('flat',8.,230.)),
 ('H2 molecule, Earth surface',2*1.67e-27/Msun_kg,1e-10/kpc_m,0.5,('pt',6.371e6/kpc_m,3.0e-6)),
 ('proton in Coma ICM at 500 kpc',1.67e-27/Msun_kg,1e-15/kpc_m,1200.,('flat',500.,1400.)),
 ('Coma ICM + galaxies (whole)',1.5e14,2000.,1200.,('amb',)),
 ('GMC 1e6, R 30 pc, MW 8 kpc',1e6,0.030,7.,('flat',8.,230.)),
 ('small cloud 1e3, R 2 pc, MW 8 kpc',1e3,0.002,1.,('flat',8.,230.)),
 ('DDO 154 (gas-rich dwarf, isolated)',4e8,8.,10.,('amb',)),
 ('Cloud-9 (gas 1.4e6, R 1.4, 90 kpc from M94)',1.4e6,1.4,5.1,('flat',90.,180.)),
 ('Cloud-9 at 3D distance 200 kpc',1.4e6,1.4,5.1,('flat',200.,180.)),
]
print('%-44s %9s %9s %5s | %7s %7s %3s | %9s %8s %3s'%('system','R','r_J','H?','s','v_e(b)','B?','M_nb(<R)','v_e(b+nb)','B*?'))
for n,m,R,s,h in cases:
    if h[0]=='flat': rj=rJ_flat(m,h[1],h[2])
    elif h[0]=='pt': rj=h[1]*(m/(3*h[2]))**(1/3)
    else: rj=rJ_amb(m)
    ve=np.sqrt(G*m/R); mnb,sg=nb_in(m,min(R,rj)); ve2=np.sqrt(G*(m+mnb)/R)
    print('%-44s %9.2e %9.2e %5s | %7.1f %7.2e %3s | %9.2e %8.2e %3s'%(n,R,rj,'yes' if rj>R else 'no',s,ve,'yes' if s<ve else 'no',mnb,ve2,'yes' if s<ve2 else 'no'))
print('(R and r_J in kpc; speeds km/s. Molecule/proton rows: nb columns are formal, the law is not meant to apply at particle scale.)')
