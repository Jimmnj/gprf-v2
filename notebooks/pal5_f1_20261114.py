# pal5_f1_20261114.py [UNRATIFIED] K4 (JM 2026-11-14 21:04 "rerun pal 5"): Palomar 5 nb fraction inside r_h under DF1-F1 (Plummer a = projected r_h).
# RULE 3a: canon sec 6 R1, R1/R6 (ambient excess), R6, NEST-1, ONE-FIELD, A' read at boot (19:16). Nothing re-typed; machinery exec'd verbatim:
#  (A) RECORD METHOD: test5_gc_nb_fraction_20261018.py exec'd with its GCS line replaced by the Baumgardt values of the TEST5 addendum (2026-10-19:
#      M 1.34e4, r_h 20.42 pc projected half-light, R_gc 17.27 kpc) and lrc 7.42 added. Must reproduce the addendum's Pal 5 0.25 / 0.62 at 10^7.45 / 10^8.0.
#      That method adds the host ambient boosted by e^2 (the cap) as a worst case. Variant F1: "a=rh/1.305" -> "a=rh" (the only change).
#  (B) DWARF PROTOCOL: n9_dwarfs_20261003.py up to 'def run(' and df1_trail_20261114.py sig_one, exec'd verbatim (as DF-1). First reproduce the DF-1 record
#      (DF4 old inputs 1.5e8, rh 1.6, host 250 km/s at 80 kpc: 17.31 km/s, record 17.3). A': Pal 5 solved in its own phi_b with its own sigma; the host
#      enters only as ambient density (MW v_c 200 at R_gc, as TEST5) + rho_bg; no e^2 boost (NEST-1: exp(phi_body/sigma_host^2) = 1.000 for Pal 5, tex 2.3).
#      n9 omits the saturation S (x_b < 0.5 here; A.2 states the same caveat for the dSphs).
# OBSERVABLE (TEST5): nb/M* inside the 3-D radius r_h. (B) reports the body's excess (M_tot - M_b) and excess + ambient pass-through (4/3 pi r_h^3 rho_amb).
# Bar 0.1 (TEST5); sigma-equivalent = nb/M* / 0.179 (catalogue mass error 0.24/1.34, as the addendum's "CORRECTION 2").
import io, contextlib, numpy as np
OUT='pal5_f1_20261114_report.txt'; open(OUT,'w').close()
def say(s): print(s); open(OUT,'a').write(s+'\n')
say('=== (A) record method, test5 verbatim with the addendum inputs')
src=open('test5_gc_nb_fraction_20261018.py').read()
old="GCS=[('47 Tuc',7e5,4.2e-3,7.4),('omega Cen',3.5e6,7.6e-3,6.4),('NGC 2419',9e5,2.4e-2,90.),('Pal 5',1.5e4,2.0e-2,18.)]"
assert src.count(old)==1 and src.count('for lrc in (7.45,8.0):')==1 and src.count('a=rh/1.305')==1
srcA=src.replace(old,"GCS=[('Pal 5',1.34e4,2.042e-2,17.27)]").replace('for lrc in (7.45,8.0):','for lrc in (7.42,7.45,8.0):')
for tag,code in (('P (record convention, a = r_h/1.305)',srcA),('C (DF1-F1, a = r_h)',srcA.replace('a=rh/1.305','a=rh'))):
    buf=io.StringIO()
    with contextlib.redirect_stdout(buf): exec(code,{})
    say('-- '+tag); [say('   '+l) for l in buf.getvalue().strip().split('\n')]
say('=== (B) dwarf protocol (n9 + df1 sig_one verbatim)')
ns={}; exec(open('n9_dwarfs_20261003.py').read().split('def run(')[0],ns)
d1=open('df1_trail_20261114.py').read(); i=d1.index('def sig_one('); j=d1.index('def stars_only('); exec(d1[i:j],ns)
G,a0k,rho_bg,host_rho,solve=(ns[k] for k in ('G','a0k','rho_bg','host_rho','solve'))
v=ns['sig_one']((1.5e8,1.6/1.305,1.6,host_rho(250.,80.)+rho_bg,7.45))
say('REPRODUCE: DF4 old inputs at 10^7.45, host: %.2f km/s (record 17.3; DF-1 rep 17.31)'%v); assert abs(v-17.31)<0.02
Mb,rh,R=1.34e4,2.042e-2,17.27; amb=host_rho(200.,R)+rho_bg; s2=0.5*np.sqrt(G*a0k*Mb)
say('Pal 5: M_b %.2e, r_h %.1f pc, R_gc %.2f kpc, sigma %.2f km/s, rho_amb 10^%.2f'%(Mb,rh*1e3,R,np.sqrt(s2),np.log10(amb)))
def frac(a,lrc):
    rmax=max(8*1.305*a,3*np.sqrt(s2/(2*np.pi*G*amb))); s=solve(Mb,a,s2,amb,lrc,rmax)
    Ms=Mb*rh**3/(rh**2+a**2)**1.5; exc=s.sol(rh)[1]-Ms; return exc/Ms,(exc+4/3*np.pi*rh**3*amb)/Ms,Ms
for tag,a in (('P (a = r_h/1.305)',rh/1.305),('C (DF1-F1, a = r_h)',rh)):
    x0=G*Mb/a/s2; say('-- %s: x_b at centre %.2f'%(tag,x0))
    for lrc in (7.42,7.45,8.0):
        e,t,Ms=frac(a,lrc); say('   rho_c 10^%.2f: M*(<r_h) %.0f | nb/M* excess %.3f, + ambient %.3f | %.1f sigma-equiv'%(lrc,Ms,e,t,t/0.179))
    lo,hi=5.,9.
    for _ in range(40):
        m=0.5*(lo+hi); lo,hi=(m,hi) if frac(a,m)[1]<0.1 else (lo,m)
    say('   bar nb/M* = 0.1 at rho_c 10^%.2f; at 1 sigma-equiv (0.179): see the row values'%m)
