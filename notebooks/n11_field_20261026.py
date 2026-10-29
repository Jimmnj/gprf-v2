# n11_field_20261026.py -- N11 pass 1 (mean field), per PREREG_N11_20261026.md (stamped 13:30; pre-run corrections sec 8).
# Method (a): each galaxy's isothermal excess truncated where it meets the conditional mean nb density around it;
# FLOOR-1 reading (ii): rho_f + H = rho_bar (global mean of ALL nb mass). Units: kpc, km/s, Msun.
import numpy as np, sys
G=4.30091e-6; Mpc=1000.0; kpc_m=3.0857e19
a0=1.2e-10*kpc_m/1e6            # (km/s)^2/kpc
rho_bar=35.5; rho_m=42.8; h=0.7
rho_c=10**7.45
def run(r0h=5.5, gam=1.85, Mcut=1e7, lumdep=False, L_xi=40*Mpc, verbose=True):
    r0=r0h/h*Mpc
    # GSMF Baldry+2012 (h=0.7), per dex per Mpc^3; M_b = 1.3 M*
    lMs=np.arange(np.log10(Mcut),12.0+1e-9,0.05); Ms=10**lMs; dlog=0.05
    Mstar=10**10.66; p1,a1,p2,a2=3.96e-3,-0.35,0.79e-3,-1.47
    phi=np.log(10)*np.exp(-Ms/Mstar)*(p1*(Ms/Mstar)**(a1+1)+p2*(Ms/Mstar)**(a2+1))
    n=phi*dlog/Mpc**3               # per kpc^3
    Mb=1.3*Ms; s2=0.5*np.sqrt(G*a0*Mb); sig=np.sqrt(s2)
    rc=3*sig/np.sqrt(4*np.pi*G*rho_c)
    # bias for luminosity-dependent clustering (Zehavi+2011 Table 2 thresholds, Mr -> log M* via L* <-> 10.66)
    if lumdep:
        Mr_t=np.array([-18.0,-18.5,-19.0,-19.5,-20.0,-20.5,-21.0,-21.5,-22.0]); r0_t=np.array([4.10,4.48,4.86,5.11,5.54,5.60,5.98,7.27,10.71])
        lM_t=10.66-0.4*(Mr_t+20.44)
        b=(np.interp(lMs,lM_t,r0_t)/r0h)**(gam/2)
    else: b=np.ones_like(Ms)
    def Xi(u,Lin=0.0):              # int_Lin^u s xi(s) ds, xi=(s/r0)^-gam for s<L_xi
        u=np.clip(u,Lin,L_xi); return r0**gam*(u**(2-gam)-Lin**(2-gam))/(2-gam)
    def excess(j,rl):               # returns t-grid, e(t), M_ex
        s2j=s2[j]; rtp2=s2j/(2*np.pi*G*rl)-rc[j]**2
        if rtp2<=0: return None
        rt=np.sqrt(rtp2); t=np.logspace(np.log10(rc[j]*1e-3),np.log10(rt),400)
        e=s2j/(2*np.pi*G*(t**2+rc[j]**2))-rl; e=np.clip(e,0,None)
        Mex=np.trapz(4*np.pi*t**3*e,np.log(t))
        return t,e,Mex,rt
    def conv(t,e,r,Lin=0.0):        # (xi conv e)(r) = 2pi/r int t e(t)[Xi(r+t)-Xi(|r-t|)] dt
        r=np.atleast_1d(r)[:,None]
        integ=t*e*(Xi(r+t,Lin)-Xi(np.abs(r-t),Lin))*t   # extra t for d ln t
        return 2*np.pi/r[:,0]*np.trapz(integ,np.log(t),axis=1)
    rl=np.full(len(Ms),rho_bar*1.0)
    for it in range(200):
        prof=[excess(j,rl[j]) for j in range(len(Ms))]
        rt=np.array([p[3] for p in prof])
        C0=np.zeros(len(Ms))       # sum_j n_j b_j (xi0 conv e_j)(rt_i)
        for j,p in enumerate(prof):
            C0+=n[j]*b[j]*conv(p[0],p[1],rt)
        new=rho_bar+b*C0
        if np.max(np.abs(np.log(new/rl)))<1e-5: break
        rl=np.sqrt(rl*new)
    Mex=np.array([p[2] for p in prof]); H=np.sum(n*Mex); rho_f=rho_bar-H
    out=dict(it=it,H=H,rho_f=rho_f,lMs=lMs,Mb=Mb,rl=rl,rt=rt,Mex=Mex,n=n,prof=prof,b=b,conv=conv,excess=excess,s2=s2)
    if verbose:
        w=n*Mex
        print(f'r0={r0h} h^-1 Mpc ({r0/Mpc:.2f} Mpc) gamma={gam} Mcut={Mcut:.0e} lumdep={lumdep}: iter {it}')
        print(f'  H (mean excess) = {H:.1f} Msun/kpc^3 = {H/rho_bar:.2f} rho_bar ; rho_f = {rho_f:.1f} ({rho_f/rho_bar:+.2f} rho_bar)  K1 {"PASS" if rho_f>=0 else "FAIL"}')
        print(f'  mass-weighted rho_loc = {np.sum(w*rl)/np.sum(w):.0f} ({np.sum(w*rl)/np.sum(w)/rho_bar:.1f} rho_bar)')
    return out
def isolated(o,lMb_bins=(10.10,10.66,10.96,11.29),L=4/h*Mpc):
    res=[]
    for lMb in lMb_bins:
        Mb=10**lMb; s2=0.5*np.sqrt(G*a0*Mb); sig=np.sqrt(s2); rc=3*sig/np.sqrt(4*np.pi*G*rho_c)
        rl=rho_bar
        for it in range(200):
            rt=np.sqrt(s2/(2*np.pi*G*rl)-rc**2)
            C=0.0
            for j,p in enumerate(o['prof']): C+=o['n'][j]*o['b'][j]*o['conv'](p[0],p[1],rt,Lin=L)[0]
            new=rho_bar+C
            if abs(np.log(new/rl))<1e-6: break
            rl=np.sqrt(rl*new)
        res.append((lMb,rl,rt))
    return res
if __name__=='__main__':
    print('=== N11 pass 1, FIDUCIAL ===')
    o=run()
    print('  per mass (log M*, M_b, rho_loc, r_t Mpc, M_ex/M_b, share of H):')
    w=o['n']*o['Mex']
    for k in range(0,len(o['lMs']),10):
        print(f'   {o["lMs"][k]:5.2f} {o["Mb"][k]:.2e} {o["rl"][k]:8.0f} {o["rt"][k]/Mpc:6.2f} {o["Mex"][k]/o["Mb"][k]:7.1f} {w[k]/w.sum():.3f}')
    print('  K2 isolated (xi=0 inside 5.7 Mpc): log M_b, rho_loc, r_t [Mpc]  (bound 160-355 fixed amplitude)')
    for lMb,rl,rt in isolated(o): print(f'   {lMb:5.2f} {rl:7.1f} {rt/Mpc:5.2f}  K2 {"PASS" if rl<160 else ("MARGINAL" if rl<355 else "FAIL")}')
    print('\n=== SENSITIVITY (prereg sec 3 / 8a) ===')
    for kw in [dict(r0h=4.7),dict(r0h=6.3),dict(gam=1.75),dict(gam=1.95),dict(Mcut=1e8),dict(lumdep=True),dict(L_xi=20*Mpc)]:
        run(**kw)
