# df1_gc_20261114.py [UNRATIFIED] DF-1 follow-up (JM 2026-11-14 16:03 "even if the rho_c is a fitted number does the rest of the curve fit"; 16:04 "df2 4 9").
# With rho_c FITTED to the stellar dispersion inside R_e (df1_trail main, convention C, FIELD), what does the same solution predict further out?
# (i) enclosed nb vs stellar mass at 1, 2, 4 R_e, 10 kpc and at r_t; (ii) GC dispersion: number-weighted projected sigma_los of a Plummer tracer
# (a_GC = 1.0, 1.5, 2.0 R_e) inside R <= R_max, in the same potential (n9 sigma_los verbatim, tracer scale passed as its 'a').
# GC data [T1]: DF2 7.8 +5.2/-2.2 (10 GCs, vD+2018 RNAAS) and 10.6 +3.9/-2.3 (7 GCs + 3 PNe, Emsellem+2019); GCs to R = 7.6 kpc at 20 Mpc.
#               DF4 4.2 +4.4/-2.2 (7 GCs, vD+2019); GCs to >= 7 kpc at 19.9 Mpc. DF9: no GC dispersion used.
import sys, warnings, numpy as np; warnings.filterwarnings('ignore')
sys.argv=['x','none']; exec(open('df1_trail_20261114.py').read().split("if STAGE=='rep':")[0])
def prof(Mb,a,lrc):
    s2=0.5*np.sqrt(G*a0k*Mb); amb=FIELD; rt=np.sqrt(s2/(2*np.pi*G*amb)); rmax=max(8*1.305*a,3*rt)
    s=solve(Mb,a,s2,amb,lrc,rmax); return s,s2,rt,rmax
Ms=lambda Mb,a,r: Mb*r**3/(r*r+a*a)**1.5
OBJ=[('DF2',2.0e8*0.88**2,2.0*0.88,7.6*17.6/20,[(5.913,'fit'),(6.415,'+1sd bound')],'7.8 +5.2/-2.2 (vD18) | 10.6 +3.9/-2.3 (Emsellem19)'),
     ('DF4',1.5e8*0.94**2,1.6*0.94,7.0*18.8/19.9,[(5.865,'fit'),(6.478,'+1sd bound')],'4.2 +4.4/-2.2 (vD19)'),
     ('DF9',1.42e8*(18.7/20.6)**2,1.1*18.7/20.6,np.nan,[(6.735,'+1sd bound')],'none used')]
say('=== GC / outer-profile check (df1_gc_20261114.py) %s: rho_c fitted to the stars inside R_e; convention C; FIELD ambient'%__import__('time').strftime('%Y-%m-%d %H:%M'))
for nm,Mb,Re,Rmax,lrcs,gcobs in OBJ:
    sN,s2,rt,rmaxN=prof(Mb,Re,-99)
    for lrc,lab in lrcs:
        s,s2,rt,rmax=prof(Mb,Re,lrc); r0=np.sqrt(9*s2/(4*np.pi*G*10**lrc))
        Mt=lambda r: s.sol(min(r,rmax))[1]
        rows=[]
        for r in (Re,2*Re,4*Re,10.,rt):
            mt=Mt(r); ms=Ms(Mb,Re,r); rows.append('%.3g kpc: nb %.2e / stars %.2e (x%.1f)'%(r,mt-ms,ms,(mt-ms)/ms))
        say('%s log rho_c %.2f (%s): core radius r0 = %.1f kpc, r_t = %.0f kpc | '%(nm,lrc,lab,r0,rt)+' ; '.join(rows))
        if np.isfinite(Rmax):
            gc=[]; 
            for f in (1.0,1.5,2.0):
                g=sigma_los(s,f*Re,Rmax,rmax); g0=sigma_los(sN,f*Re,Rmax,rmaxN); gc.append('a_GC %.1f R_e: %.1f (stars only %.1f)'%(f,g,g0))
            say('      GC sigma_los inside R <= %.1f kpc: '%Rmax+' | '.join(gc)+' || observed '+gcobs)
        say('      v_c: '+', '.join('%.0f kpc %.1f'%(r,np.sqrt(G*Mt(r)/r)) for r in (Re,5.,10.,30.))+' km/s')
