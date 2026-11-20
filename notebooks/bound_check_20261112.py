# bound_check_20261112.py -- BOUND-1 (JM 2026-11-12 23:43: "0289 is very close to 10.5 if it is arbitrary"). [UNRATIFIED]
# Origin of the bound: bounds=(4.5,10.5) first appears in pass2_table2_full_20261022.py l.21 with no stated rationale -> treat as ARBITRARY.
# Solver + chi^2 verbatim via errprop_globalfit_20261112.chi_fn (rule 3a: same chi_fn that reproduced 147/147).
# EXPECTATIONS (written ~23:44, before any run; stamp corrected):
#  E1 NGC 0289: interior minimum near 10.37; chi^2 rises past it; nothing beyond 10.5 beats 13.3.
#  E2 NGC 5907: chi^2 keeps falling past 10.5, then either a minimum by ~12 or a flat floor (saturation).
#  E3 all 147: <= 3 galaxies find a lower chi^2 beyond 10.5 (scan 10.5-13.0 at 0.1).
import os, sys, time, numpy as np, pandas as pd
MODE=sys.argv[1] if len(sys.argv)>1 else 'two'   # captured before exec (errprop_20261112 overwrites sys.argv)
os.environ['NOFIT']='1'
exec(open('errprop_globalfit_20261112.py').read())
EXT=np.round(np.arange(4.5,13.0001,0.1),2)
if __name__=='__main__':
    mode=MODE
    if mode=='two':
        for g in ['NGC0289','NGC5907']:
            f=chi_fn(g); t=time.time(); c=[f(l) for l in EXT]
            print(g,'(%.1f s)'%(time.time()-t)); print(' '.join('%.1f:%.1f'%(l,x) for l,x in zip(EXT,c) if (l*10)%5==0 or l>=10))
    if mode=='all':
        HI=np.round(np.arange(10.5,13.0001,0.1),2); G=pd.read_csv('errprop_globalfit_20261112_part.csv').set_index('gal')
        PART='bound_check_20261112_part.csv'; done=set(pd.read_csv(PART).gal) if os.path.exists(PART) else set(); t0=time.time()
        for g in gals:
            if g in done: continue
            if time.time()-t0>140: break
            f=chi_fn(g); c=np.array([f(l) for l in HI]); i=int(np.argmin(c)); x,cx=HI[i],c[i]
            if cx<G.loc[g,'chi_g']-1e-9:
                r=minimize_scalar(f,bounds=(max(10.5,x-0.1),min(13.0,x+0.1)),method='bounded',options={'xatol':0.005})
                if r.fun<cx: x,cx=r.x,r.fun
            pd.DataFrame([dict(gal=g,lrc_g=G.loc[g,'lrc_g'],chi_g=G.loc[g,'chi_g'],lrc_hi=x,chi_hi=cx,beats=int(cx<G.loc[g,'chi_g']-1e-9),edge13=int(i==len(HI)-1))]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False)
        print('progress %d / %d (%.0f s)'%(len(pd.read_csv(PART)),len(gals),time.time()-t0))
# OUTPUT two (run 23:44-23:47; earlier stamp guessed, corrected): NGC0289 interior min ~10.4 (chi2 13.3), rises to 20.1 at 12.3, eases to 19.3 at 13.0; low branch 7.0 (22.4). E1 HOLDS.
#   NGC5907 min ~10.85 (chi2 5.6) -- BEYOND the 10.5 bound; rises to 9.5 at 12.8; low branch 7.5 (14.5). E2 holds (interior minimum by 12).
# OUTPUT all (run 23:44-23:47; stamp corrected): 3 galaxies find lower chi^2 in 10.5-13.0: NGC5907 10.50->10.84 (6.13->5.61); NGC2998 8.05->11.35 (18.03->16.26); UGC11914 8.63->12.21 (25.68->16.93). One galaxy's 10.5-13 scan minimum sits at the 13.0 edge (does not beat its global chi^2). E3 holds at the limit (3) but TWO MID-RANGE FITS JUMP by 3.3-3.6 dex.
# E4 (written before the far run; stamp corrected): beyond 13 (to 16) chi^2 flattens (saturation) and no galaxy finds a lower minimum than its 4.5-13 one.
    if mode=='far':
        B=pd.read_csv('bound_check_20261112_part.csv').set_index('gal'); FAR=np.round(np.arange(12.0,16.0001,0.25),2)
        for g in list(B.index[B.edge13==1])+['NGC0289','NGC2998','UGC11914','NGC5907']:
            f=chi_fn(g); c=[f(l) for l in FAR]; print(g,'best(4.5-13) %.2f @ %.2f |'%(min(B.loc[g,'chi_g'],B.loc[g,'chi_hi']),B.loc[g,'lrc_hi'] if B.loc[g,'beats'] else B.loc[g,'lrc_g']),' '.join('%.2f:%.1f'%(l,x) for l,x in zip(FAR,c)))
# OUTPUT far (run 23:44-23:47; stamp corrected): E4 holds for the 5 scanned (no lower minimum in 12-16), but chi^2 is NOT flat: it oscillates in log rho_c.
#   NGC0289 max ~12.3, min ~14.5; NGC5907 max ~12.8, min ~15.0; NGC2998 max ~13.5; UGC11914 max ~14.25. Second dips shallower than the first (damped).
#   Spread stats, all 147: grid 4.5-10.5 median 7.415 robust 0.497 std 0.628 | grid 4.5-13 median 7.415 robust 0.497 std 0.805 (movers all above p84 7.985).
# HYPOTHESIS (written before the fine scan, ~23:47; stamp corrected) [T3]: the oscillation is the isothermal sphere's homology oscillation. Lane-Emden: psi ~ ln(xi^2/2) + A xi^(-1/2) cos((sqrt7/2) ln xi + d)
#   -> period 4 pi/sqrt7 = 4.750 in ln r = 2.063 dex in r; r0 = sigma/sqrt(4 pi G rho_c) ~ rho_c^(-1/2) -> PREDICTED PERIOD 4.13 dex in log rho_c.
    if mode=='fine':
        F=np.round(np.arange(9.0,17.0001,0.05),2)
        for g in ['NGC0289','NGC5907','NGC2998','UGC11914']:
            f=chi_fn(g); c=np.array([f(l) for l in F])
            mn=[F[k] for k in range(1,len(c)-1) if c[k]<c[k-1] and c[k]<=c[k+1]]; mx=[F[k] for k in range(1,len(c)-1) if c[k]>c[k-1] and c[k]>=c[k+1]]
            print(g,'minima',mn,'maxima',mx)
# OUTPUT fine: min-to-min spacing NGC0289 10.35->14.40 (4.05); NGC5907 10.85->14.90 (4.05); NGC2998 11.35->15.50 (4.15); UGC11914 12.20->16.35 (4.15);
#   max-to-max NGC2998 9.40->13.45 (4.05), UGC11914 10.20->14.30 (4.10). PREDICTED 4.13. HYPOTHESIS SUPPORTED [T3]. Pre-asymptotic first gap shorter (0289 7.0->10.35, 5907 7.5->10.85: 3.35).
#   Wiggles in 0289 at 16.25-16.65: numerical, ignore.
# FOUND: marg_ns_20261030.py arm A prior on log rho_c is flat [6,10] -> arm A excludes the high branch for all four galaxies.
    if mode=='first':   # V3: first (lowest-rho_c) local minimum, for galaxies with >1 grid minimum in [4.5,10.5] or a lower minimum in 10.5-13
        B=pd.read_csv('bound_check_20261112_part.csv').set_index('gal'); G=pd.read_csv('errprop_globalfit_20261112_part.csv').set_index('gal')
        sel=[g for g in B.index if G.loc[g,'n_localmin']>1 or B.loc[g,'beats']==1]; out=[]
        for g in sel:
            f=chi_fn(g); c=np.array([f(l) for l in GRID]); k=next(k for k in range(1,len(c)-1) if c[k]<c[k-1] and c[k]<=c[k+1])
            r=minimize_scalar(f,bounds=(GRID[k]-0.05,GRID[k]+0.05),method='bounded',options={'xatol':0.005}); out.append(dict(gal=g,lrc_first=r.x,chi_first=r.fun))
        pd.DataFrame(out).to_csv('bound_check_first_20261112.csv',index=False); print(pd.DataFrame(out).round(3).to_string())
