# NP-17 side check, 2026-11-11. Input: five columns transcribed from notebooks/dms_fb_20261105.csv (ratified run).
# [UNRATIFIED] Run by Claude in the cloud workspace 2026-11-11 13:19 EST (pure Python, no solver). Post-hoc check, not pre-registered.
# Expectation declared BEFORE running (Claude, 13:2x): like-for-like MOND median F_b 0.66-0.74 (RAR function);
# MOND at the observed V_c (its own masses) median 0.75-0.85.
# Columns below: gal, hR [kpc], Fb_dms (V_b/V_c, DMS), Fb_pred (V_b^2/V_c^2, eq. 1), V22 (model V_tot at 2.2 hR, km/s).
# A version that reads the csv directly should replace this one if the check is ratified.
import math, statistics as st
a0k = 1.2e-10*3.0857e19/1e6   # (km/s)^2/kpc
rows = """448 3.9 0.57 0.50170 159.571
463 3.8 0.72 0.58053 193.644
1081 3.1 0.63 0.51301 130.276
1087 3.2 0.48 0.48690 126.710
1529 3.6 0.55 0.51403 154.193
1635 2.9 0.45 0.47176 106.989
1862 1.4 0.60 0.59708 77.738
1908 4.9 0.68 0.53945 229.070
3091 3.6 0.41 0.32275 128.800
3140 3.5 0.68 0.62725 205.533
3701 3.6 0.55 0.41397 120.197
3997 5.5 0.48 0.29842 136.310
4036 4.3 0.58 0.45221 157.326
4107 3.2 0.56 0.50139 129.794
4256 4.7 0.73 0.51903 210.129
4368 3.2 0.72 0.58504 168.743
4380 5.0 0.46 0.38565 156.361
4458 9.0 0.68 0.39509 245.831
4555 4.1 0.59 0.47946 159.505
4622 7.6 0.52 0.36313 214.447
6903 4.2 0.34 0.27099 101.094
6918 1.2 0.64 0.78428 105.416
7244 3.9 0.56 0.26617 126.817
7917 8.5 0.56 0.36880 226.688
8196 4.9 1.09 0.73332 373.776
9177 7.0 0.51 0.36007 192.452
9837 5.8 0.44 0.21423 117.394
9965 3.5 0.53 0.44144 145.258
11318 4.5 0.66 0.53175 202.974
12391 3.9 0.70 0.50837 174.112"""
nu = {'rar': lambda y: 1/(1-math.exp(-math.sqrt(y))),
      'simple': lambda y: 0.5+math.sqrt(0.25+1/y),
      'standard': lambda y: math.sqrt(0.5+math.sqrt(0.25+1/y**2))}
def gN_from_gobs(x, f):   # solve nu(y)*y = x for y by bisection
    lo, hi = 1e-8, x
    for _ in range(200):
        mid = 0.5*(lo+hi)
        if f(mid)*mid > x: hi = mid
        else: lo = mid
    return 0.5*(lo+hi)
R=[]
for l in rows.splitlines():
    g,hR,Fd,Fp2,V22 = l.split(); hR,Fd,Fp2,V22 = map(float,(hR,Fd,Fp2,V22))
    r = 2.2*hR; vb2 = Fp2*V22**2; yN = vb2/r/a0k
    Vc_obs = math.sqrt(vb2)/Fd; x_obs = Vc_obs**2/r/a0k
    d = dict(gal=g, Frel=math.sqrt(Fp2), Fd=Fd, yN=yN, xobs=x_obs)
    for k,f in nu.items():
        d['like_'+k] = 1/math.sqrt(f(yN))                       # same DMS masses in, F_b out
        d['own_'+k]  = math.sqrt(gN_from_gobs(x_obs,f)/x_obs)   # V_b,N/V_c MOND needs at the observed V_c
    R.append(d)
med = lambda k,S=R: st.median(x[k] for x in S)
rms = lambda k,S=R: math.sqrt(sum((x[k]-x['Fd'])**2 for x in S)/len(S))
print("N =",len(R)," a0 = %.1f (km/s)^2/kpc"%a0k)
print("transcription check: relation median F %.3f, median(F-Fdms) %+.3f, rms %.3f   [file: 0.66, +0.091, 0.130]"%(med('Frel'), st.median(x['Frel']-x['Fd'] for x in R), rms('Frel')))
print("DMS median F_b %.2f ; median g_N/a0 at 2.2hR (DMS masses) %.2f ; median g_obs/a0 %.2f"%(med('Fd'),med('yN'),med('xobs')))
for k in nu:
    print("%-9s like-for-like: median F %.3f, median(F-Fdms) %+.3f, rms %.3f | own masses: median F %.3f, median(F-Fdms) %+.3f"%(
        k, med('like_'+k), st.median(x['like_'+k]-x['Fd'] for x in R), rms('like_'+k), med('own_'+k), st.median(x['own_'+k]-x['Fd'] for x in R)))
S=[x for x in R if x['gal']!='8196']
print("excluding UGC 8196:")
for k in nu:
    print("%-9s like-for-like: median(F-Fdms) %+.3f, rms %.3f | own: median F %.3f"%(k, st.median(x['like_'+k]-x['Fd'] for x in S), rms('like_'+k,S), med('own_'+k,S)))
print("median (F_rel - F_like_rar) per galaxy: %+.3f"%st.median(x['Frel']-x['like_rar'] for x in R))
print("eta implied by F=1 (1/median Fdms) = %.2f ; F implied by eta 1.34 = %.2f, by eta 1.41 = %.2f"%(1/med('Fd'), 1.34*med('Fd'), 1.41*med('Fd')))
# OUTPUT 2026-11-11 13:19 EST:
# N = 30  a0 = 3702.8 (km/s)^2/kpc
# transcription check: relation median F 0.695, median(F-Fdms) +0.091, rms 0.130   [file: 0.66, +0.091, 0.130]
# DMS median F_b 0.56 ; median g_N/a0 at 2.2hR (DMS masses) 0.32 ; median g_obs/a0 1.07
# rar       like-for-like: median F 0.659, median(F-Fdms) +0.078, rms 0.106 | own masses: median F 0.727, median(F-Fdms) +0.127
# simple    like-for-like: median F 0.656, median(F-Fdms) +0.075, rms 0.105 | own masses: median F 0.719, median(F-Fdms) +0.120
# standard  like-for-like: median F 0.725, median(F-Fdms) +0.155, rms 0.159 | own masses: median F 0.854, median(F-Fdms) +0.250
# excluding UGC 8196:
# rar       like-for-like: median(F-Fdms) +0.078, rms 0.102 | own: median F 0.723
# simple    like-for-like: median(F-Fdms) +0.076, rms 0.099 | own: median F 0.716
# standard  like-for-like: median(F-Fdms) +0.157, rms 0.160 | own: median F 0.851
# median (F_rel - F_like_rar) per galaxy: +0.005
# eta implied by F=1 (1/median Fdms) = 1.77 ; F implied by eta 1.34 = 0.76, by eta 1.41 = 0.80
