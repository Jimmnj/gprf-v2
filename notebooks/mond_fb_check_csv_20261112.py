# mond_fb_check_csv_20261112.py -- DM-6 (ruled JM 2026-11-11 23:15): mond_fb_check_20261111.py rewritten to read notebooks/dms_fb_20261105.csv directly
# (ratified run output, md5 e6d35873...), replacing the five transcribed columns. Logic from "nu = {" on is VERBATIM from the 11/11 script. [UNRATIFIED until JM stamps]
# Pure Python + csv; no solver; nothing touches the law. Post-hoc check, not pre-registered.
# EXPECTATION, written 2026-11-12 01:26 EST (before the 01:26:48 run): every printed number equals the 11/11 footer to the printed precision
# (relation 0.695 / +0.091 / 0.130; RAR like-for-like 0.659 / +0.078 / 0.106; simple 0.656 / +0.075 / 0.105; standard 0.725 / +0.155 / 0.159;
# own masses 0.727 / 0.719 / 0.854; per-galaxy +0.005). The csv carries full precision, so a last-digit change is possible but none is expected.
import math, statistics as st, csv
a0k = 1.2e-10*3.0857e19/1e6   # (km/s)^2/kpc
with open("dms_fb_20261105.csv") as fh:
    rows = "\n".join("%s %s %s %s %s"%(d["gal"],d["hR"],d["Fb_dms"],d["Fb_pred"],d["V22"]) for d in csv.DictReader(fh))
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
