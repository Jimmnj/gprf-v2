# gate2b_axisym_20261004.py  [UNRATIFIED]
# GATE 2b: axisymmetric canon-law Poisson-Boltzmann solve vs SPARC.
# Law (SNAP_20260922_001 sec 0c; same form as pb_solve_canon_20260925.py):
#   rho_nb = rho_c exp(-psi/sigma^2) - rho_amb,  psi(R,z) = Phi_N(R,z) - Phi_N(0) >= 0 (Newtonian sign)
# Baryons (ASSUMPTIONS): stellar disk Sigma = Yd*SBdisk(SPARC), exponential vertical, h_z = 0.196 Rd^0.633
#   (Bershady+2010); gas = exponential disk, M = 1.33 MHI, Rg fitted to Vgas (Freeman), h_g = 0.2 kpc;
#   bulge spherical from Vbul. Phi_b(R,z) by ring sum (complete elliptic integral K), softening eps = h/2.
# nb: monopole by outward Heun shooting of the angle-averaged law; l = 2..LMAX (even) by Green's function,
#   damped fixed-point iteration, warm-started across rho_c evaluations.
# Fit: in-plane v^2 = Vbar^2(SPARC) + R dPhi_nb/dR (z=0). One free parameter per galaxy (log rho_c), 10% weights.
# Modes: C0   spherical Phi_b from SPARC Vbar, monopole -> must reproduce pb_solve_canon (ballot G1 V0)
#        C1   real flattened Phi_b, nb monopole only
#        full real flattened Phi_b, l <= 12
# Kill criteria fixed BEFORE running: notebooks/BALLOT_20260926.md block G (G1).
# Usage (from notebooks/): python gate2b_axisym_20261004.py C0,C1,full [ideal|vflat] [--workers 16] [--gal A,B] [--amb 1e3]
import sys, os, argparse, time, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.special import ellipk, eval_legendre, i0e, i1e, k0e, k1e
from scipy.optimize import minimize_scalar
from scipy.integrate import cumulative_trapezoid as ctz
import gate2c_robustness_20260925 as g2
from gate2c_robustness_20260925 import load, boot, G, a0k, kpc
from pb_solve_canon_20260925 import table_col

YD, YB = 0.5, 0.7
NMU = 16
_x, _w = np.polynomial.legendre.leggauss(2*NMU)
MU, WMU = _x[NMU:], _w[NMU:]                 # positive half-range nodes; weights sum to 1
MU_ALL = np.concatenate([MU, [0.0, 1.0]])    # + probes: midplane (index NMU), axis (NMU+1)
LMAX_FULL = 12
NR, NRC = 400, 120                           # fine shooting grid, coarse ring-sum grid

def rotmod(g):
    return np.atleast_2d(np.loadtxt(os.path.join(g2.SPARC_DIR, g+'_rotmod.dat'), comments='#'))

def freeman_v2(R, M, Rs):
    y = np.clip(R/(2*Rs), 1e-6, None)
    return 2*G*M/Rs * y**2 * (i0e(y)*k0e(y) - i1e(y)*k1e(y))

def ring_sources(Rs, sig, h, nz=8):
    e = np.sqrt(Rs[1:]*Rs[:-1]); e = np.concatenate([[0], e, [Rs[-1]*np.sqrt(Rs[-1]/Rs[-2])]])
    dm = sig*np.pi*(e[1:]**2 - e[:-1]**2)
    zk = -h*np.log(1 - (np.arange(nz)+0.5)/nz)
    return (np.repeat(Rs, 2*nz), np.tile(np.concatenate([zk, -zk]), len(Rs)),
            np.repeat(dm/(2*nz), 2*nz), dm.sum())

def ring_phi(Rf, zf, Rs, zs, ms, eps, chunk=300):
    out = np.empty(len(Rf))
    for a in range(0, len(Rf), chunk):
        Rq = Rf[a:a+chunk, None]; zq = zf[a:a+chunk, None]
        d2 = (Rq+Rs)**2 + (zq-zs)**2 + eps**2
        out[a:a+chunk] = -G*np.sum(ms*2*ellipk(4*Rq*Rs/d2)/(np.pi*np.sqrt(d2)), axis=1)
    return out

def prep(d, raw, mode, Rd, L, MHI, R_out):
    r = d.r.values
    rr = np.concatenate([[0], np.geomspace(0.005, R_out, NR-1)])
    diag = dict(mratio=np.nan, vbratio=np.nan, Rg=np.nan)
    if mode == 'C0':
        vb2 = np.clip(d.vbar2.values, 0, None)
        v2 = np.where(rr <= r[-1], np.interp(rr, np.r_[0, r], np.r_[0, vb2]), vb2[-1]*r[-1]/np.maximum(rr, 1e-9))
        integ = np.where(rr > 0, v2/np.maximum(rr, 1e-12), 0); integ[0] = integ[1]
        psi = ctz(integ, rr, initial=0)
        return rr, np.repeat(psi[:, None], len(MU_ALL), axis=1), diag
    R0, sb, vg, vbul = raw[:, 0], raw[:, 6], raw[:, 3], raw[:, 5]
    hz = 0.196*Rd**0.633
    Rs = np.geomspace(0.005*Rd, max(10*Rd, 1.5*R0[-1]), 160)
    sbq = np.interp(Rs, R0, np.clip(sb, 0, None)); pos = np.where(sb > 0)[0]
    if len(pos):
        jl = pos[-1]; b = Rs > R0[jl]; sbq[b] = sb[jl]*np.exp(-(Rs[b]-R0[jl])/Rd)
    sbq[Rs < R0[0]] = max(sb[0], 0)
    Rd_, zd_, md_, Mdisk = ring_sources(Rs, YD*sbq*1e6, hz)
    Mg = 1.33*MHI*1e9; ok = vg > 0
    if ok.sum() >= 3 and Mg > 0:
        Rg = np.exp(minimize_scalar(lambda lr: np.sum((np.sqrt(freeman_v2(R0[ok], Mg, np.exp(lr)))-vg[ok])**2),
                                    bounds=(np.log(0.2*Rd), np.log(30*Rd)), method='bounded').x)
    else: Rg = 2*Rd
    Rgs = np.geomspace(0.005*Rg, max(10*Rg, 1.5*R0[-1]), 160)
    Rg_, zg_, mg_, _ = ring_sources(Rgs, Mg/(2*np.pi*Rg**2)*np.exp(-Rgs/Rg), 0.2)
    rc = np.concatenate([[0], np.geomspace(0.01, R_out, NRC)])
    Rf = (rc[:, None]*np.sqrt(1-MU_ALL**2)).ravel(); zf = (rc[:, None]*MU_ALL).ravel()
    phi = (ring_phi(Rf, zf, Rd_, zd_, md_, hz/2) + ring_phi(Rf, zf, Rg_, zg_, mg_, 0.1)).reshape(len(rc), len(MU_ALL))
    psic = phi - phi[0, 0]
    psib = np.column_stack([np.interp(rr, rc, psic[:, j]) for j in range(len(MU_ALL))])
    if np.any(vbul > 0):                              # spherical bulge, M(<r) = Yb vbul^2 r / G
        vb2b = YB*vbul**2
        v2 = np.where(rr <= R0[-1], np.interp(rr, np.r_[0, R0], np.r_[0, vb2b]), vb2b[-1]*R0[-1]/np.maximum(rr, 1e-9))
        integ = np.where(rr > 0, v2/np.maximum(rr, 1e-12), 0); integ[0] = integ[1]
        psib = psib + ctz(integ, rr, initial=0)[:, None]
    else: vb2b = 0*vbul
    v2m = rc*np.gradient(phi[:, NMU], rc)             # disk+gas in-plane v^2 (model)
    spc = YD*raw[:, 4]**2 + np.sign(vg)*vg**2
    m = np.interp(R0, rc, v2m); sel = spc > 0
    diag = dict(mratio=Mdisk/(YD*L*1e9) if L > 0 else np.nan,
                vbratio=float(np.median(m[sel]/spc[sel])) if sel.any() else np.nan, Rg=Rg)
    return rr, psib, diag

def multipoles(rr, rho, ls):
    phis, dphis = [], []
    for l in ls:
        rl = (2*l+1)*(rho @ (WMU*eval_legendre(l, MU)))
        Iin = ctz(rl*rr**(l+2), rr, initial=0)
        b = np.zeros_like(rr); b[1:] = rl[1:]*rr[1:]**(1.0-l)
        Iout = -ctz(b[::-1], rr[::-1], initial=0)[::-1]        # int_r^Rout, accumulated outside-in
        c = -4*np.pi*G/(2*l+1); ph = np.zeros_like(rr); dph = np.zeros_like(rr); s = rr > 0
        ph[s] = c*(Iin[s]*rr[s]**(-l-1.0) + rr[s]**l*Iout[s])
        dph[s] = c*(-(l+1)*Iin[s]*rr[s]**(-l-2.0) + l*rr[s]**(l-1.0)*Iout[s])
        phis.append(ph); dphis.append(dph)
    return np.array(phis), np.array(dphis)

def solve_halo(rr, psib, s2, lrc, amb, lmax, Phil=None, iters=25, damp=0.5):
    rc = 10**lrc; nr = len(rr); ls = np.arange(2, lmax+1, 2)
    Pl = np.array([eval_legendre(l, MU_ALL) for l in ls]) if len(ls) else None
    if Phil is None or not len(ls): Phil = np.zeros((len(ls), nr))
    dPhil = np.zeros_like(Phil); it = 0
    for it in range(iters if len(ls) else 1):
        base = psib + (Phil.T @ Pl if len(ls) else 0)
        bm = base[:, :NMU]
        f = lambda i, p: 4*np.pi*rr[i]**2*(rc*np.dot(WMU, np.exp(np.clip(-(bm[i]+p)/s2, -700, 50))) - amb)
        M = np.zeros(nr); p0 = np.zeros(nr)
        for i in range(nr-1):
            h = rr[i+1]-rr[i]; fi = f(i, p0[i]); gi = G*M[i]/rr[i]**2 if rr[i] > 0 else 0.0
            Mp = M[i]+h*fi; pp = p0[i]+h*gi
            M[i+1] = M[i]+0.5*h*(fi+f(i+1, pp)); p0[i+1] = p0[i]+0.5*h*(gi+G*Mp/rr[i+1]**2)
        psi = base + p0[:, None]
        if not len(ls): break
        new, dnew = multipoles(rr, rc*np.exp(np.clip(-psi[:, :NMU]/s2, -700, 50)), ls)
        delta = np.max(np.abs(new-Phil)); dPhil = dnew
        if delta < 1e-3*s2: Phil = new; break
        Phil = Phil + damp*(new-Phil)
    return M, p0, Phil, dPhil, psi, it

def fit_gal(args):
    g, d, mode, smode, amb = args
    t0 = time.time()
    r = d.r.values; vobs = d.vobs.values; Rd = d.Rd.iloc[0]; L = d.L.iloc[0]; MHI = d.MHI.iloc[0]
    Mb = float(0.5*L*1e9 + 1.33*MHI*1e9)
    s2 = np.sqrt(a0k*G*Mb)/2 if smode == 'ideal' else d.Vf.iloc[0]**2/2
    R_out = 2.0*r[-1]
    rr, psib, diag = prep(d, rotmod(g), mode, Rd, L, MHI, R_out)
    lmax = LMAX_FULL if mode == 'full' else 0
    ls = np.arange(2, lmax+1, 2); P0 = eval_legendre(ls, 0.0) if len(ls) else np.zeros(0)
    st = {'Phil': None}
    def model(lr):
        M, p0, Phil, dPhil, psi, it = solve_halo(rr, psib, s2, lr, amb, lmax, st['Phil'])
        st['Phil'] = Phil; st['last'] = (M, p0, Phil, dPhil, psi, it)
        vnb2 = G*np.interp(r, rr, M)/r
        if len(ls): vnb2 = vnb2 + r*np.interp(r, rr, (dPhil*P0[:, None]).sum(0))
        return d.vbar2.values + vnb2
    w = 1/(0.1*vobs)**2
    chi = lambda lr: np.sum(w*(np.sqrt(np.clip(model(lr), 1, None))-vobs)**2)
    lr = minimize_scalar(chi, bounds=(4, 11), method='bounded', options={'xatol': 0.02}).x
    v2 = model(lr); M, p0, Phil, dPhil, psi, it = st['last']
    r0 = np.sqrt(9*s2/(4*np.pi*G*10**lr))
    out = dict(gal=g, mode=mode, lrhoc=lr, r0=r0, Rd=Rd, Mb=Mb, s2=s2, chi2=chi(lr)/len(r), iters=it,
               secs=time.time()-t0, **diag)
    for c in (2.0, 3.2, 5.0):                        # isodensity (= isopotential) axis ratios
        R0 = c*Rd
        for tag, ps in (('q', psi), ('qb', psib)):
            if R0 < rr[-1] and mode != 'C0':
                tgt = np.interp(R0, rr, ps[:, NMU]); ax = np.maximum.accumulate(ps[:, NMU+1])
                out[f'{tag}{c:g}'] = np.interp(tgt, ax, rr)/R0 if tgt <= ax[-1] else np.nan
            else: out[f'{tag}{c:g}'] = np.nan
    dd = d.copy(); dd['lgPB'] = np.log10(np.clip(v2, 1e-3, None)/r*1e6/kpc); dd['mode'] = mode
    return out, dd

def summarize(Pp, C, mode):
    Pp = Pp.reset_index(drop=True); Pp['resPB'] = Pp.lgo - Pp.lgPB
    clean = (Pp.Q == 1) & (Pp.inc <= 80) & (Pp.rdr >= 4)
    print(f'\n== {mode}: N {len(Pp)} gal {Pp.gal.nunique()} scatter {Pp.resPB.std():.3f} mean {Pp.resPB.mean():+.3f}')
    for k, m in {'ALL': Pp.Q > 0, 'Q1 clean': clean}.items():
        g2.rng = np.random.default_rng(1)
        e, lo, hi, n, ng = boot(Pp[m].reset_index(drop=True), ['lrRd'], y='resPB', nb=200)
        print(f'  {k:9s} r/Rd slope {e[0]:+.3f} [{lo[0]:+.3f},{hi[0]:+.3f}]  (N={n}, gal={ng})')
    ok = (C.lrhoc > 4.05) & (C.lrhoc < 10.95)
    x = np.log10(C.Mb[ok].values); y = np.log10(10**C.lrhoc[ok].values*C.r0[ok].values/1e6)
    rng = np.random.default_rng(1); bs = []
    for _ in range(1000):
        i = rng.integers(0, len(x), len(x)); bs.append(np.polyfit(x[i], y[i], 1)[0])
    print(f'  F23 log(rho_c r0 /Msun pc^-2) vs log Mb: slope {np.polyfit(x, y, 1)[0]:+.3f} '
          f'[{np.percentile(bs, 2.5):+.3f},{np.percentile(bs, 97.5):+.3f}]  median {np.median(y):.2f}  (gal {ok.sum()}, at bound {(~ok).sum()})')
    print(f'  median r0/Rd {np.median(C.r0/C.Rd):.2f}  median sec/gal {C.secs.median():.1f}  max iters {C.iters.max()}')
    if mode != 'C0':
        print('  baryon model checks: median disk-mass ratio %.2f ; median in-plane v^2(disk+gas) model/SPARC %.2f [IQR %.2f-%.2f]'
              % (C.mratio.median(), C.vbratio.median(), C.vbratio.quantile(.25), C.vbratio.quantile(.75)))
        print('  isodensity q (nb) at 2/3.2/5 Rd: ' + ' '.join(f'{C[f"q{c:g}"].median():.2f}' for c in (2, 3.2, 5)) +
              ' ; baryon-only isopotential qb: ' + ' '.join(f'{C[f"qb{c:g}"].median():.2f}' for c in (2, 3.2, 5)))

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('modes'); ap.add_argument('smode', nargs='?', default='ideal')
    ap.add_argument('--workers', type=int, default=1); ap.add_argument('--gal', default=None)
    ap.add_argument('--amb', type=float, default=1e3); a = ap.parse_args()
    P = load(); P['Vf'] = P.gal.map(table_col(15))
    if a.smode != 'ideal': P = P[P.Vf > 0].reset_index(drop=True)
    if a.gal: P = P[P.gal.isin(a.gal.split(','))].reset_index(drop=True)
    print('Gate 2b | smode', a.smode, '| rho_amb', a.amb, '| N', len(P), 'gal', P.gal.nunique(), '| workers', a.workers)
    for mode in a.modes.split(','):
        jobs = [(g, d, mode, a.smode, a.amb) for g, d in P.groupby('gal')]
        with Pool(a.workers) as pool: res = pool.map(fit_gal, jobs, chunksize=1)
        C = pd.DataFrame([r[0] for r in res]); Pp = pd.concat([r[1] for r in res])
        tag = f'gate2b_{mode}_{a.smode}' + ('_test' if a.gal else '')
        C.to_csv(tag+'_gals.csv', index=False); Pp.to_csv(tag+'_points.csv', index=False)
        if a.gal: print(C[['gal', 'lrhoc', 'r0', 'Rd', 'chi2', 'iters', 'secs', 'mratio', 'vbratio', 'q2', 'qb2']].to_string())
        else: summarize(Pp, C, mode)
