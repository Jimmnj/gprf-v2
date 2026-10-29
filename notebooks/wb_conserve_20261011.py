#!/usr/bin/env python3
"""
wb_conserve_20261011.py -- conserving floor-clip family for the wide-binary boost (GPRF-H). 2026-10-11 [UNRATIFIED]
Builds on wb_clip_gamma_20261011.py (same physics, same validated rho-kernel nested numerics).
Result: RESULT_WB_CLIP_GAMMA_20261011.md sec 8.

Conservation family (GPRF-H has no transport law, D6 -> bracket, not solve):
  MOND phantom rho_ph. Floor deficit d = max(0, -RHO_HOST - rho_ph) (medium the collar cannot give up).
  Net inflow from far is kept at MOND's value. Inside a sphere of radius R about the binary COM, the deficit is removed
  from the positive phantom in the same sphere, uniformly in proportion:
     rho_nb = rho_ph + d - f_R * rho_ph^+ * [r < R],    f_R = D_R / P_R,  D_R = int_{r<R} d,  P_R = int_{r<R} rho_ph^+
  Outside R the deficit is filled from far (one-shot). R -> 0: one-shot clip (upper bound). R large: full conservation.
  VALID ONLY for R > s/2 + ~1 (sphere must contain both stars and their cocoons); smaller R is ill-posed.
Forces are linear in density -> one pass gives every R.
Usage: python wb_conserve_20261011.py single
       python wb_conserve_20261011.py binary s theta W N b Nf   (sweep: (1 6 192 0 192) (2 8 192 0 192) (4 16 192 1.6 128))
"""
import sys
import numpy as np

GE = 1.2
RHO_HOST = 0.01 / 2.52e4
RS = (1.0, 2.0, 4.0, 8.0, 16.0)


def nu_prime(y):
    u = np.sqrt(y)
    e = np.exp(-u)
    return -e / ((1.0 - e) ** 2 * 2.0 * u)


def fields(P, masses, ge, h):
    g = np.zeros_like(P)
    g[2] -= ge
    T = np.zeros((3, 3) + P.shape[1:])
    for m, xs in masses:
        d = P - xs[:, None, None, None]
        r2 = np.maximum((d ** 2).sum(0), (0.25 * h) ** 2)
        r = np.sqrt(r2)
        r3 = r2 * r
        r5 = r3 * r2
        g -= m * d / r3
        for a in range(3):
            for b in range(3):
                T[a, b] -= m * ((1.0 if a == b else 0.0) / r3 - 3.0 * d[a] * d[b] / r5)
    return g, T, np.sqrt((g ** 2).sum(0))


def channels(masses, evalpts, centre, half, N, exclude=(), ge=GE):
    """Forces at evalpts and masses for channels: full, deficit, and per R: deficit_inR, positive_inR."""
    h = 2 * half / N
    c = -half + h * (np.arange(N) + 0.5)
    dV = h ** 3
    names = ["full", "def"] + [f"defR{R}" for R in RS] + [f"posR{R}" for R in RS]
    F = {k: np.zeros((len(evalpts), 3)) for k in names}
    M = {k: 0.0 for k in names}
    for k0 in range(0, N, 16):
        nz = min(16, N - k0)
        gx, gy, gz = np.meshgrid(c + centre[0], c + centre[1], c[k0:k0 + nz] + centre[2], indexing="ij")
        P = np.stack([gx, gy, gz])
        g, T, y = fields(P, masses, ge, h)
        gTg = np.einsum("aijk,abijk,bijk->ijk", g, T, g)
        rho = np.nan_to_num(-nu_prime(y) * gTg / (4 * np.pi * y))
        keep = np.ones(rho.shape, bool)
        for (ctr, b) in exclude:
            keep &= ~(np.abs(P - ctr[:, None, None, None]) < b).all(0)
        rho = np.where(keep, rho, 0.0)
        dfc = np.maximum(0.0, -RHO_HOST - rho) * keep
        pos = np.maximum(rho, 0.0)
        rr = np.sqrt((P ** 2).sum(0))
        dens = {"full": rho, "def": dfc}
        for R in RS:
            inR = rr < R
            dens[f"defR{R}"] = dfc * inR
            dens[f"posR{R}"] = pos * inR
        for k, d in dens.items():
            M[k] += d.sum() * dV
        for i, xe in enumerate(evalpts):
            Rv = xe[:, None, None, None] - P
            R2 = (Rv ** 2).sum(0)
            ok = R2 > (0.5 * h) ** 2
            K = (-Rv) / np.where(ok, R2 * np.sqrt(R2), 1.0) * ok * dV
            for k, d in dens.items():
                F[k][i] += (K * d).sum(axis=(1, 2, 3))
    return F, M


def setup(s, th, b):
    thr = np.radians(th)
    n = np.array([np.sin(thr), 0.0, np.cos(thr)])
    x1, x2 = 0.5 * s * n, -0.5 * s * n
    if s <= 2.0:
        bb = 0.5 * s + 1.2
        fine = [(np.zeros(3), bb)]
    else:
        fine = [(x1, b), (x2, b)]
    return n, x1, x2, fine


def run(masses, ev, W, N, fine, Nf):
    F, M = channels(masses, ev, np.zeros(3), W, N, exclude=fine)
    for ctr, bb in fine:
        Ff, Mf = channels(masses, ev, ctr, bb, Nf)
        for k in F:
            F[k] += Ff[k]
            M[k] += Mf[k]
    return F, M


def model_force(F, M, R):
    if R is None:      # MOND
        return F["full"]
    if R == 0:         # one-shot clip
        return F["full"] + F["def"]
    f = M[f"defR{R}"] / M[f"posR{R}"]
    return F["full"] + F["def"] - f * F[f"posR{R}"], f


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "binary":
        s, th, W, N, b, Nf = (float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]),
                              float(sys.argv[6]), int(sys.argv[7]))
        n, x1, x2, fine = setup(s, th, b)
        F, M = run(((0.5, x1), (0.5, x2)), [x1, x2], W, N, fine, Nf)
        aN = np.array([-0.5 * n / s ** 2, 0.5 * n / s ** 2])
        out = [f"{s:5.2f} {th:5.1f}"]
        for R in (None, 0) + RS:
            r = model_force(F, M, R)
            Fm, f = (r if isinstance(r, tuple) else (r, None))
            a = aN + Fm
            g = -((a[0] - a[1]) @ n) * s ** 2
            out.append(f"{'MOND' if R is None else 'R'+str(R)}={g:.4f}" + (f"(f={f:.2f})" if f is not None else ""))
        print(" ".join(out), flush=True)
    elif mode == "single":
        x0 = np.array([0.013, 0.007, 0.011])
        F, M = run(((1.0, x0),), [x0], 16.0, 192, [(x0, 1.5)], 160)
        out = ["single star self-accel z (a0):"]
        for R in (None, 0) + RS:
            r = model_force(F, M, R)
            Fm, f = (r if isinstance(r, tuple) else (r, None))
            out.append(f"{'MOND' if R is None else 'R'+str(R)}={Fm[0][2]:+.4f}" + (f"(f={f:.2f})" if f is not None else ""))
        print(" ".join(out), flush=True)
