#!/usr/bin/env python3
# wb_conserve_r1pc_20261015.py = wb_conserve_20261011.py with RS extended to 24, 32 r_M (~1 pc) + timing. 2026-10-15. Physics unchanged. Result: RESULT_WB_R1PC_20261015.md
import sys
import numpy as np

GE = 1.2
RHO_HOST = 0.01 / 2.52e4
RS = (1.0, 2.0, 4.0, 8.0, 16.0, 24.0, 32.0)


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
    if R is None:
        return F["full"]
    if R == 0:
        return F["full"] + F["def"]
    f = M[f"defR{R}"] / M[f"posR{R}"]
    return F["full"] + F["def"] - f * F[f"posR{R}"], f


if __name__ == "__main__":
    import time
    t0 = time.time()
    s, th, W, N, b, Nf = (float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]),
                          float(sys.argv[5]), int(sys.argv[6]))
    n, x1, x2, fine = setup(s, th, b)
    F, M = run(((0.5, x1), (0.5, x2)), [x1, x2], W, N, fine, Nf)
    aN = np.array([-0.5 * n / s ** 2, 0.5 * n / s ** 2])
    out = [f"s={s:.2f} th={th:.1f} W={W:.0f} N={N} b={b} Nf={Nf}"]
    for R in (None, 0) + RS:
        r = model_force(F, M, R)
        Fm, f = (r if isinstance(r, tuple) else (r, None))
        a = aN + Fm
        g = -((a[0] - a[1]) @ n) * s ** 2
        out.append(f"{'MOND' if R is None else 'R'+str(R)}={g:.4f}" + (f"(f={f:.2f})" if f is not None else ""))
    out.append(f"Mdef={M['def']:.3f} t={time.time()-t0:.0f}s")
    print(" ".join(out), flush=True)
