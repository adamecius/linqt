#!/usr/bin/env python3
"""Stage-5 tolerance table (handoff/0005 codex.md section 2), called by scripts/check_goldens.sh 5.

usage: figures/stage_5_table.py <run_dir> <goldens/stage0/bases>

<run_dir> holds the stage5_bases outputs of one thread count laid out by the gate as
<model>/W0/{real,bloch,blochx,eigen}/{dos,kuboxx}/ and <model>/W1/{real,bloch}/{dos,kuboxx}/
(real, bloch: random phase, 32 states, seed 12345; blochx, eigen: exact trace).

Prints one line per check: "PASS" or "FAIL" for a gated row (counted by the gate),
"INFO" for a reported-only row. Definitions:

  nrel(a, b) = max_i |a_i - b_i| / max_i |b_i|      (sup-norm relative; the gated quantity)
  prel(a, b) = max_i |a_i - b_i| / max(|a_i|, |b_i|) (pointwise, as scripts/golden_diff.py)

Moments are read with lsquant._core.ReadChebmom. Curves are recomputed at full double
precision from the moment files with lsquant.reconstruct.dos_curve / kubo_greenwood at
the run's broadening (the written .dat curves carry 6 significant digits, too few for the
1e-8 row); the written .dat is checked against the recomputed curve (INFO).

G9 (hr_exactdiag.py, 48 x 48 shifted mesh (i + 1/2)/48, Gaussian sigma = pi a / 256) is
interpolated linearly onto the KPM energies. Units: DOS_kpm / N_k (hr_exactdiag DOS is per
cell); sigma_kpm * 2 pi / (A_hr N_k) with A_hr = |a1_x a2_y| as hr_exactdiag cmd_sigma uses it
(the legacy Kubo-Greenwood curve is pi Tr[v delta(E-H) v delta(E-H)], hr_exactdiag's is
(2 pi^2 / A_hr) (1/N_k) sum_k sum_nm |v_nm|^2 g(E-e_n) g(E-e_m)); the least-squares factor is
printed next to the analytic one. Energy masks: graphene |E| < 2.5, Haldane 0.5 < |E| < 2.5.
"""
import glob
import json
import math
import os
import sys

import numpy as np

from lsquant import _core
from lsquant.reconstruct import dos_curve, kubo_greenwood

NK = 48 * 48
M = 256
TOL_MOM, TOL_BE, TOL_DOS, TOL_SIG = 1e-9, 1e-8, 2e-2, 5e-2
MASK = {"graphene": lambda E: np.abs(E) < 2.5,
        "haldane": lambda E: (np.abs(E) > 0.5) & (np.abs(E) < 2.5)}


def nrel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


def prel(a, b):
    d = np.abs(a - b)
    s = np.maximum(np.abs(a), np.abs(b))
    r = np.where(s > 0, d / np.where(s > 0, s, 1), 0.0)
    return float(np.max(r))


def params(gold, m):
    p = {}
    for line in open(os.path.join(gold, m, "G9_params.txt")):
        k, v = line.split()
        p[k] = float(v)
    return p


def broadening_meV(p, W):
    """The run's broadening (check_goldens.sh stage5): pi * a / M with the W-widened window."""
    if W == 0:
        return float("%.7f" % (p["sigma_eV"] * 1000))
    a = (p["half_width"] / 1.10 + 0.5 * W) * 1.10
    return float("%.7f" % (math.pi * a / M * 1000))


def load(d, word, bmev):
    """(moments array, energies, full-precision curve, written .dat) of one run directory."""
    mom = glob.glob(os.path.join(d, "*.chebmom"))
    dat = glob.glob(os.path.join(d, "*.dat"))
    if len(mom) != 1 or len(dat) != 1:
        return None
    mu = _core.ReadChebmom(mom[0])
    E, y = (dos_curve if word == "dos" else kubo_greenwood)(mu, bmev)
    return np.asarray(mu.data), np.asarray(E), np.asarray(y), np.loadtxt(dat[0])


def row(ok, text):
    print(("PASS  " if ok else "FAIL  ") + text)


def info(text):
    print("INFO  " + text)


def g9(gold, m, word):
    f = os.path.join(gold, m, "G9_dos.json" if word == "dos" else "G9_sigma_xx.json")
    d = json.load(open(f))
    return np.array(d["energy_eV"]), np.array(d["dos" if word == "dos" else "sigma_xx_e2h"])


def area_hr(gold, m):
    lat = np.array([[float(x) for x in l.split()] for l in open(os.path.join(gold, m, "model", m + ".uc")) if l.split()][:3])
    return abs(lat[0, 0] * lat[1, 1]), abs(lat[0, 0] * lat[1, 1] - lat[0, 1] * lat[1, 0])


def gamma_gauss(gold, m, word, E, sig):
    """hr_exactdiag's Gaussian formula (dos / sigma --comp xx, same prefactors) evaluated on the
    commensurate Gamma mesh k = n/48 instead of its shifted mesh, at the energies E; plus the
    van Hove energies |E| at the three M points. Uses wannier2sparse/tools/hr_exactdiag.py."""
    sys.path.insert(0, os.environ.get("W2S_TOOLS", os.path.expanduser("~/wannier2sparse/tools")))
    from hr_exactdiag import read_hr, opk
    d = os.path.join(gold, m, "model")
    iR, H = read_hr(os.path.join(d, m + "_hr.dat"))
    n = np.arange(48) / 48.0
    K = np.stack(np.meshgrid(n, n, indexing="ij"), -1).reshape(-1, 2)
    kf = np.c_[K, np.zeros(len(K))]
    Hk = opk(iR, H, kf); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
    w, U = np.linalg.eigh(Hk)
    g = np.exp(-((E[None, None, :] - w[:, :, None]) ** 2) / (2 * sig ** 2)) / (math.sqrt(2 * math.pi) * sig)
    if word == "dos":
        out = g.sum((0, 1)) / len(kf)
    else:
        iRa, va = read_hr(os.path.join(d, m + "_vx_hr.dat"))
        Va = np.einsum("kmi,kmn,knj->kij", np.conj(U), opk(iRa, va, kf), U)
        A_hr = area_hr(gold, m)[0]
        out = np.einsum("knm,kne,kme->e", np.abs(Va) ** 2, g, g) * (2 * math.pi ** 2 / A_hr) / len(kf)
    hm = opk(iR, H, np.array([[0.5, 0, 0], [0, 0.5, 0], [0.5, 0.5, 0]]))
    evh = sorted(set(round(abs(x), 3) for x in np.linalg.eigvalsh(0.5 * (hm + np.conj(np.transpose(hm, (0, 2, 1))))).ravel()))
    return out, evh


def oracle_moments(gold, m, word, lo, hi):
    """Exact Chebyshev moments of the 48 x 48 Gamma-mesh system from hr_exactdiag's H(k) and
    the wannier2sparse velocity (independent of lsquant's model code), in the layout of the
    .chebmom files (kpm_weighted 1: mu_m (2 - delta_m0), mu_mn (2 - delta_m0)(2 - delta_n0),
    normalised by the dimension N_k W)."""
    sys.path.insert(0, os.environ.get("W2S_TOOLS", os.path.expanduser("~/wannier2sparse/tools")))
    from hr_exactdiag import read_hr, opk
    d = os.path.join(gold, m, "model")
    iR, H = read_hr(os.path.join(d, m + "_hr.dat"))
    n = np.arange(48) / 48.0
    K = np.stack(np.meshgrid(n, n, indexing="ij"), -1).reshape(-1, 2)
    kf = np.c_[K, np.zeros(len(K))]
    Hk = opk(iR, H, kf); Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
    w, U = np.linalg.eigh(Hk)
    x = (w - 0.5 * (lo + hi)) / (0.5 * (hi - lo))
    T = np.cos(np.arange(M)[None, None, :] * np.arccos(np.clip(x, -1, 1))[:, :, None])
    h = np.where(np.arange(M) == 0, 1.0, 2.0)
    if word == "dos":
        return T.sum((0, 1)) / x.size * h
    iRa, va = read_hr(os.path.join(d, m + "_vx_hr.dat"))
    V = np.einsum("kmi,kmn,knj->kij", np.conj(U), opk(iRa, va, kf), U)
    return np.einsum("kam,kab,kbn->mn", T, np.abs(V) ** 2, T) / x.size * np.outer(h, h)


def main():
    run, gold = sys.argv[1], sys.argv[2]
    for m in ("graphene", "haldane"):
        p = params(gold, m)
        A_hr, A_true = area_hr(gold, m)
        for word in ("dos", "kuboxx"):
            tag = "%s %s" % (m, "DOS" if word == "dos" else "sigma_xx")
            b0, b1 = broadening_meV(p, 0), broadening_meV(p, 1)
            R = {b: load(os.path.join(run, m, "W0", b, word), word, b0) for b in ("real", "bloch", "blochx", "eigen")}
            D = {b: load(os.path.join(run, m, "W1", b, word), word, b1) for b in ("real", "bloch")}
            missing = [b for b, v in R.items() if v is None] + ["W1 " + b for b, v in D.items() if v is None]
            if missing:
                row(False, "%s: outputs missing for %s" % (tag, ", ".join(missing)))
                continue
            # written .dat vs the full-precision recomputation (6 significant digits expected)
            for b, v in list(R.items()) + [("W1 " + b, v) for b, v in D.items()]:
                mu, E, y, dat = v
                if dat.shape[0] != E.size or np.max(np.abs(dat[:, 0] - E)) > 1e-5 * np.max(np.abs(E)):
                    row(False, "%s %s: written .dat grid differs from the recomputed curve grid" % (tag, b))
                    continue
                info("%s %-8s written .dat vs recomputed curve: nrel %.2e (6-digit text)" % (tag, b, nrel(dat[:, 1], y)))
            # W = 0: real vs Bloch moments (same random-phase vectors through B)
            r = nrel(R["real"][0], R["bloch"][0])
            row(r <= TOL_MOM, "%s W=0 real vs bloch moments (random 32, seed 12345): nrel %.3e (tol %.0e); prel %.3e"
                % (tag, r, TOL_MOM, prel(R["real"][0], R["bloch"][0])))
            info("%s W=0 real vs bloch curve: nrel %.3e" % (tag, nrel(R["real"][2], R["bloch"][2])))
            # Bloch vs eigen, same kernel, exact trace on both sides
            r = nrel(R["blochx"][2], R["eigen"][2])
            row(r <= TOL_BE, "%s W=0 bloch(exact) vs eigen curve: nrel %.3e (tol %.0e); moments nrel %.3e"
                % (tag, r, TOL_BE, nrel(R["blochx"][0], R["eigen"][0])))
            info("%s W=0 bloch(random 32) vs eigen curve: nrel %.3e (stochastic trace error)" % (tag, nrel(R["bloch"][2], R["eigen"][2])))
            info("%s W=0 real(random 32) vs eigen curve: nrel %.3e (stochastic trace error)" % (tag, nrel(R["real"][2], R["eigen"][2])))
            # any basis vs G9
            Eg, yg = g9(gold, m, word)
            E = R["eigen"][1]
            msk = MASK[m](E)
            ref = np.interp(E, Eg, yg)
            if word == "dos":
                fac, tol = 1.0 / NK, TOL_DOS
            else:
                fac, tol = 2 * math.pi / (A_hr * NK), TOL_SIG
                k = R["eigen"][2][msk]
                fit = float(np.dot(ref[msk], k) / np.dot(k, k))
                info("%s units: analytic factor 2 pi/(A_hr N_k) = %.6e (A_hr = %.6f, true cell area %.6f), least squares %.6e, ratio %.5f"
                     % (tag, fac, A_hr, A_true, fit, fit / fac))
            for b in ("real", "bloch", "blochx", "eigen"):
                y = R[b][2] * fac
                r = nrel(y[msk], ref[msk])
                i = int(np.argmax(np.abs(y[msk] - ref[msk])))
                row(r <= tol, "%s W=0 %-6s vs G9 on %s: nrel %.3e at E=%+.3f (tol %.0e); prel %.3e"
                    % (tag, b, "|E|<2.5" if m == "graphene" else "0.5<|E|<2.5", r, E[msk][i], tol, prel(y[msk], ref[msk])))
            info("%s W=0 eigen vs G9 over the whole KPM grid |E| <= %.3f: nrel %.3e"
                 % (tag, np.max(np.abs(E)), nrel(R["eigen"][2] * fac, ref)))
            # localisation of the G9 row: the exact moments of the same finite system from the
            # oracle's own H(k) and velocity, then mesh (shifted vs Gamma) and kernel (Jackson vs Gaussian)
            mo = oracle_moments(gold, m, word, p["window_lower"], p["window_upper"])
            info("%s W=0 independent moments (hr_exactdiag H(k) + wannier2sparse VX, Gamma mesh n/48, exact): "
                 "eigen nrel %.3e, bloch(exact) nrel %.3e, real(random 32) nrel %.3e, bloch(random 32) nrel %.3e"
                 % (tag, nrel(R["eigen"][0].reshape(mo.shape), mo), nrel(R["blochx"][0].reshape(mo.shape), mo),
                    nrel(R["real"][0].reshape(mo.shape), mo), nrel(R["bloch"][0].reshape(mo.shape), mo)))
            gam, evh = gamma_gauss(gold, m, word, E, p["sigma_eV"])
            info("%s W=0 mesh effect: G9 (mesh (i+1/2)/48) vs the same Gaussian formula on the Gamma mesh n/48: nrel %.3e on the mask"
                 % (tag, nrel(ref[msk], gam[msk])))
            info("%s W=0 kernel effect: eigen (Jackson, M=256) vs Gaussian on the Gamma mesh n/48: nrel %.3e on the mask"
                 % (tag, nrel(R["eigen"][2][msk] * fac, gam[msk])))
            far = msk.copy()
            for e in evh:
                far &= np.abs(np.abs(E) - e) > 0.15
            j = int(np.argmax(np.abs(R["eigen"][2][far] * fac - ref[far])))
            info("%s W=0 eigen vs G9 on the mask without |(|E| - E_vH)| <= 0.15 eV (E_vH = %s, M-point energies): nrel %.3e "
                 "at E=%+.3f; real %.3e, bloch %.3e"
                 % (tag, ", ".join("%.3f" % e for e in evh), nrel(R["eigen"][2][far] * fac, ref[far]), E[far][j],
                    nrel(R["real"][2][far] * fac, ref[far]), nrel(R["bloch"][2][far] * fac, ref[far])))
            # W = 1: real vs Bloch with the disorder sandwich B^+ V B
            r = nrel(D["real"][0], D["bloch"][0])
            row(r <= TOL_MOM, "%s W=1 (seed 7) real vs bloch moments (random 32, seed 12345): nrel %.3e (tol %.0e); prel %.3e"
                % (tag, r, TOL_MOM, prel(D["real"][0], D["bloch"][0])))
            info("%s W=1 real vs bloch curve: nrel %.3e" % (tag, nrel(D["real"][2], D["bloch"][2])))
            info("%s W=1 eigen: skipped as codex.md 2.4 directs (no per-k eigenbasis with disorder; N_k W = %d)"
                 % (tag, NK * 2))


if __name__ == "__main__":
    main()
