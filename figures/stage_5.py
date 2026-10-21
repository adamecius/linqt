#!/usr/bin/env python3
"""Stage-5 figure (handoff 0005, codex.md section 3): graphene DOS and sigma_xx in three bases.

Reads the stage5_bases outputs kept by the gate (`S5_KEEP=<dir> scripts/check_goldens.sh 5`,
layout <dir>/T16/graphene/W{0,1}/<basis>/{dos,kuboxx}/) and the committed G9 goldens
(goldens/stage0/bases/graphene/G9_dos.json, G9_sigma_xx.json, G9_params.txt). Curves are
recomputed at full precision from the moment files with lsquant.reconstruct (dos_curve,
kubo_greenwood) at the run's broadening, exactly as figures/stage_5_table.py does.

Panels: (a) DOS, W = 0: real (random phase 32), Bloch (random phase 32), eigen (exact) as
lines, hr_exactdiag G9 as markers, residual basis - G9 below; (b) DOS, W = 1 eV (seed 7):
real and Bloch as lines, the clean G9 as markers for reference, residual real - Bloch below;
(c), (d) the same for sigma_xx. Units: DOS per cell (KPM / N_k); sigma_xx in e^2/h per
cell area with the true honeycomb cell A = |a1 x a2| (hr_exactdiag cmd_sigma divides by
|a1_x a2_y|, half the honeycomb cell, so its values are rescaled by that ratio here).

Outputs figures/out/stage_5/fig_stage_5.png and fig_stage_5.caption.txt.
Usage (repository root, on TCN):
    MPLBACKEND=Agg PYTHONPATH=<build>/lsquant_core:python python3 figures/stage_5.py --runs <dir>/T16 \
        [--table <gate table file>]
"""
import argparse
import glob
import json
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lsquant import _core
from lsquant.reconstruct import dos_curve, kubo_greenwood

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GOLD = os.path.join(REPO, "goldens", "stage0", "bases", "graphene")
NK, M = 48 * 48, 256
COL = {"real": "#0072B2", "bloch": "#E69F00", "eigen": "#009E73", "g9": "0.15"}
LBL = {"real": "real space (random 32)", "bloch": "Bloch (random 32)", "eigen": "eigen (exact)"}


def params():
    return {l.split()[0]: float(l.split()[1]) for l in open(os.path.join(GOLD, "G9_params.txt"))}


def bmev(p, W):
    a = p["half_width"] if W == 0 else (p["half_width"] / 1.10 + 0.5 * W) * 1.10
    return float("%.7f" % (p["sigma_eV"] * 1000)) if W == 0 else float("%.7f" % (math.pi * a / M * 1000))


def curve(d, word, b):
    mu = _core.ReadChebmom(glob.glob(os.path.join(d, "*.chebmom"))[0])
    E, y = (dos_curve if word == "dos" else kubo_greenwood)(mu, b)
    return np.asarray(E), np.asarray(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, help="<S5_KEEP>/T16")
    ap.add_argument("--table", default=None, help="gate table output (stage_5_table.py) for the caption")
    ap.add_argument("--out", default=os.path.join(REPO, "figures", "out", "stage_5"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    p = params()
    lat = np.array([[float(x) for x in l.split()] for l in open(os.path.join(GOLD, "model", "graphene.uc")) if l.split()][:3])
    A_hr = abs(lat[0, 0] * lat[1, 1])
    A = abs(lat[0, 0] * lat[1, 1] - lat[0, 1] * lat[1, 0])
    unit = {"dos": 1.0 / NK, "kuboxx": 2 * math.pi / (A * NK)}
    g9 = {"dos": json.load(open(os.path.join(GOLD, "G9_dos.json"))),
          "kuboxx": json.load(open(os.path.join(GOLD, "G9_sigma_xx.json")))}
    g9 = {"dos": (np.array(g9["dos"]["energy_eV"]), np.array(g9["dos"]["dos"])),
          "kuboxx": (np.array(g9["kuboxx"]["energy_eV"]), np.array(g9["kuboxx"]["sigma_xx_e2h"]) * A_hr / A)}

    plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "lines.linewidth": 1.0, "axes.unicode_minus": False})
    fig = plt.figure(figsize=(7.0, 6.4))
    gs = fig.add_gridspec(5, 2, height_ratios=[3, 1.1, 0.45, 3, 1.1], hspace=0.08, wspace=0.28)
    stats = {}
    letters = iter("abcd")
    for r, word in enumerate(("dos", "kuboxx")):
        for c, W in enumerate((0, 1)):
            ax = fig.add_subplot(gs[3 * r, c]); rx = fig.add_subplot(gs[3 * r + 1, c], sharex=ax)
            Eg, yg = g9[word]
            bases = ("real", "bloch", "eigen") if W == 0 else ("real", "bloch")
            C = {b: curve(os.path.join(a.runs, "graphene", "W%d" % W, b, word), word, bmev(p, W)) for b in bases}
            for b in bases:
                E, y = C[b]
                ax.plot(E, y * unit[word], color=COL[b], lw=1.6 if b == "real" else 1.0,
                        ls={"real": "-", "bloch": "--", "eigen": ":"}[b], label=LBL[b])
            step = max(1, len(Eg) // 330)
            ax.plot(Eg[::step], yg[::step], "o", ms=1.6, mfc="none", mec=COL["g9"], mew=0.4,
                    label="hr_exactdiag G9" + (" (clean)" if W else ""))
            E = C["real"][0]
            if W == 0:
                ref = np.interp(E, Eg, yg)
                for b in bases:
                    rx.plot(E, C[b][1] * unit[word] - ref, color=COL[b], lw=0.8,
                            ls={"real": "-", "bloch": "--", "eigen": ":"}[b])
                msk = np.abs(E) < 2.5
                stats[(word, W)] = {b: float(np.max(np.abs(C[b][1][msk] * unit[word] - ref[msk])) / np.max(np.abs(ref[msk]))) for b in bases}
                rx.set_ylabel("basis $-$ G9", fontsize=7)
            else:
                d = (C["real"][1] - C["bloch"][1]) * unit[word]
                rx.plot(E, d, color=COL["real"], lw=0.8)
                stats[(word, W)] = {"real-bloch": float(np.max(np.abs(d)) / np.max(np.abs(C["bloch"][1] * unit[word])))}
                rx.set_ylabel("real $-$ Bloch", fontsize=7)
            rx.axhline(0, color="0.6", lw=0.5)
            for x in (-2.5, 2.5):
                ax.axvline(x, color="0.85", lw=0.5, zorder=0)
            ax.set_xlim(-3.4, 3.4)
            plt.setp(ax.get_xticklabels(), visible=False)
            if r == 1:
                rx.set_xlabel("$E$ (eV)")
            ax.set_ylabel(r"DOS (eV$^{-1}$ cell$^{-1}$)" if word == "dos" else r"$\sigma_{xx}$ ($e^2/h$)")
            ax.set_title("(%s) %s, $W$ = %d eV" % (next(letters), "DOS" if word == "dos" else r"$\sigma_{xx}$", W),
                         loc="left", fontsize=8)
            ax.legend(loc="upper right", fontsize=5.5, frameon=False, handlelength=2.2)
            rx.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
            rx.yaxis.get_offset_text().set_fontsize(6)
    png = os.path.join(a.out, "fig_stage_5.png")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    tab = ""
    if a.table and os.path.isfile(a.table):
        txt = open(a.table).read()

        def num(pat):
            mm = re.search(pat, txt)
            return mm.group(1) if mm else "n/a"
        parts = []
        for w, lab in (("DOS", "DOS"), ("sigma_xx", "sigma_xx")):
            parts.append("%s: real vs Bloch moments %s (W = 0), %s (W = 1), tol 1e-9; Bloch (exact trace) vs eigen curve %s, tol 1e-8; "
                         "eigen moments vs exact moments from hr_exactdiag's H(k) on the same Gamma mesh %s"
                         % (lab,
                            num(r"graphene %s W=0 real vs bloch moments[^:]*: nrel (\S+)" % w),
                            num(r"graphene %s W=1 \(seed 7\) real vs bloch moments[^:]*: nrel (\S+)" % w),
                            num(r"graphene %s W=0 bloch\(exact\) vs eigen curve: nrel (\S+)" % w),
                            num(r"graphene %s W=0 independent moments[^:]*: eigen nrel ([0-9.e+-]+)," % w)))
        tab = (" Gate (figures/stage_5_table.py, T = 16; identical at T = 1), graphene: " + "; ".join(parts) +
               ". The G9 rows fail the codex tolerances (2e-2 DOS, 5e-2 sigma_xx): shifted vs Gamma mesh and "
               "Jackson vs Gaussian kernel (tests.md section 4.3).")
    s0d, s0s, s1d, s1s = stats[("dos", 0)], stats[("kuboxx", 0)], stats[("dos", 1)], stats[("kuboxx", 1)]
    caption = (
        "FIG. 5. Graphene (wannier2sparse examples/gen_models.py, t = -1 eV, nearest-neighbour distance 1 A, "
        "supercell 48 x 48 x 1 = 4608 orbitals, Bloch mesh k = n/48) evaluated by lsquant.apps.stage5_bases in three "
        "bases from one model file, M = 256 Chebyshev moments, Jackson kernel. Top: DOS per cell; bottom: sigma_xx "
        "(Kubo-Greenwood) in e^2/h with the true cell area %.4f A^2. Left column W = 0: real space (CSR, random phase, "
        "32 states, seed 12345, solid), Bloch (block-diagonal H(k), the same 32 vectors mapped by B, dashed), eigen "
        "(exact trace in the eigenbasis of each H(k), dotted), window [%.2f, %.2f] eV from the eigenvalues padded "
        "10 %%, broadening %.4f meV (= pi a / M, the Jackson resolution of M = 256, a = %.2f eV); markers: G9, "
        "tools/hr_exactdiag.py (wannier2sparse 3f553db) dos and sigma --comp xx on the 48 x 48 mesh (i + 1/2)/48 with "
        "Gaussian sigma = %.4f meV; residual basis - G9 below. Right column W = 1 eV Anderson disorder (uniform in "
        "[-0.5, 0.5] eV, numpy default_rng seed 7): real space (CSR + diag V) and Bloch (H(k) + B V B^+ through FFTW), "
        "window [-3.85, 3.85] eV, broadening %.4f meV; clean G9 markers for reference; residual real - Bloch below "
        "(eigen skipped: no per-k eigenbasis with disorder). Max |basis - G9| / max |G9| on |E| < 2.5 eV, W = 0: DOS "
        "real %.2e, Bloch %.2e, eigen %.2e; sigma_xx real %.2e, Bloch %.2e, eigen %.2e. W = 1: max |real - Bloch| / "
        "max |Bloch| on the curve: DOS %.1e, sigma_xx %.1e.%s"
        % (A, -p["half_width"], p["half_width"], bmev(p, 0), p["half_width"], p["sigma_eV"] * 1000, bmev(p, 1),
           s0d["real"], s0d["bloch"], s0d["eigen"], s0s["real"], s0s["bloch"], s0s["eigen"],
           s1d["real-bloch"], s1s["real-bloch"], tab))
    with open(os.path.join(a.out, "fig_stage_5.caption.txt"), "w") as f:
        f.write(caption + "\n")
    print(png)
    print(caption)


if __name__ == "__main__":
    main()
