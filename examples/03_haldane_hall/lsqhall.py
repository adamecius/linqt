#!/usr/bin/env python3
"""lsqhall - Hall and longitudinal conductivity of the Haldane model, over LSQUANT.

The heavy work runs through the unified run-file engine (`lsquant.runfile.execute_runfile`),
which evaluates the two-operator Kubo words

    Tr[VX delta(E-H) VY dG-(E)/dE]
    Tr[VX delta(E-H) VX delta(E-H)]

and reconstructs `bastin`, `sea`, `surf` (from VX-VY) and `greenwood` (from VX-VX).

The Haldane model (t1 = -1, t2 = 0.15, phi = pi/2) is a Chern insulator with
C = +1. Units follow the golden: e = hbar = 1, so e^2/h = 1/(2 pi).
"""

import os
import sys
import glob
import tempfile

import numpy as np

# Shared publication (APS/PRL) plot style: examples/prl_style.py.
_EX_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_REPO_ROOT = os.path.dirname(_EX_DIR)
for _p in (
    _EX_DIR,
    os.path.join(_REPO_ROOT, "python"),
    os.path.join(_REPO_ROOT, "build", "lsquant_core"),
    os.path.join(_REPO_ROOT, "build"),
):
    if os.path.isdir(_p) and _p not in sys.path:
        sys.path.insert(0, _p)
import prl_style as prl
from lsquant.runfile import execute_runfile

A_CELL = 3 * np.sqrt(3) / 2          # honeycomb unit-cell area
E2_OVER_H = 1.0 / (2 * np.pi)        # conductance quantum in e = hbar = 1 units


def _read_bounds(opdir, label):
    desc_path = os.path.join(opdir, f"{label}.HAM.desc")
    if os.path.isfile(desc_path):
        vals = {}
        with open(desc_path, encoding="utf-8") as fh:
            for line in fh:
                if ":" in line and not line.lstrip().startswith("#"):
                    k, v = line.split(":", 1)
                    vals[k.strip()] = v.strip()
        if "spectral_min" in vals and "spectral_max" in vals:
            return [float(vals["spectral_min"]), float(vals["spectral_max"])]
    return [-4.5, 4.5]


def haldane_conductivities(opdir, label="haldane", M=128):
    """Return dim and the conductivity curves as {name: (E, sigma)}.

    Each curve keeps its OWN energy grid: the Kubo-Bastin route (sigma_xy, sea, surf) uses a
    fine 30*M grid, while Kubo-Greenwood (sigma_xx) uses an M-point grid, so they are not the
    same length.
    """
    dim = int(open(os.path.join(opdir, f"{label}.HAM.CSR")).readline().split()[0])
    win = _read_bounds(opdir, label)
    with tempfile.TemporaryDirectory() as wd:
        execute_runfile(
            {
                "model": {
                    "kind": "csr",
                    "dir": os.path.abspath(opdir),
                    "label": label,
                    "operators": ["VX", "VY"],
                },
                "word": "Tr[VX delta(E-H) VY dG-(E)/dE]",
                "words": [
                    "Tr[VX delta(E-H) VY dG-(E)/dE]",
                    "Tr[VX delta(E-H) VX delta(E-H)]",
                ],
                "expansion": {
                    "moments": int(M),
                    "window": win,
                    "trace": "exact",
                },
                "reconstruction": {
                    "kernel": "jackson",
                    "broadening_meV": 10.0,
                    "formulas": ["bastin", "sea", "surf", "greenwood"],
                },
                "output": {
                    "dir": wd,
                    "legacy_names": True,
                    "legacy_layout": True,
                },
            }
        )
        d_xy = np.loadtxt(sorted(glob.glob(os.path.join(wd, f"KuboBastin_NonEqOpVX-VY{label}*.dat")))[-1])
        d_sea = np.loadtxt(sorted(glob.glob(os.path.join(wd, f"KuboBastinSea_NonEqOpVX-VY{label}*.dat")))[-1])
        d_surf = np.loadtxt(sorted(glob.glob(os.path.join(wd, f"KuboBastinSurf_NonEqOpVX-VY{label}*.dat")))[-1])
        d_xx = np.loadtxt(sorted(glob.glob(os.path.join(wd, f"KuboGreenWood_NonEqOpVX-VX{label}*.dat")))[-1])
    return dim, {
        "sigma_xy": (d_xy[:, 0], d_xy[:, 1]),
        "sea": (d_sea[:, 0], d_sea[:, 1]),
        "surf": (d_surf[:, 0], d_surf[:, 1]),
        "sigma_xx": (d_xx[:, 0], d_xx[:, 1]),
    }


def chern_from_plateau(dim, E, sigma_xy, gap_halfwidth=0.5):
    """Chern number from the gap plateau: C = (plateau / A) * 2 pi.

    Matches hall_haldane.sh: A = (dim/2) * A_cell, e = hbar = 1. The Kubo-Bastin sigma_xy is
    area-extensive, so the plateau divided by A and multiplied by 2 pi reads off the integer.
    """
    area = (dim / 2) * A_CELL
    plateau = np.mean(sigma_xy[np.abs(E) < gap_halfwidth])
    return (plateau / area) * 2 * np.pi


def figure(opdir, label="haldane", M=128):
    import matplotlib.pyplot as plt
    dim, c = haldane_conductivities(opdir, label, M)
    norm = (dim / 2) * A_CELL * E2_OVER_H          # divide to read off e^2/h units
    Exy, sxy = c["sigma_xy"]
    C = chern_from_plateau(dim, Exy, sxy)

    prl.use("web", height=5.8)
    fig, (top, bot) = plt.subplots(2, 1, sharex=True)

    for i, (key, lab) in enumerate([("sigma_xy", "Bastin (total)"), ("sea", "Fermi sea"),
                                    ("surf", "Fermi surface")]):
        E, s = c[key]
        top.plot(E, s / norm, label=lab, **prl.trace(i))
    top.axhline(1.0, color="0.6", ls=":", lw=0.8)
    top.set_ylabel(r"$\sigma_{xy}$ ($e^2/h$)")
    prl.panel(top, "(a)", loc="upper right")
    top.legend(loc="upper left", title=r"$C=%+.3f$ (gap plateau)" % C)

    Exx, sxx = c["sigma_xx"]
    bot.plot(Exx, sxx / norm, **prl.trace(3))
    bot.set_xlabel(r"energy $E$ (eV)")
    bot.set_ylabel(r"$\sigma_{xx}$ ($e^2/h$)")
    prl.panel(bot, "(b)")

    prl.save(fig, "fig_haldane_hall")
    print(f"Chern from gap plateau: C = {C:+.3f}")


def demo():
    opdir = os.environ.get("HALDANE_OPDIR", "operators")
    if not os.path.exists(os.path.join(opdir, "haldane.HAM.desc")):
        import make_haldane                                 # stage the golden operators + .desc
        make_haldane.build()
    figure(opdir, label="haldane", M=128)


if __name__ == "__main__":
    demo()
