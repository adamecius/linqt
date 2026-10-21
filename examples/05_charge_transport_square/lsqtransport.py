"""
lsqtransport.py -- charge transport on the 2D square Anderson lattice.

For each disorder strength W it builds the model, runs the KPM mean-square
displacement word `Tr[[X, U(t)]^+ delta(E-H) [X, U(t)]]` and DOS word
`Tr[delta(E-H)]` via `lsquant.runfile.execute_runfile`, and reconstructs
DeltaX^2(E,t) at the working energy E. From those it forms the running diffusion
coefficient D(E,t) = DeltaX^2/(2t), the semiclassical maximum D_max = max_t D,
the KPM density of states rho(E), the Einstein conductivity sigma = e^2 rho D_max,
and the mean free path ell = 2 D_max / v_F. It writes two figures:

  fig_diffusion.png  -- D(E,t) for every W (ballistic -> diffusive -> localized)
  fig_transport.png  -- sigma(W) and ell(W)

Everything is toy-scale and illustrative, not converged production numbers.
"""
import os
import sys
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import make_square

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

L = 200                       # linear lattice size (L*L sites)
E = 2.0                       # working energy (away from E=0 van Hove and band edges)
W_LIST = [0.5, 1.0, 2.0, 4.0]
M, NTIME, TMAX = 256, 160, 150.0
ETA_MEV = 80                  # reconstruction broadening
SEED = "1"


def _read_bounds():
    with open("BOUNDS", encoding="utf-8") as fh:
        b_lo, b_hi = (float(x) for x in fh.readline().split()[:2])
    return [b_lo, b_hi]


def diffusion(label):
    """Run the MSD word via `execute_runfile`, reconstruct DeltaX^2(E,t)."""
    execute_runfile(
        {
            "model": {
                "kind": "csr",
                "dir": os.path.abspath("operators"),
                "label": label,
                "operators": ["VX"],
            },
            "word": "Tr[[X, U(t)]^+ delta(E-H) [X, U(t)]]",
            "expansion": {
                "moments": M,
                "window": _read_bounds(),
                "trace": "random_phase",
                "states": 1,
                "seed": int(SEED),
            },
            "time": {
                "steps": NTIME,
                "tmax": TMAX,
                "bessel": "legacy",
            },
            "reconstruction": {
                "kernel": "jackson",
                "broadening_meV": ETA_MEV,
                "fermi_energies": [E],
                "formulas": ["time_correlation"],
                "legacy_text_round_trip": True,
            },
            "output": {
                "dir": os.path.abspath("."),
                "legacy_names": True,
                "legacy_layout": True,
            },
        }
    )
    dat = sorted(glob.glob("mean*%s*EF%f*dat" % (label, E)))[-1]
    t, dx2 = np.loadtxt(dat, unpack=True)
    good = t > 0
    return t[good], dx2[good] / (2.0 * t[good])


def dos_at(label, norb):
    """KPM DOS rho(E) per site for this disordered sample."""
    execute_runfile(
        {
            "model": {
                "kind": "csr",
                "dir": os.path.abspath("operators"),
                "label": label,
            },
            "word": "Tr[delta(E-H)]",
            "expansion": {
                "moments": 1024,
                "window": _read_bounds(),
                "trace": "random_phase",
                "states": 1,
                "seed": int(SEED),
            },
            "reconstruction": {
                "kernel": "jackson",
                "broadening_meV": ETA_MEV,
                "formulas": ["dos"],
                "legacy_text_round_trip": True,
            },
            "output": {
                "dir": os.path.abspath("."),
                "legacy_names": True,
                "legacy_layout": True,
            },
        }
    )
    dat = sorted(glob.glob("meanSpectralOp1%s*JACKSON.dat" % label))[-1]
    en, dos = np.loadtxt(dat, unpack=True)
    return float(np.interp(E, en, dos)) / norb


def fermi_velocity():
    """v_F = <|grad_k E|> over the constant-energy contour of the clean band."""
    k = np.linspace(-np.pi, np.pi, 600, endpoint=False)
    kx, ky = np.meshgrid(k, k)
    Ek = -2.0 * (np.cos(kx) + np.cos(ky))
    vmag = 2.0 * np.sqrt(np.sin(kx) ** 2 + np.sin(ky) ** 2)
    shell = np.abs(Ek - E) < 0.05
    return float(vmag[shell].mean())


def main():
    vF = fermi_velocity()
    curves, sigma, ell = {}, [], []
    for W in W_LIST:
        label = make_square.build_square(L, W, int(SEED))
        norb = L * L
        t, D = diffusion(label)
        Dmax = float(D.max())               # semiclassical maximum
        rho = dos_at(label, norb)
        sigma.append(rho * Dmax)            # sigma = e^2 rho D_max (e^2 = 1)
        ell.append(2.0 * Dmax / vF)
        curves[W] = (t, D)
        print("W=%.1f  D_max=%6.3f  rho(E)=%6.4f  sigma=%6.3f  ell=%6.2f"
              % (W, Dmax, rho, sigma[-1], ell[-1]))

    # figure 1: the ballistic -> diffusive -> localized crossover, in magnitude and shape.
    # An ORDERED disorder family -> a perceptually-uniform sequential colormap (viridis,
    # monotonic luminance so it also reads in greyscale) plus a colorbar.
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(W_LIST)))
    prl.use("web", height=5.8)
    fig, (axA, axB) = plt.subplots(2, 1, sharex=True)
    for (W, (t, D)), c in zip(curves.items(), colors):
        axA.plot(t, D, color=c, lw=1.6, label=r"$W=%g$" % W)
        axB.plot(t, D / D.max(), color=c, lw=1.6)
    axA.set_ylabel(r"$D(E,t)=\Delta X^2/2t$")
    prl.panel(axA, "(a)")
    axA.legend(title=r"disorder $W/t$", ncol=2, loc="upper right")
    axB.set_ylabel(r"$D(E,t)/D_{\max}$")
    axB.set_xlabel(r"time $t$ (fs)")
    axB.axhline(1.0, color="0.7", lw=0.6, zorder=0)
    prl.panel(axB, "(b)")
    prl.save(fig, "fig_diffusion")

    # figure 2: conductivity and mean free path vs disorder (log-y to span the decades).
    cS, cL = prl.COLORS[1], prl.COLORS[2]
    prl.use("web", aspect=0.52)
    fig, ax = plt.subplots()
    ax.semilogy(W_LIST, sigma, "o-", color=cS, label=r"$\sigma=e^2\rho D_{\max}$")
    ax.set_xlabel(r"disorder $W/t$")
    ax.set_ylabel(r"conductivity $\sigma$ (arb.)", color=cS)
    ax.tick_params(axis="y", labelcolor=cS)
    ax2 = ax.twinx()
    ax2.semilogy(W_LIST, ell, "s--", color=cL, label=r"$\ell=2D_{\max}/v_F$")
    ax2.axhline(1.0, color=cL, lw=0.7, ls=":")
    ax2.set_ylabel(r"mean free path $\ell$ ($a$)", color=cL)
    ax2.tick_params(axis="y", labelcolor=cL)
    ax2.spines["top"].set_visible(False)
    prl.save(fig, "fig_transport")


def figure_fermi():
    """Opening figure: the square-lattice band and the E=2t Fermi contour."""
    k = np.linspace(-np.pi, np.pi, 400); KX, KY = np.meshgrid(k, k)
    Ek = -2.0 * (np.cos(KX) + np.cos(KY))
    prl.use("web", aspect=0.85)
    fig, ax = plt.subplots(constrained_layout=True)
    cf = ax.contourf(KX, KY, Ek, levels=20, cmap="coolwarm", alpha=0.9)
    ax.contour(KX, KY, Ek, levels=[E], colors="k", linewidths=2.0)
    ax.set_xlabel(r"$k_x$"); ax.set_ylabel(r"$k_y$"); ax.set_aspect("equal")
    fig.colorbar(cf, ax=ax, label=r"$E_{\mathbf{k}}$ (eV)")
    fig.savefig("fig_fermi.png"); print("wrote fig_fermi.png")


if __name__ == "__main__":
    figure_fermi()
    main()
