"""Random-phase (stochastic-trace) spin-relaxation evaluation for tutorial 06,
and its validation against the exact trace and against Elliott-Yafet theory.

Method. S_a(E,t) = Tr[P U(t)^+ S_a U(t) delta(E-H) P] / Tr[P U(t)^+ 1 U(t) delta(E-H) P].
Evaluated and reconstructed via `lsquant.runfile.execute_runfile` using run-file
configurations modelled on `runs/goldens/tutorial06.json`.
"""
import os
import sys
import glob
import numpy as np

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

from lsquant.runfile import execute_runfile

E = 2.0


def _read_bounds():
    if os.path.isfile("BOUNDS"):
        with open("BOUNDS", encoding="utf-8") as fh:
            b_lo, b_hi = (float(x) for x in fh.readline().split()[:2])
        return [b_lo, b_hi]
    return "gershgorin"


def _run_pair(label, N, op, proj, M, ntime, tmax, eta_mev, trace, states, seed):
    for pat in (
        f"TimeProj{op}-{proj}{label}*chebmomTD",
        f"TimeProj1-{proj}{label}*chebmomTD",
        f"mean*{op}-{proj}{label}*EF*dat",
        f"mean*1-{proj}{label}*EF*dat",
    ):
        for f in glob.glob(pat):
            os.remove(f)
    num_word = f"Tr[P:{proj} U(t)^+ {op} U(t) delta(E-H) P:{proj}]"
    den_word = f"Tr[P:{proj} U(t)^+ 1 U(t) delta(E-H) P:{proj}]"
    execute_runfile(
        {
            "model": {
                "kind": "csr",
                "dir": os.path.abspath("operators"),
                "label": label,
                "operators": [op, proj],
            },
            "basis": {"name": "real"},
            "metric": {"policy": "orthogonal"},
            "word": num_word,
            "words": [num_word, den_word],
            "expansion": {
                "moments": int(M),
                "window": _read_bounds(),
                "trace": str(trace),
                "states": int(states),
                "seed": int(seed),
            },
            "time": {
                "steps": int(ntime),
                "tmax": float(tmax),
                "bessel": "legacy",
            },
            "reconstruction": {
                "kernel": "jackson",
                "broadening_meV": float(eta_mev),
                "fermi_energies": [float(E)],
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
    num_dat = sorted(glob.glob("mean*%s-%s%s*EF%f*dat" % (op, proj, label, E)))[-1]
    den_dat = sorted(glob.glob("mean*1-%s%s*EF%f*dat" % (proj, label, E)))[-1]
    t, num = np.loadtxt(num_dat, unpack=True)
    _, den = np.loadtxt(den_dat, unpack=True)
    return t, num, den


def exact_spin_vs_time(label, N, op, proj, M, ntime, tmax, eta_mev=140):
    """Reference: exact trace over all N basis states via run-file execution."""
    t, num, den = _run_pair(
        label, N, op, proj, M, ntime, tmax, eta_mev,
        trace="exact", states=N, seed=12345,
    )
    good = t > 0
    s = num[good] / den[good].mean()
    return t[good], s / s[0]


def rp_spin_vs_time(label, N, op, proj, M, ntime, tmax, Nr=20, eta_mev=140, threads=16):
    """Stochastic trace over Nr random-phase vectors (seed = 1..Nr) via run-file execution."""
    NUM = DEN = t = None
    for r in range(1, Nr + 1):
        tt, num, den = _run_pair(
            label, N, op, proj, M, ntime, tmax, eta_mev,
            trace="random_phase", states=1, seed=r,
        )
        for f in glob.glob(f"TimeProj*-{proj}{label}*chebmomTD"):
            os.remove(f)
        NUM = num if NUM is None else NUM + num
        DEN = den if DEN is None else DEN + den
        t = tt
    good = t > 0
    s = NUM[good] / DEN[good].mean()
    return t[good], s / s[0]


def fit_rate(t, s, lo=0.2, hi=0.85):
    m = (s > lo) & (s < hi)
    if m.sum() < 3:
        return np.nan
    slope = np.polyfit(t[m], np.log(s[m]), 1)[0]
    return -1.0 / slope
