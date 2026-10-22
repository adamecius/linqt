#!/usr/bin/env python3
"""Item 6: measure tables/GEMM versus FFT at identical nodes, reporting residuals.

The tester records timings, hashes and all residuals; this driver sets no tolerance,
does not certify a golden, and writes only to its explicitly supplied output path.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("moments")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--broadening-meV", type=float, default=10.0)
    parser.add_argument("--formula", choices=("bastin", "greenwood"), default="bastin")
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("positive --repeats required")
    import numpy as np
    from lsquant import _core
    from lsquant.reconstruct import kubo_fft_grid

    args.out.mkdir(parents=True, exist_ok=True)
    mu = _core.ReadChebmom(args.moments)
    report = {"engine_profile": _core.build_profile, "source_profile": mu.build_profile,
              "formula": args.formula, "broadening_meV": args.broadening_meV, "backends": {}}
    results = {}
    for backend in ("tables", "fft"):
        kubo_fft_grid(mu, args.broadening_meV, args.formula, backend=backend)
        timings = []
        for _ in range(args.repeats):
            start = time.perf_counter()
            energy, curve = kubo_fft_grid(mu, args.broadening_meV, args.formula, backend=backend)
            timings.append(time.perf_counter() - start)
        results[backend] = energy, curve
        report["backends"][backend] = {"wall_seconds": timings,
            "curve_sha256": hashlib.sha256(memoryview(np.ascontiguousarray(curve)).cast("B")).hexdigest()}
        np.savez(args.out / f"{backend}.npz", energy=energy, values=curve)
    a, b = results["tables"][1], results["fft"][1]
    difference = np.abs(a - b)
    report.update(grid_identical=bool(np.array_equal(results["tables"][0], results["fft"][0])),
                  max_absolute_difference=float(np.max(difference)),
                  relative_to_max_norm=float(np.max(difference) / max(float(np.max(np.abs(a))), 1e-300)),
                  max_pointwise_relative_difference=float(np.max(difference / np.maximum(np.abs(a), 1e-300))))
    (args.out / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
