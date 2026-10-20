# Handoff 0004 — tester report (stage 4 gate: G6 and G7 from Python; tutorials 02, 04, 05)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`), checkout
`~/lsquant-unified`. Branch `handoff/0004-codex` was cut from `cpp-unified` (`552de6e`), and
`origin/handoff/0004-agi` (`be54b79`) was merged into it as a fast-forward. Nothing under
`lsquant/`, `python/lsquant/`, `tests/dev/`, the legacy tree, `goldens/`, a seed or a tolerance
was edited. This handoff edited only `scripts/check_goldens.sh` (stage 4 added),
`figures/stage_4.py` (new) and `figures/out/stage_4/` (logs and figure). Every run used an
explicit `OMP_NUM_THREADS`. TCN was loaded throughout (load average 60 to 115 on 64 cores,
with other users' jobs). The developer's own check job in `~/lsquant-agi-build` also ran part
of the time. All wall times below are affected by this load.

**Result: GATE PASS (stage 4): 65 of 65 files IDENTICAL, all four greps empty.** That covers
stages 1 to 3 rerun, plus G6 (moments and 8 curves) and G7 (3 moments) at
OMP_NUM_THREADS 16 and 1.

**Tutorials from Python: 04 is 3 of 3 IDENTICAL. 02 is 0 of 193 and 05 is 0 of 16.** The
apps `tutorial02` and `tutorial05` exit 1 before computing anything, so no files are
produced (defect D1, section 5). With only the generator entry point aliased in a scratch
wrapper (not committed), the rest of both pipelines is byte-identical to the goldens:
02 is 193 of 193 and 05 is 16 of 16 (section 4.2).

## 1. Build (fresh)

```
$ rm -rf build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
cmake exit: 0
$ cmake --build build -j 32
[100%] Built target _core
build exit: 0
```

Full log: `figures/out/stage_4/build.log`. Warnings: **0 under `lsquant/`**. All 101
warnings are in the legacy tree, as in stages 1 to 3: 62 in `include/chebyshev_moments.hpp`,
32 in `include/units.hpp` and 7 in `src/chebyshev_momentsLocal.cpp`. `_core` links
`GSL::gsl` (`lsquant/CMakeLists.txt:51`). `lsquant_core` links only Eigen and OpenMP.

## 2. Developer self-check (`tests/dev/stage_4_check.py`, OMP_NUM_THREADS = 16)

```
$ PYTHONPATH=build/lsquant_core:python OMP_NUM_THREADS=16 python3 -u tests/dev/stage_4_check.py
(a) G7 chain1d_mag Gershgorin window = [-2.31, 2.31], BandWidth = 4.62
(a) G7 TimeProjSX-PX array_equal vs golden: True
(a) G7 TimeProjSY-PX array_equal vs golden: True
(a) G7 TimeProjSZ-PX array_equal vs golden: True
(b) G6 CorrelationVXchain1d array_equal vs golden: True
(b) G6 EF=+0.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+0.5 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.0 curve differing lines: 0 (total lines: 64)
(b) G6 EF=+1.5 curve differing lines: 0 (total lines: 64)
(c) compat GSL Bessel vs SciPy jv max_abs_diff = 5.551115e-16
(d) time_correlation_plan chain1d VX-VX (M=64, steps=16, tmax=5, exact) array_equal vs inline_compute-kpm-CorrOp: True
wall 655.12 s
exit: 0
```

Log: `figures/out/stage_4/stage_4_check.log`. The lines match the developer's report, which
ran at OMP_NUM_THREADS = 1. Item (d) ran the legacy binary from `build/` in this session.

## 3. Gate (`scripts/check_goldens.sh 4 build`)

Harness extension (this handoff): stage 4 runs stages 1, 2 and 3 in full, with their
comparisons and greps. It then runs `python3 -m lsquant.apps.stage4_time` at
`OMP_NUM_THREADS=16` and again at `1`, with `PYTHONPATH=build/lsquant_core:python`:

1. **G6 moments.** Command: `--operators goldens/stage0/chain1d_ballistic/operators --label
   chain1d --velocity VX --window -2 2 --trace exact --states 256 --moments 256 --steps 64
   --tmax 20 --word msd`. The produced `.chebmomTD` is compared with `cmp`.
2. **G6 curves.** Command: `--moments <file> --broadening-meV 10 --fermi-energy 0 0.5 1.0
   1.5`. The curves are reconstructed twice: from the committed golden `.chebmomTD`, and from
   the `.chebmomTD` produced in step 1. Both sets are compared with
   `goldens/stage0/chain1d_ballistic/mean...EF{0.000000,0.500000,1.000000,1.500000}JACKSON.dat`.
   This also checks that `--fermi-energy` formats `E_F` the way `std::to_string` does.
3. **G7 moments.** Command: `--operators goldens/stage0/chain1d_mag/operators --label
   chain1d_mag --operator SX|SY|SZ --projector PX --window gershgorin --trace exact --states
   128 --moments 32 --steps 96 --tmax 95 --word projected`. The three `.chebmomTD` files are
   compared with `cmp`.

   About the window: the brief's section 6 resolves it as the legacy `SpectralBounds` without
   a `BOUNDS` file. The developer identifies that as branch 3, the Gershgorin enclosure with
   pad 0.10, which gives `[-2.31, 2.31]`, BandWidth 4.62, matching the golden header
   `128 4.62 0 96 1`. The gate therefore passes `--window gershgorin`, not the literal
   numbers.
4. **Grep.** `grep -rnwE "time|velocity|spin|projector|bessel|gsl" lsquant/core lsquant/io`
   matches whole lowercase words. Comment lines (first non-blank characters `//`, `*` or `/*`)
   are excluded from the verdict and listed. The identifiers `TimeAxis`, `EvolveRegister`,
   `EvolveDerived` and `ApplyDerivedSeries` cannot match: the pattern is case-sensitive and
   whole-word.

Command: `OMP_NUM_THREADS=16 bash scripts/check_goldens.sh 4 build`. Wall time: 2458 s. Log:
`figures/out/stage_4/check_goldens_stage4.log`. The stage 1 to 3 part (41 files and three
greps) is identical, line for line in its verdict lines, to the stage-3 final run.

### 3.1 Stage-4 part, verbatim (the "printed:" lines are shortened here; the full text is in the log)

```
== stage 4: python3 -m lsquant.apps.stage4_time, PYTHONPATH=/home/ICN2/jgarcia/lsquant-unified/build/lsquant_core:python, vs goldens/stage0/{chain1d_ballistic,chain1d_mag}
          [T=16] G6 msd moments wall 140 s
IDENTICAL [T=16] chain1d_ballistic/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD (G6 moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=16] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat (G6 curve from new moments)
          [T=16] G7 SX-PX moments wall 45 s
          [T=16] G7 SY-PX moments wall 52 s
          [T=16] G7 SZ-PX moments wall 36 s
IDENTICAL [T=16] chain1d_mag/TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
IDENTICAL [T=16] chain1d_mag/TimeProjSY-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
IDENTICAL [T=16] chain1d_mag/TimeProjSZ-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
          [T=1] G6 msd moments wall 23 s
IDENTICAL [T=1] chain1d_ballistic/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD (G6 moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat (G6 curve from golden moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF0.500000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.000000JACKSON.dat (G6 curve from new moments)
IDENTICAL [T=1] chain1d_ballistic/meanCorrelationVXchain1dKPM_M256_state20exactEF1.500000JACKSON.dat (G6 curve from new moments)
          [T=1] G7 SX-PX moments wall 4 s
          [T=1] G7 SY-PX moments wall 4 s
          [T=1] G7 SZ-PX moments wall 4 s
IDENTICAL [T=1] chain1d_mag/TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
IDENTICAL [T=1] chain1d_mag/TimeProjSY-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
IDENTICAL [T=1] chain1d_mag/TimeProjSZ-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)
GREP      grep -rnwE "time|velocity|spin|projector|bessel|gsl" lsquant/core lsquant/io, code lines: empty
          excluded comment line: (27 comment lines, listed in the log; e.g. plan.hpp:7 "time-evolution step descriptors", plan.hpp:26 "... recursions such as mean-square", sparse_matrix.hpp:208 "purely imaginary velocity operators", and "O(N) time" complexity notes)
== 65 files compared; failures: 0
GATE PASS (stage 4)
gate exit: 0 (wall 2458 s)
```

Per-file count: stage 1 has 11, stage 2 has 6, stage 3 has 24, and stage 4 has 24
(12 at T=16 and 12 at T=1). Total 65, all IDENTICAL. Greps: `basis|std::cout|getenv` (core)
empty; `delta|DOS|Kubo|basis|std::cout|getenv` (core, io, code lines) empty; `compat/`
(core, io) empty; stage-4 physics words (core, io, code lines) empty. No `DIFFERS`, so the
localisation checklist of codex.md section 2 was not needed for the gate.

## 4. Figure (`figures/stage_4.py`)

Command: `MPLBACKEND=Agg python3 figures/stage_4.py --build build --threads 16`.

The script runs `python3 -m lsquant.apps.tutorial04|tutorial02|tutorial05 --out <tmp>`, each
in its own temporary directory, timed, with `KPM_SEED=12345` in the environment. Each app
sets its own trace seeds, as the legacy tutorial scripts do.

Matching rule for the byte comparison: every golden `captured/NNNN_<name>` is compared
(`filecmp`, `shallow=False`) with the new `captured/NNNN_<name>`, that is, the same index and
the same name. Tutorial 02 has 16 captures per name (the 16 disorder realisations), and each
is compared with its own counterpart. If a capture differs, it is also checked against every
new capture of the same name, to tell an ordering difference from a value difference. The
golden top-level `*.dat` files are compared by name.

A failing app does not stop the script: its files count as `DIFFERS (not produced)` and its
panel shows the golden only, with the error printed in the residual axis. The script exits
1 if any file differs.

Outputs:
- `figures/out/stage_4/fig_stage_4.png`
- `fig_stage_4.caption.txt`
- `fig_stage_4.cmp.txt` (the per-file cmp list, 212 entries)
- run log `figure_run.log` (exit 1)

### 4.1 Production run (the deliverable)

```
tutorial04: wall 217.9 s, 3 lines printed, 3 of 3 IDENTICAL (3 new captures)
tutorial02: wall 1.0 s, 0 lines printed, 0 of 193 IDENTICAL (0 new captures); tutorial02 exited 1: tutorial02 error: module 'lsquant_tutorial02_make_disordered_chain' has no attribute 'build'
tutorial05: wall 1.1 s, 0 lines printed, 0 of 16 IDENTICAL (0 new captures); tutorial05 exited 1: tutorial05 error: module 'lsquant_tutorial05_make_square' has no attribute 'build'
max |new - golden|: (a) 0.0e+00  (b) n/a (no new data)  (c) n/a (no new data); (a) max |golden - analytic| 3.997e-15
exit: 1
```

Tutorial 04 cmp list: `captured/0000_TimeProjSX-PX...chebmomTD`, `0001_...SY...` and
`0002_...SZ...` are all IDENTICAL. Tutorial 02: all 192 captures and the top-level
`meanCorrelationVXchain1d_dis_N512_W1KPM_M256_state600defaultEF0.000000JACKSON.dat` are
DIFFERS (not produced). Tutorial 05: the 8 captures and the 8 top-level `.dat` files are
DIFFERS (not produced).

Caption (verbatim, `fig_stage_4.caption.txt`):

> FIG. 4. Stage 4 (handoff 0004): tutorials 04, 02 and 05 reproduced from Python by python3 -m lsquant.apps.tutorial04|tutorial02|tutorial05 (lsquant package over the pybind11 module lsquant._core, OMP_NUM_THREADS = 16), each run in a temporary directory; time evolution U(dt) = exp(-i H dt/hbar) by the Bessel-Chebyshev series (GSL J_n via lsquant._core.compat, |J_n| <= 1e-15 cut, hbar = 0.6582119624 eV fs), dt = tmax/(steps - 1). (a) Tutorial 04, Larmor precession on chain1d_mag (N = 128, t = -1, J_ex = 0.1 eV): m = 0 row of TimeProj{SX,SY,SZ}-PX (projector PX = (1 + S_x)/2, exact trace over the 128 site states, M = 32, 96 time steps, tmax = 95 fs, Gershgorin window [-2.31, 2.31] eV), normalised by <S_x(0)> as lsqspin.py does; new as filled markers, golden as open grey markers (every 4th step), cos(wt), -sin(wt), 0 with w = 2 J_ex/hbar as thin lines; residual new - analytic below (golden - analytic max 4.00e-15, round-off level); max |new - golden| = 0.0e+00. (b) Tutorial 02, localisation in the 1D Anderson chain (box disorder W = 1 eV, disorder seed and KPM_SEED = r for r = 1..16, one random-phase vector, window from BOUNDS, M = 256, 128 steps, tmax = 600 fs, MSD(E,t) at 10 meV): xi(E) = cal sqrt(<MSD>) over the last 30% of the time window, averaged over the 16 realisations, N = 256, 1024, 4096, E = 0, +-1, 1.9 eV, cal fixed once against the Born length 24(4 - E^2)/W^2 at N = 4096, E = 1 eV (lsqloc.py); new as lines, golden as open markers, Born dashed; residual new - golden below (max n/a (no new data) sites). (c) Tutorial 05, 2D square lattice L = 200 (N = 40000), disorder W = 0.5, 1, 2, 4 eV (seed 1, KPM_SEED = 1, one random-phase vector): D(E,t) = DeltaX^2/(2t) at E = 2 eV from the MSD (M = 256, 160 steps, tmax = 150 fs, 80 meV; the app also writes the DOS, M = 1024, 80 meV); new as lines, golden as open markers (every 8th step); residual new - golden below (max n/a (no new data)). Byte comparison with goldens/stage0/tutorials/0N (captured/NNNN_<name> against the capture with the same index and name, top-level .dat by name; cmp criterion): tutorial 04 3 of 3 IDENTICAL; tutorial 02 0 of 193 IDENTICAL (app failed: tutorial02 exited 1: tutorial02 error: module 'lsquant_tutorial02_make_disordered_chain' has no attribute 'build'); tutorial 05 0 of 16 IDENTICAL (app failed: tutorial05 exited 1: tutorial05 error: module 'lsquant_tutorial05_make_square' has no attribute 'build'). Wall times (TCN, load 60 to 100 on 64 cores): tutorial 04 217.9 s, tutorial 02 1.0 s (failed), tutorial 05 1.1 s (failed). AT LEAST ONE FILE DIFFERS OR IS NOT PRODUCED.

### 4.2 Diagnostic run with the generator entry point aliased (localisation of D1; not a gate result)

Purpose: find out whether D1 is the only thing wrong with `tutorial02` and `tutorial05`.

Method: a scratch wrapper, `/tmp/s4_shim/shim_run.py` on TCN, not committed and not part of
the harness. It imports the app module unchanged and wraps its `_load_make_*()` so that the
loaded generator module gains one alias:

- `build(dir, n, W, seed) = build_chain(n, W, seed, outdir=dir)`
  (`make_disordered_chain.py`)
- `build(dir, L, W, seed) = build_square(L, W, seed, outdir=dir)` (`make_square.py`)

It then calls `main(["--out", <dir>])`. Nothing else is changed. The outputs were fed back
through the same figure script with `--reuse 04=... --reuse 02=... --reuse 05=... --name
fig_stage_4_diag_alias`.

Runs:
- tutorial05 at OMP_NUM_THREADS = 16: wall 624.4 s, concurrent with the gate.
- tutorial02 at OMP_NUM_THREADS = 1: wall 94.9 s. A first attempt at 16 threads had done
  only 28 of 192 realisations after about 10 minutes and was stopped (see S1).
- tutorial04 at OMP_NUM_THREADS = 1, unaliased, for comparison: wall 9 s.

```
tutorial04: reused /tmp/s4_shim/out04_T1, 3 of 3 IDENTICAL
tutorial02: reused /tmp/s4_shim/out02, 193 of 193 IDENTICAL (192 new captures)
tutorial05: reused /tmp/s4_shim/out05, 16 of 16 IDENTICAL (8 new captures)
max |new - golden|: (a) 0.0e+00  (b) 0.0e+00  (c) 0.0e+00
exit: 0
```

Files: `figures/out/stage_4/fig_stage_4_diag_alias.png`, `.caption.txt`, `.cmp.txt` and
`fig_stage_4_diag_alias_run.log`. That caption says OMP_NUM_THREADS = 16 because the figure
script was run with `--threads 16`, but the reused tutorial 02 and 04 outputs came from
1-thread runs, as listed above. The physics is visible in this figure:
- (a) `<S_x>`, `<S_y>` and `<S_z>` follow cos, -sin and 0, with golden - analytic at about
  4e-15.
- (b) xi(E) follows the Born parabola for N = 1024 and 4096, while N = 256 is finite-size
  limited at E = 0.
- (c) D(E,t) goes ballistic, then diffusive, then localised as W grows.

All of it is byte-identical to the goldens.

## 5. Defects found (localised; not fixed)

- **D1 (blocking for the tutorial deliverable): `tutorial02.py` and `tutorial05.py` call a
  generator function that does not exist.**
  - `python/lsquant/apps/tutorial02.py`, `_run_one_chain`, calls
    `make_mod.build(stage_dir, n_sites, w_val, disorder_seed)`.
  - `python/lsquant/apps/tutorial05.py`, `run_tutorial05`, calls
    `make_square.build(stage_dir, L_SIDE, w_val, DISORDER_SEED)`.
  - The example generators on `cpp-unified` define only
    `build_chain(N, W, seed, outdir=".")` (`examples/02_localization_1d_chain/make_disordered_chain.py:108`)
    and `build_square(L, W, seed, outdir=".")` (`examples/05_charge_transport_square/make_square.py:100`).
    Both take the directory as the keyword argument `outdir`, last, and return the label.
  - Both apps exit 1 with `AttributeError ... has no attribute 'build'` before any
    computation.
  - Section 5.2 shows that with the call alone corrected, both apps are byte-identical to
    the goldens: 193 of 193 and 16 of 16.
  - The developer's report shows only a `tutorial04` smoke run, so tutorials 02 and 05 were
    apparently never executed.

## 6. Smells

- **S1. Thread count and speed.** At OMP_NUM_THREADS = 16 on the loaded TCN, the
  time-evolution path is 6 to 25 times slower than at 1 thread:

  | Run | 16 threads | 1 thread |
  |---|---|---|
  | G7, per tensor | 36 to 52 s | 4 s |
  | G6 | 140 s | 23 s |
  | tutorial04 | 108 to 324 s (four runs) | 9 s |
  | tutorial02 (aliased) | under 3/192 per minute | 94.9 s total |

  - Cause: `SparseMatrix` SpMV (`lsquant/core/algebra/sparse_matrix.cpp:126`, stage 1) opens
    a `#pragma omp parallel for` on every product, with no minimum size. The time path
    performs many SpMVs on vectors of 128 to 4096 entries (the series terms plus `M` moments
    per step, per state).
  - Results are unaffected (byte-identical at 1 and 16).
  - The self-check took 655 s at 16 threads.
  - This is not a correctness issue, but it will matter for stage 8, and it makes
    "OMP_NUM_THREADS = 16" the slow setting for every small gate case.
- **S2. Physics word in the engine.** The identifier `TimeEvolutionStep`
  (`using TimeEvolutionStep = std::variant<EvolveRegister, EvolveDerived, CopyRegister>`,
  `lsquant/core/chebyshev/plan.hpp:150`, used in `moment_engine.cpp:190` and `:388`) is not
  in the allowed list (`TimeAxis`, `EvolveRegister`, `EvolveDerived`, `ApplyDerivedSeries`).
  It passes the codex.md grep only because that grep matches lowercase whole words.
  - The metadata fields `ChebyshevMoments::timeStep` and `numberOfTimeSteps` are named by
    the brief (section 1) and are not counted here.
  - Comment lines in `plan.hpp` speak of "time-evolution step descriptors" (line 7),
    "recursions such as mean-square [displacement]" (line 26) and "Time-evolution step ..."
    (lines 114, 125, 140).
- **S3. GSL included from `core/`: no.** `grep -rni gsl lsquant/core lsquant/io` is empty.
  `#include <gsl/gsl_sf_bessel.h>` appears only in `lsquant/compat/legacy_kernels.cpp`, and
  only `_core` links `GSL::gsl`.
- **S4. Clenshaw where forward accumulation was required: no.** `RunApplySeries` and
  `RunApplyDerivedSeries` (`moment_engine.cpp:43` to `83`) accumulate forward,
  `phi += c_n T_n` (`D_n`), with the recursion advancing, and apply the phase as a final
  `scal`. The derived series zeroes the `n = 0` term and follows the `Conv_Evolve` order.
  - Minor: both loops compute one recursion step after the last coefficient that is never
    used. This is one wasted SpMV, or two for the derived series, per series application.
- **S5. scipy Bessel on the gate path: no.** `time_evolution_series(..., use_legacy_bessel=True)`
  is the default, and `stage4_time`, `tutorial02/04/05` and the self-check (a), (b) and (d)
  all call it with the default. That routes through `legacy_bessel_coefficients`, which is
  `_core.compat.bessel_series_coefficients` (GSL `gsl_sf_bessel_Jn`). scipy `jv` is used only
  by the sympy `TimeEvolution` path and in self-check (c), where the max difference is
  5.6e-16.
- **S6. `tutorial04` depends on the legacy tree.** By default it reads
  `test/models/chain1d_mag/operators` (the legacy `test/` tree, to be deleted in stage 7),
  not `goldens/stage0/chain1d_mag/operators`. The five CSR files are byte-identical in both
  places (md5 checked), so the result is the same. This is the same pattern that stage 3
  review 2 removed from `tutorial03`.
- **S7. The `MSD` commutator operator is rescaled like `H`** (`RescaledOperator(VX, window)`,
  adding `ShiftFactor * I`). The developer records this as a reproduced legacy defect, and it
  is harmless here because all stage-4 windows are symmetric (BandCenter 0). This note keeps
  it visible for AGENTS.md section 10.
- **S8. `stage4_time` writes a second file.** Next to every legacy `.chebmomTD` it also writes
  the new-format `.chebmom`. The gate ignores it. `TimeAxis` ignores `--batch`: the time path
  has no batching.
- **No golden, seed or tolerance was touched.** `goldens/` is unchanged. The harness compares
  only with `cmp`.

## 7. Files committed on `handoff/0004-codex`

- `scripts/check_goldens.sh`: stage 4, which is `run_s4`, `stage4` and `stage4_grep`, plus
  the usage text.
- `figures/stage_4.py`
- `figures/out/stage_4/`:
  - `build.log`
  - `stage_4_check.log`
  - `check_goldens_stage4.log`
  - `figure_run.log`
  - `fig_stage_4.{png,caption.txt,cmp.txt}`
  - `fig_stage_4_diag_alias.{png,caption.txt,cmp.txt}`
  - `fig_stage_4_diag_alias_run.log`
- `handoff/0004-time-evolution/tests.md` (this file)
