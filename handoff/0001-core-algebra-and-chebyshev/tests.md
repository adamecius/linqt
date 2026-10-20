# Handoff 0001 tests (Codex role, Opus subagent) — stage 1 gate on TCN

Machine: TCN (`fatnode-tacn`), checkout `~/lsquant-unified`, branch `handoff/0001-codex` =
`cpp-unified` (90b0bae, stage 0 closed) + `git merge --no-edit origin/handoff/0001-agi`
(merge commit a3a3d8e, no conflicts). g++ 11.4.0, OpenMP 4.5, Python 3.10, numpy 2.2.6,
matplotlib 3.10.9. Nothing under `lsquant/`, `tests/dev/`, the legacy tree, `goldens/`
or `scripts/` was modified.

## Summary

| Item | Result |
|---|---|
| Build (fresh `build/`) | exit 0; 0 warnings under `lsquant/` and `tests/dev/` |
| `stage_1_check` at 1 and 8 threads | dos 2.47e-15, kubo 1.14e-13, 1 vs 8 bit-identical YES |
| **Gate `scripts/check_goldens.sh 1 build`** | **did not run: `ERROR: stage1_moments not found`, exit 2** (section 3) |
| Same commands as the gate, run by hand (informational) | G1 and G2 IDENTICAL at T = 1, 8, 64; random-phase pair IDENTICAL; grep empty |
| Legacy `ctest` | 30/31, 30/31, 31/31 in three runs; the failure is a clock-dependent defect of the legacy test `graphene_reproduce.sh` (section 6), not the new core |
| Figure | `figures/out/stage_1/fig_stage_1.png`, script exit 0 |

The numerical result is clean. The gate itself does not pass as delivered, because a
fresh build does not put the executables where codex.md, the user instructions and the
harness expect them (finding F1). Fable decides whether that is fixed in the developer's
CMake or in the harness lookup; I did not change either.

## 1. Build

```
cd ~/lsquant-unified && rm -rf build
cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_INPUT_TOOLS=ON > cmake.log 2>&1
cmake --build build -j 32 > build.log 2>&1        # BUILD_EXIT=0
```

`build.log`: 101 lines containing `warning`, all in the legacy tree:
62 in `include/chebyshev_moments.hpp`, 32 in `include/units.hpp`, 7 in
`src/chebyshev_momentsLocal.cpp`. **Warnings under `lsquant/` and `tests/dev/`: none**
(`lsquant_core` and both executables are compiled with `-Wall -Wextra`).

Executable location after a fresh build: `build/stage1_moments` and `build/stage_1_check`.
`build/lsquant_core/` contains only CMake files (`CMakeFiles`, `Makefile`,
`cmake_install.cmake`). See F1.

## 2. `stage_1_check` (log: `figures/out/stage_1/stage_1_check.log`)

```
### OMP_NUM_THREADS=1 build/stage_1_check
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
exit=0
### OMP_NUM_THREADS=8 build/stage_1_check
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
exit=0
```

Identical to the developer's report. Both errors are below the 1e-11 of brief section 5.

## 3. The gate (log: `figures/out/stage_1/check_goldens_stage1.log`)

```
$ OMP_NUM_THREADS=16 bash scripts/check_goldens.sh 1 build 2>&1 | tee figures/out/stage_1/check_goldens_stage1.log
ERROR: stage1_moments not found under /home/ICN2/jgarcia/lsquant-unified/build/lsquant_core or /home/ICN2/jgarcia/lsquant-unified/build/lsquant
exit=2
```

**The gate did not pass: it did not run.** Localisation: the root `CMakeLists.txt` sets
`CMAKE_RUNTIME_OUTPUT_DIRECTORY` to `${CMAKE_BINARY_DIR}`, so `add_subdirectory(lsquant
lsquant_core)` links `stage1_moments` and `stage_1_check` into `build/`, not into
`build/lsquant/` (codex.md section 1) or `build/lsquant_core/` (the developer's report,
section 3, and the commands I was given). The developer's report says so and adds that
"symlinks in `build/lsquant_core/` are also present on TCN", i.e. made by hand in
`~/lsquant-agi-build`; a fresh build has none.

Two workarounds were tried and refused by the session's permission policy, so neither is
in this branch: (1) a fallback `EXE=$BUILD/stage1_moments` in the harness lookup (the
lookup only, not the criterion); (2) symlinks in the untracked `build/lsquant_core/`.
The harness and the build tree are as they were.

## 4. Same commands as the gate, by hand (informational, not a gate result)

To have the outputs the figure needs (and evidence for Fable), I ran exactly the
commands of codex.md section 2 with `build/stage1_moments`, `KPM_SEED=12345`, on
`goldens/stage0/chain1d/operators`, and compared with `cmp -s`. Script output verbatim
(`figures/out/stage_1/stage1_moments_manual_section2.log`):

```
[T=1] dos exit=0 : stage1_moments: word=dos dim=2048 states=2048 (exact) M=64 -> /home/ICN2/jgarcia/s1_manual/T1/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
[T=1] kubo exit=0 : stage1_moments: word=kubo dim=2048 states=2048 (exact) M=64x64 -> /home/ICN2/jgarcia/s1_manual/T1/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
IDENTICAL [T=1] SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL [T=1] NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
[T=8] dos exit=0 : ...
[T=8] kubo exit=0 : ...
IDENTICAL [T=8] SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL [T=8] NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
[T=64] dos exit=0 : ...
[T=64] kubo exit=0 : ...
IDENTICAL [T=64] SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
IDENTICAL [T=64] NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
[random 1] dos exit=0
[random 1] kubo exit=0
[random 2] dos exit=0
[random 2] kubo exit=0
IDENTICAL [random run2 vs run1] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom
IDENTICAL [random run2 vs run1] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D
IDENTICAL [random run2 vs run1] SpectralOp1chain1dKPM_M64_statedefault.chebmom
IDENTICAL [random run2 vs run1] SpectralOp1chain1dKPM_M64_statedefault.chebmom1D
IDENTICAL [random T=1 vs T=64] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom
IDENTICAL [random T=1 vs T=64] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D
IDENTICAL [random T=1 vs T=64] SpectralOp1chain1dKPM_M64_statedefault.chebmom
IDENTICAL [random T=1 vs T=64] SpectralOp1chain1dKPM_M64_statedefault.chebmom1D
GREP: 0 hits
```

(the `...` lines are the same `stage1_moments:` message with `T8`/`T64` in the path; the
log has them in full.)

| Golden | T = 1 | T = 8 | T = 64 |
|---|---|---|---|
| G1 `SpectralOp1chain1dKPM_M64_stateexact.chebmom1D` | IDENTICAL | IDENTICAL | IDENTICAL |
| G2 `NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D` | IDENTICAL | IDENTICAL | IDENTICAL |

No `cmp` differed, so no section-3 localisation was needed.

Random-phase determinism (`--trace random --seed 12345 --states 4`, T = 8, two runs):
all four outputs (legacy layout and `.chebmom`, both words) IDENTICAL. Additionally
T = 1 vs T = 64: IDENTICAL (see F3 for why this is weaker than it looks).

`grep -rn "basis\|std::cout\|getenv" lsquant/core`: empty.

The figure script (section 5) repeats the byte comparison independently
(`filecmp.cmp(shallow=False)`) at T = 1, 8, 64 and also finds every file IDENTICAL and
the `.chebmom` tensors bit-identical across thread counts.

## 5. Figure

`figures/stage_1.py` (run: `MPLBACKEND=Agg python3 figures/stage_1.py --threads 1 8 64`,
exit 0, log `figures/out/stage_1/stage_1_figure.log`). It runs `stage1_moments` itself
(lookup `build/lsquant_core`, `build/lsquant`, `build/`) into a temporary directory at
each thread count, byte-compares with `goldens/stage0/chain1d`, reads the new tensors
from the `.chebmom` files and the goldens from the legacy layouts, and computes the
residuals from the new legacy-layout files (same 15-digit text as the goldens).

Output: `figures/out/stage_1/fig_stage_1.png`, caption
`figures/out/stage_1/fig_stage_1.caption.txt`:

> FIG. 1. Stage 1 of the unified C++ core (handoff 0001): Chebyshev moments of the periodic tight-binding chain chain1d, N = 2048 sites, spectral window [-2, 2], exact trace over all 2048 site states, M = 64 moments, KPM_SEED = 12345, computed by build/stage1_moments on goldens/stage0/chain1d/operators at OMP_NUM_THREADS = 1, 8, 64. (a) Rank-1 moments mu_m (G1, word dos): legacy golden (open circles), new core (dashed line) and the analytic reference delta_m0 (triangles); right axis: |new - golden| (dash-dot), max 0.0e+00. (b) |new - golden| over the m, n <= 10 block of the rank-2 moments mu_mn (G2, word kubo, KPM weights and Hermitian symmetrisation applied), max 0.0e+00. (c) The closed-form block B.15, Mref(m, n) of test/chain1d_analytic_test.cpp times the diagonal plateau mean(mu_mm, m = 2..9) = 4 of the new tensor; max |Re mu_mn / plateau - Mref| over m, n < 10 is 0.0e+00, max |Im mu_mn| over the block is 0.0e+00. Byte comparison of the legacy-layout files with goldens/stage0/chain1d (cmp criterion): T=1 SpectralOp1chain1dKPM_M64_stateexact.chebmom1D: IDENTICAL; T=1 NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D: IDENTICAL; T=8 SpectralOp1chain1dKPM_M64_stateexact.chebmom1D: IDENTICAL; T=8 NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D: IDENTICAL; T=64 SpectralOp1chain1dKPM_M64_stateexact.chebmom1D: IDENTICAL; T=64 NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D: IDENTICAL. The new-format .chebmom tensors are bit-identical across the thread counts and agree with the 15-digit legacy-layout files to 0.0e+00. All files IDENTICAL.

Panel layout: codex.md asks for two panels; the second (the G2 block) is drawn as two
heat maps side by side, tagged (b) and (c). Style follows `examples/prl_style.py`
(double column, no title, tags inside the axes, Okabe-Ito colours with distinct
linestyles and markers).

## 6. Legacy `ctest` with the new core on

```
OMP_NUM_THREADS=16 KPM_SEED=12345 ctest --test-dir build --output-on-failure
```

| Run | Summary | Log |
|---|---|---|
| 1 | `97% tests passed, 1 tests failed out of 31` — `6 - graphene_regression (Failed)` | `figures/out/stage_1/ctest_run1.log` |
| 2 | `97% tests passed, 1 tests failed out of 31` — same test | `ctest_run2.log` |
| 3 | `100% tests passed, 0 tests failed out of 31` (83.98 s) | `ctest_run3.log` |

Failure output: `FAIL(bounds): estimate |E|max=58.682 out of sane range (1,50)` (run 1),
`|E|max=53.694` (run 2). In isolation the test passed every time (13 single-test runs
with and without `KPM_SEED` and `OMP_NUM_THREADS=16`, plus one `ctest -I 1,6`).

Localisation: `test/graphene_reproduce.sh` line 51 takes the maximum of every decimal
number in the solver log, `grep -oE '\-?[0-9]+\.[0-9]+'`. The log ends with a spdlog
timestamp `[hh:mm:SS.mmm] [INFO ] ... saved moments ...`, whose `SS.mmm` matches the
pattern. Whenever the run lands at second 50 or later of a minute, `|E|max` becomes
that second (58.682, 53.694) and the test fails; the real Gershgorin value is 7.51834.
This is a legacy test defect (about 1 run in 6 fails), independent of
`LSQUANT_BUILD_NEW_CORE`: the legacy drivers do not link `lsquant_core`. It did not show
in stage 0 or in the developer's run by chance. Not fixed here (legacy tree).

## 7. Smells and findings in the implementation

Read-only review of `origin/handoff/0001-agi`; nothing changed.

- **F1 (blocks the gate as delivered). Executable location.** A fresh build puts
  `stage1_moments` and `stage_1_check` in `build/`, while codex.md section 1
  (`build/lsquant/...`), the developer's own report (`build/lsquant_core/...`) and the
  harness expect a subdirectory. The developer's runs relied on hand-made symlinks in
  `~/lsquant-agi-build`. Fix options for Fable: set `RUNTIME_OUTPUT_DIRECTORY` on the two
  targets in `lsquant/CMakeLists.txt` (developer), or add `$BUILD/stage1_moments` to the
  harness lookup (tester, needs permission). The deviation is also wider than the
  brief's "exactly two lines": `add_subdirectory(lsquant lsquant_core)` instead of
  `add_subdirectory(lsquant)`, declared in the report.
- **F2. Environment variables read by the app.** `apps/stage1_moments.cpp` lines 49-54
  read `KPM_SEED` and `BATCH_SIZE` with `std::getenv`; `--seed` and `--batch` override
  them. Brief section 2.3 says these are read by the app "as command-line options".
  Not in `core/` (grep is clean), but a stray `KPM_SEED` in the environment silently
  changes the random-phase seed when `--seed` is absent, and `BATCH_SIZE` silently
  changes the batch when `--batch` is absent. Fable to rule whether the env fallback is
  wanted.
- **F3. The thread-count checks never reach the parallel kernel.**
  `SparseMatrix::Multiply` runs the OpenMP loop only `if (rp[n] > kMinParallelNnz)`,
  `kMinParallelNnz = 30000` (`core/algebra/sparse_matrix.cpp:14`, same threshold as
  legacy `LSQ_SPMV_MIN_PARALLEL_NNZ`). chain1d `HAM` and `VX` have 4096 non-zeros
  each (at most 6144 after rescaling adds the diagonal) and the self-check 64 x 64 at
  most 4096, so at T = 1, 8, 64 every
  matvec in the gate and in `stage_1_check` is sequential. "Bit-identical across thread
  counts" is therefore vacuous for both; the per-row parallel path is untested in this
  stage. Suggest one determinism run on a matrix above the threshold in a later gate
  (the committed chain1d and graphene golden operators do not qualify: 4096 and 3328
  non-zeros; a longer chain, N above about 15000, would); I did not add it (not asked
  for by name).
- **F4. Batching is not exercised.** The gate uses the default `batchSize = M`; the
  `PaddedMoments` round-up and the block-order logic in `moment_engine.cpp` (lines
  51-53, 193, 223, 252-253) with `batchSize < M` or a non-divisor of `M` are not
  compared with the legacy anywhere.
- **F5. `ReadCsrText` Hermitian-folds every operator.** `io/csr_text.cpp:56` applies
  `HermitianFold()` (`0.5 * (A^+ + A)`) unconditionally, mirroring legacy
  `ConvertFromCSR`; this is what makes G2 byte-identical (the `-0` imaginary parts).
  It is declared in the report. Consequence: any non-Hermitian operator file is
  silently replaced by its Hermitian part. Fine for stage 1, a trap for later words
  with non-Hermitian operators.
- **F6. Rank-specific branches in `MomentEngine::Run`.** The engine dispatches to four
  hand-written paths (rank 0, rank 1 ket, rank 1 bra, rank 2; `moment_engine.cpp`
  around lines 167, 177, 221, 249) rather than one rank-generic loop. Within stage-1
  scope (rank <= 2 is all the brief allows), but noted as the place stage 3 will
  have to touch.
- **F7. Trace-weight sign-of-zero special case.** `Trace` uses `weight.real() * d` when
  `weight.imag() == 0.0` to reproduce the legacy `double * complex` signed zeros
  (declared in the report). Correct and necessary for byte identity; noted because a
  weight with a genuinely zero imaginary part and a non-real intent cannot be told
  apart.
- **F8. `GershgorinBounds` is symmetric** (`[-rho(1+pad), +rho(1+pad)]`,
  `spectral_bounds.hpp:62`), as the legacy safeguard is; for a spectrum not centred at
  zero it wastes window. Not used by the gate (the window is given).
- Checked and clean: no `TODO`/`FIXME` anywhere under `lsquant/` or `tests/dev/`; no
  `getenv`, `std::cout`, `std::cerr`, `printf` under `core/` or `io/`;
  `std::ifstream`/`ofstream` only in `io/`; no physics names under `core/` other than
  `VX` inside one doc comment (`sparse_matrix.hpp:179`, explaining the signed zeros)
  and `delta_{m0}` in doc comments of `chebyshev_moments.hpp` (the KPM weight, not a
  word); `ApplySeries` throws the required message; out-of-scope classes
  (`BlochTransform`, `MetricPolicy`, `CoefficientTable`, ...) absent; Eigen and OpenMP
  are the only link dependencies of `lsquant_core`.

## 8. What I changed

Added on `handoff/0001-codex`: `figures/stage_1.py`; `figures/out/stage_1/`
(`fig_stage_1.png`, `fig_stage_1.caption.txt`, `check_goldens_stage1.log`,
`stage1_moments_manual_section2.log`, `stage_1_check.log`, `stage_1_figure.log`,
`ctest_run1.log`, `ctest_run2.log`, `ctest_run3.log`); this file. No golden, seed,
tolerance, harness criterion or source file was changed.
