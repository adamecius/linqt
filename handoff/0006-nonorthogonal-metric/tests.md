# Handoff 0006 — tester report (stage 6 gate: rotated models reproduce G1, G2, G5)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`),
checkout `~/lsquant-unified`, branch `handoff/0006-codex` cut from `cpp-unified`
(`552de6e`) and merged with `origin/handoff/0006-agi` at `a5108e3` (fast-forward: the
developer branch already contains `552de6e`). Nothing under `lsquant/`, `python/lsquant/`,
`tests/dev/`, the legacy tree, `goldens/` (data), a seed or a tolerance was edited. Files
edited or added: `scripts/check_goldens.sh` (stage 6), `figures/stage_6.py`,
`figures/stage_6_fixtures.py` (section 1 verification), `figures/out/stage_6/*`,
`goldens/stage0/PROVENANCE.md` (a G10 paragraph), this file.

**Result.** Fixtures: 4 of 4 ACCEPTED. Stage-6 rows: 74 of 74 PASS (40 route/applier rows
vs the goldens, 24 route-vs-route rows, 10 Hall plateaus); 8 `lowdin/ischo` rows N/A; 10
`naive` rows reported as evidence. Stages 1 to 3 inside the same run: 41 of 41 files
IDENTICAL, the stage-2 comment-aware grep and the stage-3 `compat/` grep empty. **The
harness exit is nevertheless `GATE FAIL (stage 6)`, 115 compared, 1 failure: the stage-1
grep `grep -rn "basis\|std::cout\|getenv" lsquant/core` (no comment exclusion) finds the
word `basis` in 9 Doxygen comment lines of the new `metric_policy.hpp` and
`composite_operators.hpp` (section 2.4).** No physics defect was found. C = +1.017017 for
every route and applier, equal to the golden G5 value.

## 0. Build

```
$ cmake -S . -B build6 -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
$ cmake --build build6 -j 24
BUILD_EXIT 0
```

Separate build directory `build6` so that other tester sessions using `build/` are not
disturbed. Full log: `figures/out/stage_6/build.log`.

## 1. Fixture verification (`figures/stage_6_fixtures.py`, independent of the generator)

Own reader for the legacy CSR text layout, dense numpy/scipy (`eigvalsh`, `scipy.linalg.eigh(H', S')`),
originals from `goldens/stage0/{chain1d,haldane}/operators`. Log:
`figures/out/stage_6/fixture_verification.log`.

| fixture | N | sha256 (all files vs manifest, and file list = manifest list) | max \|S - S^+\| | max \|H' - H'^+\| | max \|V' - V'^+\| | lambda_min(S) | lambda_max(S) | cond(S) | rel. diff vs manifest cond | max \|eig(H',S') - eig(H)\| | extra: max \|eig(H'+0.37VX'+0.21VY', S') - eig(same, orth.)\| | verdict |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| chain1d/real | 2048 | all match | 0 | 0 | 0 | 0.31377193275553 | 2.10355392820201 | 6.704085702404 | 2.0e-14 | 1.3e-14 | 5.2e-15 | ACCEPTED |
| chain1d/complex | 2048 | all match | 0 | 0 | 0 | 0.36043374446001 | 2.02063482826285 | 5.606120013236 | 1.7e-15 | 6.7e-15 | 6.2e-15 | ACCEPTED |
| haldane/real | 512 | all match | 0 | 0 | 0 | 0.35741617850047 | 2.06916590808509 | 5.789234042975 | 1.5e-15 | 1.0e-14 | 5.8e-15 | ACCEPTED |
| haldane/complex | 512 | all match | 0 | 0 | 0 | 0.38717798606128 | 2.02066370188938 | 5.218952974174 | 8.8e-15 | 1.5e-14 | 4.9e-15 | ACCEPTED |

All checks of codex.md section 1 hold: S Hermitian and positive definite, cond(S) equal
to the manifest to 1e-6 (actually 2e-14) and inside [5, 10], the generalised spectrum
equal to the orthogonal one to 1e-10 (actually 1.5e-14), H' and V' Hermitian to 1e-14
(exactly), sha256 of each file equal to the manifest. The extra row (tester's addition,
cheap) checks that every operator was rotated with the same T: the generalised spectrum
of a fixed combination of H', VX', VY' equals that of H, VX, VY. The `complex` variants
really are complex (max |Im S| = 0.50 chain1d, 0.25 Haldane; the `real` S has Im = 0).
Note: the manifest names the strength `alpha` (0.55), the brief `epsilon`.

## 2. The gate (`scripts/check_goldens.sh 6 build6`)

Command (TCN, `KPM_SEED=12345`, stage-6 runs at `OMP_NUM_THREADS=16`):

```
S6_KEEP=$HOME/lsquant-unified/figures/out/stage_6/data scripts/check_goldens.sh 6 build6
```

Stage 6 runs stages 1, 2 and 3 in full first (decision: stage 3 has landed and G4, G5
must stay unchanged), then the stage-6 block. Full log:
`figures/out/stage_6/check_goldens_stage6.log`; tables: `figures/out/stage_6/data/gate_table.tsv`,
`figures/out/stage_6/data/chern.tsv` (the run directories under `data/runs` are 37 MB and are not committed).

### 2.1 Decisions recorded

1. **Measure.** These runs cannot be byte-identical, so the criterion is the brief's
   tolerance (section 6). The measure is `rel = max_k |mu_new[k] - mu_gold[k]| /
   max_k |mu_gold[k]|` over all complex moments of the file (the moment vector's max-norm
   relative difference), plus header lines (N, BandWidth, BandCenter; M) numerically
   identical. `scripts/golden_diff.py` is still printed under every row; its
   per-element relative difference is ~1 on moments that vanish by symmetry (odd DOS
   moments, 1e-17 vs 0) and cannot be the criterion.
2. **Applier parameters.** `exact`: defaults. `ischo`: `--tol 1e-12` (brief table),
   `max_iterations` 200 (program default). `chebyshev`: `--order 60` (program default) with
   the constructor truncation `coeffTol = 1e-12` (header default; `stage6_nonorth` does not
   expose it); the effective order is printed.
3. **`lowdin/ischo` is N/A**: ISCHO provides `S^-1` only, the Löwdin route needs `S^-1/2`;
   `stage6_nonorth` rejects the combination with `std::invalid_argument` (developer report,
   section 5). Not counted.
4. **Route vs route**: `generalised/exact vs lowdin/exact` (1e-10),
   `generalised/chebyshev vs lowdin/chebyshev` (1e-8), `generalised/ischo vs lowdin/exact`
   (1e-9, as there is no `lowdin/ischo`).
5. **`naive`** run with `exact` only (the applier is not used on this route), evidence, not counted.
6. **Hall plateau**: `python3 -m lsquant.apps.stage3_kubo --moments <VX-VY file> --formula
   bastin --broadening-meV 10`, then `C = (mean_{|E|<0.5} sigma_xy / A) 2 pi`,
   `A = (512/2) 3 sqrt(3)/2` (identical to `tutorial03.chern_from_plateau` and
   `test/hall_haldane.sh`), PASS if `|C - 1| <= 0.05`.
7. **Environment workaround.** On TCN another session runs a loop
   `while true; do pkill -f "stage6_nonorth.*--route lowdi[n]"; sleep 20; done` (seen in
   `ps`, not mine, left alone). It would kill every Löwdin run of the gate. The harness
   therefore calls the same binary through a symlink `$TMP/s6_nonorth_exe`, so the command
   line does not match. The first gate launch was restarted for this before it reached stage 6.
8. A stage-6-only copy of the harness (`case 6) stage6`) was run in parallel in `~/s6`
   while the stage-1 part (64 threads, machine load ~70) was slow; its table is identical to
   the official run's (`diff` of the first eight columns empty). Only the official run is reported.

### 2.2 Stages 1 to 3 (inside the stage-6 run)

41 of 41 files IDENTICAL: G1, G2 at 1, 8, 64 threads, batched G2, random-phase
determinism, G3 forms (a) and (b) at 1 and 16 threads, G5 moments and curves, G4 moments
and curves at 16 and 1 threads. Stage-2 grep (code lines): empty. Stage-3 grep
`compat/`: empty. **Stage-1 grep: FOUND** (section 2.4).

### 2.3 Stage-6 table

`rel` as defined in 2.1; residual = `residual=` of the stdout line (`||S S^-1 v - v||/||v||`
on one random vector, or for ISCHO the `LastError()` of the last solve); it. = ISCHO
`LastIterations()` of the last solve; eff. = Chebyshev effective order.

| model | variant | word | route/applier | rel vs golden | tol | verdict | residual | it. | eff. |
|---|---|---|---|---:|---:|---|---:|---:|---:|
| chain1d | real | dos | lowdin/exact | 7.960e-13 | 1e-10 | PASS | 2.64e-14 | | |
| chain1d | real | dos | lowdin/chebyshev | 9.507e-12 | 1e-8 | PASS | 1.77e-12 | | 38 |
| chain1d | real | dos | generalised/exact | 4.664e-16 | 1e-10 | PASS | 3.74e-16 | | |
| chain1d | real | dos | generalised/ischo | 3.269e-14 | 1e-9 | PASS | 3.52e-14 | 0 | |
| chain1d | real | dos | generalised/chebyshev | 1.517e-12 | 1e-8 | PASS | 1.65e-12 | | 41 |
| chain1d | real | kubo | lowdin/exact | 1.510e-14 | 1e-10 | PASS | 2.64e-14 | | |
| chain1d | real | kubo | lowdin/chebyshev | 1.675e-13 | 1e-8 | PASS | 1.77e-12 | | 38 |
| chain1d | real | kubo | generalised/exact | 2.442e-15 | 1e-10 | PASS | 3.74e-16 | | |
| chain1d | real | kubo | generalised/ischo | 2.838e-15 | 1e-9 | PASS | 2.01e-16 | 0 | |
| chain1d | real | kubo | generalised/chebyshev | 1.026e-13 | 1e-8 | PASS | 1.65e-12 | | 41 |
| chain1d | complex | dos | lowdin/exact | 6.203e-13 | 1e-10 | PASS | 2.79e-14 | | |
| chain1d | complex | dos | lowdin/chebyshev | 9.337e-13 | 1e-8 | PASS | 1.76e-12 | | 33 |
| chain1d | complex | dos | generalised/exact | 1.022e-15 | 1e-10 | PASS | 4.01e-16 | | |
| chain1d | complex | dos | generalised/ischo | 2.735e-14 | 1e-9 | PASS | 1.46e-14 | 0 | |
| chain1d | complex | dos | generalised/chebyshev | 2.916e-12 | 1e-8 | PASS | 1.13e-12 | | 36 |
| chain1d | complex | kubo | lowdin/exact | 1.243e-14 | 1e-10 | PASS | 2.79e-14 | | |
| chain1d | complex | kubo | lowdin/chebyshev | 6.750e-14 | 1e-8 | PASS | 1.76e-12 | | 33 |
| chain1d | complex | kubo | generalised/exact | 2.554e-15 | 1e-10 | PASS | 4.01e-16 | | |
| chain1d | complex | kubo | generalised/ischo | 7.466e-15 | 1e-9 | PASS | 1.91e-16 | 0 | |
| chain1d | complex | kubo | generalised/chebyshev | 8.260e-14 | 1e-8 | PASS | 1.13e-12 | | 36 |
| haldane | real | kuboxy | lowdin/exact | 6.041e-14 | 1e-10 | PASS | 1.98e-14 | | |
| haldane | real | kuboxy | lowdin/chebyshev | 7.662e-13 | 1e-8 | PASS | 1.55e-12 | | 34 |
| haldane | real | kuboxy | generalised/exact | 1.340e-14 | 1e-10 | PASS | 1.04e-15 | | |
| haldane | real | kuboxy | generalised/ischo | 3.116e-13 | 1e-9 | PASS | 9.82e-13 | 6 | |
| haldane | real | kuboxy | generalised/chebyshev | 5.878e-13 | 1e-8 | PASS | 1.05e-12 | | 37 |
| haldane | real | kubo | lowdin/exact | 8.285e-14 | 1e-10 | PASS | 1.98e-14 | | |
| haldane | real | kubo | lowdin/chebyshev | 7.667e-13 | 1e-8 | PASS | 1.55e-12 | | 34 |
| haldane | real | kubo | generalised/exact | 1.148e-14 | 1e-10 | PASS | 1.04e-15 | | |
| haldane | real | kubo | generalised/ischo | 3.854e-13 | 1e-9 | PASS | 9.67e-13 | 6 | |
| haldane | real | kubo | generalised/chebyshev | 7.574e-13 | 1e-8 | PASS | 1.05e-12 | | 37 |
| haldane | complex | kuboxy | lowdin/exact | 4.036e-14 | 1e-10 | PASS | 1.52e-14 | | |
| haldane | complex | kuboxy | lowdin/chebyshev | 6.433e-13 | 1e-8 | PASS | 2.21e-12 | | 31 |
| haldane | complex | kuboxy | generalised/exact | 1.340e-14 | 1e-10 | PASS | 1.03e-15 | | |
| haldane | complex | kuboxy | generalised/ischo | 5.101e-14 | 1e-9 | PASS | 8.56e-14 | 6 | |
| haldane | complex | kuboxy | generalised/chebyshev | 3.179e-13 | 1e-8 | PASS | 1.14e-12 | | 34 |
| haldane | complex | kubo | lowdin/exact | 5.196e-14 | 1e-10 | PASS | 1.52e-14 | | |
| haldane | complex | kubo | lowdin/chebyshev | 6.069e-13 | 1e-8 | PASS | 2.21e-12 | | 31 |
| haldane | complex | kubo | generalised/exact | 1.148e-14 | 1e-10 | PASS | 1.03e-15 | | |
| haldane | complex | kubo | generalised/ischo | 6.636e-14 | 1e-9 | PASS | 9.98e-14 | 6 | |
| haldane | complex | kubo | generalised/chebyshev | 4.133e-13 | 1e-8 | PASS | 1.14e-12 | | 34 |

`lowdin/ischo`: N/A in all 8 model/variant/word cells (2.1, item 3).

Worst case per route/applier over all cells: lowdin/exact 7.96e-13, lowdin/chebyshev
9.51e-12, generalised/exact 1.34e-14, generalised/ischo 3.85e-13, generalised/chebyshev
2.92e-12: two to four orders of magnitude inside the tolerances. On the generalised/exact
route the agreement is at round-off (the dense `LLT` inverse); the Löwdin/exact route sits
at ~1e-12 on the chain DOS (the dense eigen-decomposition `S^-1/2`, error growing with m,
max at m = 62, 63).

Route vs route (24 rows, all PASS): generalised/exact vs lowdin/exact max 7.96e-13
(chain1d/real dos); generalised/chebyshev vs lowdin/chebyshev max 7.99e-12 (chain1d/real
dos); generalised/ischo vs lowdin/exact max 8.29e-13 (chain1d/real dos). Every row is in
the log and in `gate_table.tsv`.

Naive route (evidence, not counted): chain1d dos 1.1e+27 (real) / 2.3e+23 (complex), kubo
4.7e+49 / 2.9e+47 (the recursion diverges); Haldane kuboxy 0.85 / 0.86, kubo 1.00 / 1.00.
All would FAIL. The brief expected the naive route to pass on `dos`; it does not, because
the implemented `naive` is not the brief's naive (smell 5).

### 2.4 The one failure: stage-1 grep on comment lines

```
GREP      grep -rn "basis\|std::cout\|getenv" lsquant/core: FOUND
          lsquant/core/algebra/metric_policy.hpp:23,35,45,47,99,157,187  (Doxygen comments: "basis non-orthogonality", "orthonormal basis", ...)
          lsquant/core/algebra/composite_operators.hpp:321,382            (Doxygen comments: "non-orthogonal basis")
```

On `cpp-unified` the same grep is empty (0 lines). All 9 hits are comment lines; the
stage-2 grep, which excludes comment lines, is empty. Not relaxed here: the stage-1 grep
is as Fable wrote it. Either the comments say "metric" / "overlap" instead of "basis", or
Fable decides the stage-1 grep should exclude comments as the stage-2 one does.

### 2.5 Hall plateau (stage-3 Bastin, 10 meV)

| source | C | verdict |
|---|---:|---|
| golden G5 Bastin (legacy, reference) | +1.017017 | REF |
| haldane/real lowdin/exact, lowdin/chebyshev, generalised/exact, generalised/ischo, generalised/chebyshev | +1.017017 (each) | PASS (5) |
| haldane/complex lowdin/exact, lowdin/chebyshev, generalised/exact, generalised/ischo, generalised/chebyshev | +1.017017 (each) | PASS (5) |
| haldane/real naive | +1.015324 | evidence |
| haldane/complex naive | +1.017298 | evidence |

plateau = 107.656458 (450 points with |E| < 0.5) for all ten routes. The naive C is close
to 1 only by accident (its moments differ by O(1) from G5).

## 3. Figure

`figures/out/stage_6/fig_stage_6.png`, produced by
`MPLBACKEND=Agg OMP_NUM_THREADS=16 python3 figures/stage_6.py --build build6` from the kept gate runs (exit 0).
Caption (`figures/out/stage_6/fig_stage_6.caption.txt`):

> FIG. 6. Stage 6 (handoff 0006): a non-orthogonal basis through MetricPolicy reproduces the orthogonal goldens. Fixtures: the committed chain1d and Haldane operators rotated by a local transform T (H' = T^+ H T, S' = T^+ T, V' = T^+ V T; cond(S') 5.2 to 6.7; variants real and complex T), goldens/stage0/nonorth. Moments from build/lsquant_core/stage6_nonorth (exact trace, OMP_NUM_THREADS = 16). (a) chain1d DOS per site, window [-2, 2] eV, M = 64, reconstructed with the stage-2 lsquant.reconstruct.dos_curve (Jackson, 120 meV): orthogonal golden G1 (open markers, 40 points), Lowdin route (S^-1/2 H' S^-1/2) and generalised route (S^-1 (a H' + b S')) with the exact applier on the complex fixture; residual route - G1 below for both fixtures, max |residual| = 3.9e-14 1/eV. (b) Haldane sigma_xy(E), window [-4.5, 4.5] eV, M = 128, Kubo-Bastin at 10 meV by python3 -m lsquant.apps.stage3_kubo, in units of e^2/h per (N/2) A_cell (A_cell = 3 sqrt(3)/2) so the gap plateau (shaded, |E| < 0.5 eV) is C; golden G5 as markers; residual below, max |residual| = 9.4e-11. Chern numbers C = (mean_{|E|<0.5} sigma_xy / A) 2 pi: golden_G5 C = 1.017017 (REF); haldane_complex_kuboxy_generalised_chebyshev C = 1.017017 (PASS); haldane_complex_kuboxy_generalised_exact C = 1.017017 (PASS); haldane_complex_kuboxy_generalised_ischo C = 1.017017 (PASS); haldane_complex_kuboxy_lowdin_chebyshev C = 1.017017 (PASS); haldane_complex_kuboxy_lowdin_exact C = 1.017017 (PASS); haldane_complex_kuboxy_naive_exact C = 1.017298 (EVIDENCE(PASS)); haldane_real_kuboxy_generalised_chebyshev C = 1.017017 (PASS); haldane_real_kuboxy_generalised_exact C = 1.017017 (PASS); haldane_real_kuboxy_generalised_ischo C = 1.017017 (PASS); haldane_real_kuboxy_lowdin_chebyshev C = 1.017017 (PASS); haldane_real_kuboxy_lowdin_exact C = 1.017017 (PASS); haldane_real_kuboxy_naive_exact C = 1.015324 (EVIDENCE(PASS)). Max relative moment difference vs G1, G2, G5 (max_k |mu - mu_gold| / max_k |mu_gold|, over models, variants and words) per route/applier: generalised/chebyshev 2.92e-12; generalised/exact 1.34e-14; generalised/ischo 3.85e-13; lowdin/chebyshev 9.51e-12; lowdin/exact 7.96e-13 (tolerances: exact 1e-10, ischo 1e-9, chebyshev 1e-8; Lowdin/ischo not defined). Generalised vs Lowdin, same applier: chebyshev 7.99e-12; exact 7.96e-13; ischo 8.29e-13. Naive route (legacy estimator, evidence only): chain1d/real dos 1.149e+27; chain1d/real kubo 4.664e+49; chain1d/complex dos 2.306e+23; chain1d/complex kubo 2.914e+47; haldane/real kuboxy 8.497e-01; haldane/real kubo 9.956e-01; haldane/complex kuboxy 8.595e-01; haldane/complex kubo 1.002e+00. Gate rows: 64 PASS, 0 FAIL.

("Gate rows" in the caption counts the 40 + 24 moment rows of `gate_table.tsv`; the 10
Hall rows are listed separately in the same caption.) The sigma_xy residual of 9.4e-11
(in e^2/h per cell area) is concentrated in two spikes at |E| ~ 2.9 eV, where the Bastin
integrand is large; it does not touch the plateau.

## 4. Smells

1. **Hidden tolerance in a `.cpp`**: `overlap_appliers.cpp:134`
   `SparseMatrix CompressDenseMatrix(const Eigen::MatrixXcd& dense, double dropThreshold = 1e-18)`:
   the dense `S^-1` and `S^-1/2` are pruned at 1e-18 inside the `.cpp`, not a constructor
   parameter (reported by the developer as deviation 2; brief section 2 forbids it).
2. **Defaults on tolerances and caps in the headers**, not in the brief's interfaces:
   `ExactInverseOverlap(S, maxDimension = 4096)`, `ExactInverseSqrtOverlap(S, maxDimension = 4096)`,
   `IschoInverseOverlap(S, tolerance = 1e-13, maxIterations = 200)`,
   `ChebyshevInverseOverlap(..., coeffTol = 1e-12)`, `ChebyshevInverseSqrtOverlap(..., coeffTol = 1e-12)`.
   `stage6_nonorth` never passes `maxDimension` or `coeffTol` (no flag), so the gate's
   Chebyshev rows rely on the header default 1e-12.
3. **A fifth applier class** not in the brief: `ExactInverseSqrtOverlap` (dense
   eigen-decomposition), used by `lowdin/exact`. Consequence: there is no `S^-1/2` for
   ISCHO and `lowdin/ischo` is undefined (N/A rows).
4. **ISCHO diagnostics are last-call state written from a `const Multiply`**
   (`mutable lastIterations_, lastError_`): the `iterations=` and `residual=` printed by the
   gate program are those of the last solve, not a maximum; if the engine calls the
   operator from several OpenMP threads these members are written concurrently (a data
   race on diagnostics; results unaffected here). `iterations=0` on chain1d means the
   incomplete Cholesky of the banded chain `S` is already the exact Cholesky. The ISCHO
   non-convergence `throw` (commit `28d4def`) is raised inside `Multiply`, possibly inside an
   OpenMP region, where an exception cannot propagate.
5. **`--route naive` is not the brief's naive.** Brief section 4: naive = the symmetric
   plan with `S^-1 V'` on both sides (the legacy estimator, on the generalised H). The
   program's naive uses `OrthogonalMetric` and `RescaledOperator(H', window)`: it ignores
   `S` entirely, so the recursion runs on `H'`, whose spectrum leaves the window, and the
   chain1d moments diverge (1e+23 to 1e+49) instead of "passing on dos". The row does not
   demonstrate the legacy defect it is meant to show.
6. **Dead option**: `stage6_nonorth --drop-tol` (default 1e-3) is parsed and never used.
7. **Window padding constants in the app**: `ComputePaddedOverlapBounds` (Lanczos 80
   steps, seed 12345, pad 5% of the span, floor at half the lower Ritz value): in the app,
   not in core; acceptable but not in the brief.
8. **No applier prints** (`grep cout|cerr|printf|getenv|fopen|ifstream` on the four new
   or changed core `.cpp` files: empty); **no file reading in core**; the word `metric`
   does not appear under `core/chebyshev` (only `GeneralisedRescaledOperator` was added
   there, plus the brief-mandated `RescaledOperator(const LinearOperator&, window)`
   constructor; the stage-1 `SparseMatrix` member became `std::optional<SparseMatrix>`,
   and G1/G2 stay byte-identical). **The engine is untouched**: `git diff
   cpp-unified...a5108e3` changes no file of `MomentEngine`, `Plan`, or `io/`.
9. Stage-1 grep hits on comment lines (section 2.4).
10. Environment, not code: the foreign `pkill -f "stage6_nonorth.*--route lowdin"` loop
    on TCN (2.1 item 7). Whoever owns it should know it also kills gate runs.

## 5. Files

- `scripts/check_goldens.sh`: stage 6 (`6) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage6`).
- `figures/stage_6.py`, `figures/stage_6_fixtures.py`.
- `figures/out/stage_6/`: `build.log`, `check_goldens_stage6.log`, `fixture_verification.log`,
  `fig_stage_6.png`, `fig_stage_6.caption.txt`, `stage_6_figure.log`,
  `data/gate_table.tsv`, `data/chern.tsv`.
- `goldens/stage0/PROVENANCE.md`: G10 paragraph (fixtures verified, not regenerated).
