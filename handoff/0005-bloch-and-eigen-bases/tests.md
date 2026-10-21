# Handoff 0005 — tester report (stage 5 gate: Bloch and eigen bases, G9; G1 to G8 unchanged)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`), worktree
`~/lsquant-unified-0005` of the clone `~/lsquant-unified`, branch `handoff/0005-codex` cut
from `origin/cpp-unified` (`8c3ae32`) and fast-forwarded to `origin/handoff/0005-agi`
(`7109c20`, no merge commit). Own build directory `~/lsquant-unified-build5`; the other
testers' branches and builds were not touched. Nothing under `lsquant/`, `python/lsquant/`,
`tests/dev/`, the legacy tree, a seed, a tolerance or an existing golden was edited. Edited or
added: `scripts/check_goldens.sh` (stage-4 harness taken verbatim from
`origin/handoff/0004-codex` `3dac161`, plus stage 5), `figures/stage_5.py`,
`figures/stage_5_table.py`, `figures/out/stage_5/`, `goldens/stage0/PROVENANCE.md` (G9
section), `goldens/stage0/bases/` (new), this file.

**Result: GATE FAIL (stage 5), 149 rows compared, 37 failures** (exit 1). Every failure is
in one of two groups (16 G9 rows at each of two thread counts = 32, plus 5 D3 rows):

1. **G9 rows (basis vs hr_exactdiag), 16 of 16 FAIL at every thread count**, by 2.9e-2 to
   8.1e-1 against the codex tolerances 2e-2 (DOS) and 5e-2 (sigma_xx). No lsquant defect is
   behind them: the eigen and exact-Bloch moments equal the exact moments of the same 48 x 48
   system computed from hr_exactdiag's own `H(k)` and the wannier2sparse velocity to
   1e-14 (section 4.3). The difference is the oracle's shifted mesh (i + 1/2)/48 against the
   supercell's Gamma mesh n/48, and the Jackson kernel against a Gaussian of width pi a / M;
   the kernel term alone (same mesh, same eigenvalues) is 3.3e-2 (graphene DOS) and 2.8e-1
   (graphene sigma_xx), above the tolerances by itself. The tolerance row cannot pass at
   N = 48, M = 256 as specified; Fable's decision.
2. **Defect D3 (5 rows): the real-space supercell and `B` use the cell index
   `n1 N2 + n2` (brief section 3); wannier2sparse writes `n1 + N1 n2`.** lsquant's
   supercell CSR differs from wannier2sparse's for any model or supercell that is not
   symmetric under n1 <-> n2: graphene 6 x 4 x 1, Haldane 6 x 4 x 1 and 48 x 48 x 1;
   `B VX_w2s B^+ = v_x` (brief section 3, 1e-10) fails for Haldane at 48 x 48 by 9.1e-1.
   Renumbering lsquant's cells n1-fastest makes all of them agree to 2e-16. It does not
   show up in the three-basis rows, because RealBasis builds its own CSR (it does not read
   wannier2sparse's), and DOS and sigma_xx are invariant under the mirror n1 <-> n2.

Everything else passes. Stages 1 to 4: **65 of 65 IDENTICAL** (the stage-4 harness of
`handoff/0004-codex` at `3dac161`; that branch has since added G4c, 71 files, not included
here), all greps empty. Real vs Bloch moments: 2.5e-14 (W = 0) and 2.2e-14 (W = 1, disorder
seed 7) against 1e-9; Bloch (exact) vs eigen curves: 5.7e-15 or better against 1e-8; the
FFTW `B` equals the dense `B` to 3.6e-15 and `B^+ B = 1` to 3.9e-16. Every stage-5 output is
byte-identical between OMP_NUM_THREADS = 1 and 16 in all three bases (16 real files counted,
32 Bloch and eigen files reported). Two reader defects (D1, D2) prevent stage5_bases from
reading the wannier2sparse example files at all; the gate feeds it a format-adapted copy of
the same model (section 2.3).

## 1. Build (fresh) and FFTW

```
$ rm -rf ~/lsquant-unified-build5
$ cmake -S . -B ~/lsquant-unified-build5 -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
-- Checking for module 'fftw3'
--   Found fftw3, version 3.3.8
cmake exit: 0
$ cmake --build ~/lsquant-unified-build5 -j 32
[100%] Built target _core
build exit: 0
```

FFTW reaches `lsquant_core` through `pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)`
in the top-level `CMakeLists.txt` and `target_link_libraries(lsquant_core PUBLIC ...
PkgConfig::FFTW3)` in `lsquant/CMakeLists.txt`; no extra flag. `fftw3.h` is included only in
`lsquant/core/algebra/bloch_transform.cpp` (the header hides the plan behind a workspace
pointer). Warnings: 101, all in the legacy tree (62 `include/chebyshev_moments.hpp`, 32
`include/units.hpp`, 7 `src/chebyshev_momentsLocal.cpp`), **0 under `lsquant/`**. One CMake
warning: `cmake/FetchWannier2Sparse.cmake:9 wannier2sparse unavailable (offline?)` (TCN has
no network for the fetch; harmless here). Log: `figures/out/stage_5/build.log`.

## 2. Developer self-check and the inputs

### 2.1 `tests/dev/stage_5_check.py`, OMP_NUM_THREADS = 16 (exit 0, 2 min 48 s)

```
[check a] BlochTransform unitarity (8,8,1) graphene: max|B B^+ x - x| = 6.206e-17, max|B^+ B x - x| = 5.004e-17
[check b] (graphene 8x8x1) max_x ||H_k x - B H_csr B^+ x||_inf = 7.109e-16
[check c] (graphene 8x8x1) max_x ||v_x x - B VX_csr B^+ x||_inf = 7.109e-16
[check b] (haldane 8x8x1) max_x ||H_k x - B H_csr B^+ x||_inf = 5.496e-16
[check c] (haldane 8x8x1) max_x ||v_x x - B VX_csr B^+ x||_inf = 6.474e-16
[check d] Anderson disorder W=1.0 (seed=12345): max|mu_real - mu_bloch| DOS(64) = 4.680e-16, VX-VX(64x64) = 5.761e-14
[check e] (16,16,1) real vs bloch: DOS max|d_mu|=8.105e-15, max_rel|d_curve|=6.915e-15; KuboXX max|d_mu|=5.898e-13, max|d_curve/N|=5.313e-14, max_rel|d_curve|=1.655e-14
[check e] (16,16,1) real vs eigen: DOS max|d_mu|=7.883e-15, max_rel|d_curve|=6.825e-15; KuboXX max|d_mu|=6.040e-13, max|d_curve/N|=4.588e-14, max_rel|d_curve|=1.429e-14
[check e] (16,16,1) bloch vs eigen: DOS max|d_mu|=1.222e-15, max_rel|d_curve|=1.203e-15; KuboXX max|d_mu|=1.279e-13, max|d_curve/N|=7.256e-15, max_rel|d_curve|=2.260e-15
[check e] Jackson KPM (M=128) vs Gaussian exact-diag: rel_L2(DOS)=2.130e-02, rel_L2(KuboXX)=4.250e-02
[check f] G1 chain1d DOS byte-identical to golden: PASS
[check f] G4 graphene Kubo byte-identical to golden: PASS
[check f] G7 chain1d_mag spin precession byte-identical to golden: PASS
stage_5_check: ALL CHECKS PASSED

real	2m48.152s
EXIT:0
```

Identical to the numbers in `report.md` section 3.1. What it does not check, against brief
section 5: (a) is unitarity only, not the FFTW `B` against a dense numpy `B` (the gate does
it, section 4.2); (b), (c) compare `B` with lsquant's own supercell CSR, not with the
wannier2sparse CSR, on a square 8 x 8 supercell, which is why D3 passes there; the models
are written by the check itself (dense `_hr.dat`, `.xyz` with a comment line), not taken from
`gen_models.py` with `N = 24` as the brief says, which is why D1 and D2 pass there; (d) draws
the disorder with the random-phase seed; (e) uses the window [-10, 10] (not the eigenvalue
window); (f) is replaced by G1/G4/G7 regressions instead of `sigma_xx` Bloch vs eigen at
E = 0.5.

### 2.2 G9 goldens (generated by the tester; `goldens/stage0/PROVENANCE.md`, section G9)

`goldens/stage0/bases/gen_g9.sh` (committed) does, for graphene and Haldane from
`wannier2sparse/examples/gen_models.py` (wannier2sparse `3f553db`; the pinned `26cab4a` has no
`tools/hr_exactdiag.py`):

- velocities: `wannier2sparse <m> 1 1 1 VX VY --mode bundle` (berry_connection, eV A),
  copied as `<m>_vx_hr.dat`, `<m>_vy_hr.dat` (the names `hr_exactdiag.py sigma` reads);
- window: eigenvalues of `H(k)` on the 48 x 48 Gamma mesh, both models [-3, 3] eV, padded
  10 %: **[-3.3, 3.3] eV, a = 3.3 eV**;
- broadening: **sigma = pi a / M = pi 3.3 / 256 = 0.0404970928 eV = 40.497 meV** (the
  Jackson resolution of M = 256 on that window); the lsquant runs use
  `--broadening-meV 40.4970928`, which keeps all 256 moments (`M_eff = min(ceil(pi/eta), M)`);
- `hr_exactdiag.py dos|sigma --comp xx <m> --nk 48 --eta 0.0404970928 --ngrid 2641
  --emin -3.3 --emax 3.3`.

Committed per model: `G9_dos.json`, `G9_sigma_xx.json`, `G9_params.txt`, `model/` (the five
input files). A second generation into a scratch tree reproduced all 16 files byte for byte.
Oracle caveats (recorded, not changed): hr_exactdiag's mesh is the shifted (i + 1/2)/48 mesh,
not the supercell's; `cmd_sigma` divides by `|a1_x a2_y|` = 1.299 A^2, half of the
honeycomb cell, so its `sigma_xx_e2h` is twice the e^2/h value for these models (its own
`cell_area_2d` documents the pitfall; `cmd_sigma` does not use it).

### 2.3 Format adapter (decision; the defects it works around are D1 and D2)

`stage5_bases` cannot read the gen_models / wannier2sparse files:

```
stage5_bases error: read_hr: expected 20 hopping lines in '.../graphene/model/graphene_hr.dat'
stage5_bases error: read_xyz: expected 2 orbital lines in '.../graphene.xyz', found 1
```

The gate therefore writes, in its temporary directory, a dense copy of `_hr.dat` (every R
with all W x W entries, zeros filled, ndegen 1) and an `.xyz` with a comment line after the
orbital count, and points `--model-dir` there. That the copy is the same model is checked
inside the gate: lsquant's supercell CSR built from the copy equals the CSR wannier2sparse
builds from the committed originals to 0 on the square 48 x 48 graphene supercell, and to
2e-16 for every case once the cell index is renumbered (D3).

## 3. Decisions recorded

- **Disorder seed.** The CLI draws the Anderson disorder with the same `--seed` as the
  random-phase states (`make_basis(..., seed=rng_seed)`); there is no `--disorder-seed`. To
  run codex's "states seed 12345, disorder seed 7", the gate runs `stage5_bases.main` through
  a 9-line wrapper that replaces `lsquant.bases.make_anderson_disorder` by the same function
  with seed 7 (no other change of the code path).
- **Window.** `--window` is required by the CLI; `BlochBasis.window()` and
  `EigenBasis.window()` are Gershgorin of the real-space CSR, against the brief. The gate
  passes the eigenvalue window of `G9_params.txt` explicitly to all three bases ([-3.3, 3.3]
  eV at W = 0; at W = 1 the half-width (3 + W/2) x 1.10 = 3.85 eV, broadening
  pi 3.85 / 256 = 47.2466 meV).
- **Bloch vs eigen at 1e-8.** A random-phase Bloch curve against the exact eigen curve
  differs by the stochastic trace error (1.5e-2 to 4.9e-2), so the 1e-8 row compares
  **Bloch with exact trace** (`--trace exact`, BlockLocalStates) against eigen; the random
  Bloch vs eigen number is reported as INFO.
- **Relative difference.** Gated quantity: `nrel(a, b) = max|a - b| / max|b|` (sup-norm
  relative). The pointwise `prel` of `scripts/golden_diff.py` is printed next to it; it is
  1 wherever the reference is 0 (odd moments of a symmetric spectrum, the Haldane gap), so it
  cannot carry a 1e-9 threshold.
- **Curves at full precision.** The written `.dat` curves have 6 significant digits (the
  legacy text format), too few for 1e-8; `figures/stage_5_table.py` recomputes every curve
  from the `.chebmom` file with `lsquant.reconstruct.dos_curve` / `kubo_greenwood` at the
  run's broadening, and checks the written `.dat` against it (2.4e-6 to 3.6e-6, the 6-digit
  rounding).
- **G9 comparison.** G9 is interpolated linearly onto the KPM energies (2.5 meV grid). DOS:
  lsquant's extensive curve / N_k. sigma_xx: lsquant curve x 2 pi / (A_hr N_k), the analytic
  ratio of the two definitions (legacy Kubo-Greenwood = pi Tr[v delta v delta];
  hr_exactdiag = (2 pi^2 / A_hr) (1/N_k) sum |v_nm|^2 g g); the least-squares ratio is
  printed (0.77 graphene, 0.80 Haldane: the kernel effect). Masks as codex: graphene
  |E| < 2.5, Haldane 0.5 < |E| < 2.5.
- **Thread counts.** All stage-5 runs at 16 and 1 threads; real-path files of the two runs
  are compared with `cmp` and counted; Bloch and eigen files are compared and reported.
- **Convention rows counted.** The FFTW `B` against the dense `B` (implemented sign), `B^+ B`,
  lsquant CSR vs wannier2sparse CSR, `B H_w2s B^+ = H_k` and `B VX_w2s B^+ = v_x` (brief
  section 3) are PASS/FAIL rows of the gate; the brief's sign and the renumbered-CSR
  comparisons are INFO (localisation).
- The stage-5 W = 1 eigen row is skipped as codex section 2.4 directs.

## 4. Gate: `OMP_NUM_THREADS=16 S5_KEEP=~/s5work/s5out bash scripts/check_goldens.sh 5 ~/lsquant-unified-build5`

Final run log: `figures/out/stage_5/check_goldens_stage5.log` (`== 149 files compared;
failures: 37`, `GATE FAIL (stage 5)`, wall 53 min 34 s); the two tables as text:
`table_T16.txt`, `table_T1.txt`. The first
run (before the convention rows were counted) is `check_goldens_stage5_run1.log`
(137 compared, 32 failures = the 32 G9 rows; every other number identical to the final run).

### 4.1 Stages 1 to 4 (the real-space path did not move)

| Stage | Files | Result |
|---|---|---|
| 1 (G1, G2 at T = 1, 8, 64; batched G2; random-phase determinism) | 11 | 11 IDENTICAL |
| 2 (G1, G3 from Python, T = 1, 16) | 6 | 6 IDENTICAL |
| 3 (G4, G5 moments and curves, T = 16, 1) | 24 | 24 IDENTICAL |
| 4 (G6 moments and 8 curves, G7 moments, T = 16, 1) | 24 | 24 IDENTICAL |
| **total** | **65** | **65 IDENTICAL** |

Greps: `basis|std::cout|getenv` in `lsquant/core`: empty; `delta|DOS|Kubo|basis|std::cout|getenv`
(code lines) in core and io: empty; `compat/`: empty; `time|velocity|spin|projector|bessel|gsl`
(whole words, code lines): empty; stage 5 `fftw` in `lsquant/core/chebyshev` and
`lsquant/io`: empty.

### 4.2 FFTW convention (verified)

```
PASS  CONV graphene 6x4x1 dense B = exp(-i k.(R+tau))/sqrt(Nk) (implemented sign, gated): RealToBloch vs B x 3.589e-15, BlochToReal vs B^+ x 3.309e-15 (tol 1e-12)
INFO  CONV graphene 6x4x1 dense B = exp(+i k.(R+tau))/sqrt(Nk) (brief.md sign, reported): RealToBloch vs B x 1.523e+00, BlochToReal vs B^+ x 1.455e+00
FAIL  CONV graphene 6x4x1 lsquant supercell CSR vs wannier2sparse CSR: H 1.174e+00, VX 9.140e-01 (tol 1e-12)
INFO  CONV graphene 6x4x1 localisation: lsquant CSR with the cell index renumbered n1-fastest (n1 + N1 n2) vs wannier2sparse CSR: H 1.182e-16, VX 1.406e-16
PASS  CONV graphene 48x48x1 max|B^+ B x - x|/max|x| = 3.883e-16 (tol 1e-12)
PASS  CONV graphene 48x48x1 lsquant supercell CSR vs wannier2sparse CSR: H 0.000e+00, VX 0.000e+00 (tol 1e-12)
INFO  CONV graphene 48x48x1 localisation: lsquant CSR with the cell index renumbered n1-fastest (n1 + N1 n2) vs wannier2sparse CSR: H 1.801e-16, VX 1.813e-16
PASS  CONV graphene 48x48x1 B H_w2s B^+ x vs H_k x: 7.058e-16 (tol 1e-10)
PASS  CONV graphene 48x48x1 B VX_w2s B^+ x vs v_x(k) x: 6.652e-16 (tol 1e-10, brief.md section 3)
PASS  CONV haldane 6x4x1 dense B = exp(-i k.(R+tau))/sqrt(Nk) (implemented sign, gated): RealToBloch vs B x 3.589e-15, BlochToReal vs B^+ x 3.309e-15 (tol 1e-12)
INFO  CONV haldane 6x4x1 dense B = exp(+i k.(R+tau))/sqrt(Nk) (brief.md sign, reported): RealToBloch vs B x 1.523e+00, BlochToReal vs B^+ x 1.455e+00
FAIL  CONV haldane 6x4x1 lsquant supercell CSR vs wannier2sparse CSR: H 1.142e+00, VX 1.069e+00 (tol 1e-12)
INFO  CONV haldane 6x4x1 localisation: lsquant CSR with the cell index renumbered n1-fastest (n1 + N1 n2) vs wannier2sparse CSR: H 1.793e-16, VX 1.938e-16
PASS  CONV haldane 48x48x1 max|B^+ B x - x|/max|x| = 3.883e-16 (tol 1e-12)
FAIL  CONV haldane 48x48x1 lsquant supercell CSR vs wannier2sparse CSR: H 4.592e-01, VX 7.165e-01 (tol 1e-12)
INFO  CONV haldane 48x48x1 localisation: lsquant CSR with the cell index renumbered n1-fastest (n1 + N1 n2) vs wannier2sparse CSR: H 2.831e-16, VX 3.274e-16
FAIL  CONV haldane 48x48x1 B H_w2s B^+ x vs H_k x: 5.076e-01 (tol 1e-10)
FAIL  CONV haldane 48x48x1 B VX_w2s B^+ x vs v_x(k) x: 9.117e-01 (tol 1e-10, brief.md section 3)
```

(`rel = max|a - b| / max|b|` over 4 random complex vectors; the dense `B` is built in numpy
from `k_mesh`, `.uc` and `.xyz`; the 6 x 4 x 1 supercell is non-square so an axis swap would
show.)

So: **`B` (RealToBloch) = `FFTW_FORWARD` (exponent sign -1) on the (N1, N2, N3) grid with
`howmany = W`, stride W, followed by the phase `exp(-i k.tau_alpha)` and 1/sqrt(N_k); `B^+`
(BlochToReal) = the conjugate phase, then `FFTW_BACKWARD` (+1).** The matrix is
`B_(k,alpha),(R,beta) = delta_ab exp(-i k.(R + tau_alpha)) / sqrt(N_k)`, unitary. This is the
opposite sign of brief section 3 (`exp(+i k.(R + tau))`, `B^+` = FORWARD); the brief's sign
fails against the dense matrix by 1.5. The implemented sign is the one consistent with the
brief's own `H(k) = sum_R exp(+i k.(R + tau_b - tau_a)) H_ab(R)` and the Wannier90 /
wannier2sparse hopping convention `H[(R1, a), (R1 + R, b)] = H_ab(R)`: with the brief's `B`,
`B H B^+` is `H(-k)`. The developer's report says the same (section 2); the brief needs the
correction.

### 4.3 The codex section 2 table (T = 16; the T = 1 table is identical in every gated row)

Real and Bloch: random phase, 32 states, seed 12345. Bloch(x) and eigen: exact trace.
Window [-3.3, 3.3] eV (W = 0), [-3.85, 3.85] eV (W = 1, disorder seed 7). M = 256.

| Model, word | real vs Bloch moments, W = 0 (tol 1e-9) | Bloch(x) vs eigen curve (tol 1e-8) | real vs Bloch moments, W = 1 (tol 1e-9) | vs G9: real / Bloch / Bloch(x) / eigen (tol 2e-2 DOS, 5e-2 sigma) |
|---|---|---|---|---|
| graphene DOS | **PASS** 2.47e-14 | **PASS** 1.93e-15 | **PASS** 2.11e-14 | **FAIL** 1.52e-1 / 1.52e-1 / 1.22e-1 / 1.22e-1 (max at E = +-1.00, van Hove) |
| graphene sigma_xx | **PASS** 2.43e-15 | **PASS** 2.89e-15 | **PASS** 3.87e-16 | **FAIL** 8.08e-1 / 8.08e-1 / 7.68e-1 / 7.68e-1 (E = +-1.00) |
| Haldane DOS | **PASS** 2.47e-14 | **PASS** 1.86e-15 | **PASS** 2.17e-14 | **FAIL** 3.67e-2 / 3.67e-2 / 2.85e-2 / 2.85e-2 (E = 2.00, 2.11) |
| Haldane sigma_xx | **PASS** 2.22e-15 | **PASS** 5.74e-15 | **PASS** 4.90e-16 | **FAIL** 5.32e-1 / 5.32e-1 / 4.77e-1 / 4.77e-1 (E = 2.00, 1.78) |

(Pointwise `prel` of the moment rows: 1.8e-12 to 1.3e-10, all at near-zero moments. Real vs
Bloch curves at W = 0 and 1: 2.3e-15 to 2.1e-14.)

All G9 rows FAIL; the localisation rows (INFO, same run):

| Model, word | independent exact moments vs eigen / Bloch(x) | mesh effect: G9 vs the same Gaussian on the Gamma mesh | kernel effect: eigen vs Gaussian on the Gamma mesh | eigen vs G9 without +-0.15 eV around the van Hove energy |E| = 1 |
|---|---|---|---|---|
| graphene DOS | 2.44e-14 / 2.44e-14 | 8.48e-2 | 3.28e-2 | 1.30e-1 (at E = -1.18) |
| graphene sigma_xx | 1.15e-14 / 1.11e-14 | 6.13e-1 | 2.76e-1 | 5.60e-1 (at E = +1.24) |
| Haldane DOS | 2.44e-14 / 2.44e-14 | 1.35e-2 | 2.23e-2 | 4.84e-2 (at E = +2.11) |
| Haldane sigma_xx | 9.79e-15 / 9.59e-15 | 1.87e-1 | 3.59e-1 | 4.77e-1 (at E = +1.78) |

sigma_xx units: analytic factor 2 pi / (A_hr N_k) = 2.0993e-3; least-squares factor
1.6087e-3 (graphene) and 1.6700e-3 (Haldane), ratios 0.766 and 0.796.

"Independent exact moments": the Chebyshev moments of the 48 x 48 Gamma-mesh system computed
in numpy from hr_exactdiag's `read_hr`/`opk` `H(k)` and the wannier2sparse bundle velocity
(no lsquant model code), in the `.chebmom` layout (`kpm_weighted 1`: `(2 - delta_m0)`
factors). They agree with lsquant's eigen and exact-Bloch moments to 1e-14 for both models
and both words, so the KPM pipeline computes exactly what it should for this finite system;
the G9 differences are (i) the mesh, (ii) the kernel. Note the codex mask for graphene
(|E| < 2.5) contains the van Hove points |E| = 1, and the Haldane mask (0.5 < |E| < 2.5)
contains the band edges at |E| ~ 0.78. The clean-crystal Kubo-Greenwood sigma_xx is
dominated by the intraband delta^2 term, which scales with the inverse kernel width and
with the mesh discreteness; it does not converge in N at fixed 40 meV, hence the large
sigma_xx numbers. The stochastic error of 32 random vectors on N = 4608 is 1.5e-2 (DOS) to
4.9e-2 (sigma_xx) in nrel against the exact curve.

### 4.4 Thread counts

Every stage-5 output file (moments and curves, 2 models x 2 words x 6 runs = 48 files) is
**byte-identical between OMP_NUM_THREADS = 1 and 16**: the 16 real-space files are counted
(IDENTICAL); the 16 Bloch-random, 8 Bloch-exact and 8 eigen files are reported identical.
FFTW runs single-threaded (no `fftw_init_threads`; the OpenMP loops are over the phase
multiplications and over blocks, each element written by one thread), so the Bloch
reduction order does not depend on the thread count.

### 4.5 Wall times (TCN shared, load average 75 to 105 on 64 cores during the run)

Seconds per `stage5_bases` run (process wall time, including Python start-up and model
build), final run, T = 16 / T = 1:

| Run | graphene DOS | graphene kuboxx | Haldane DOS | Haldane kuboxx |
|---|---|---|---|---|
| W = 0 real (random 32) | 2.4 / 2.2 | 31.6 / 20.4 | 3.0 / 2.8 | 33.5 / 21.7 |
| W = 0 Bloch (random 32) | 2.4 / 2.1 | 28.1 / 19.7 | 2.8 / 2.2 | 27.9 / 19.5 |
| W = 0 Bloch (exact, 2 states) | 1.9 / 1.8 | 2.8 / 2.2 | 1.9 / 1.7 | 2.3 / 2.5 |
| W = 0 eigen (exact, 2 states) | 1.9 / 1.9 | 2.7 / 2.2 | 1.9 / 1.7 | 2.9 / 2.7 |
| W = 1 real (random 32) | 2.7 / 2.5 | 28.7 / 21.2 | 2.6 / 3.2 | 34.0 / 21.8 |
| W = 1 Bloch (random 32) | 8.7 / 4.1 | 49.6 / 23.8 | 8.1 / 4.2 | 54.0 / 22.7 |

The node was shared (load average 75 to 106 on 64 cores) and every run is faster at 1
thread than at 16; the numbers compare bases, not scaling. Exact-trace Bloch and eigen runs
use W = 2 block-local states instead of 32 random vectors. With disorder each Bloch product
adds two FFTs and four OpenMP loops over 4608 elements: W = 1 Bloch is 1.6 to 1.7 times
the real-space time at 16 threads, 1.04 to 1.12 at 1 thread.

## 5. Figure: `figures/stage_5.py`

`MPLBACKEND=Agg PYTHONPATH=~/lsquant-unified-build5/lsquant_core:python python3 figures/stage_5.py --runs ~/s5work/s5out/T16 --table figures/out/stage_5/table_T16.txt`
(exit 0) -> `figures/out/stage_5/fig_stage_5.png`, `fig_stage_5.caption.txt` (log
`figure_run.log`). The script reads the stage-5 outputs the gate kept (`S5_KEEP`) and
recomputes the curves from the moment files, as the table does. What the figure shows: the
three bases lie on one another in every panel (residuals of real - Bloch at W = 1 are
1e-14); at W = 0 the 48 x 48 Gamma-mesh DOS and sigma_xx carry a regular finite-size
ripple (period ~0.1 eV) that the G9 markers on the shifted mesh do not have, and the clean
sigma_xx of G9 sits below the Jackson curves (the kernel effect of section 4.3); at W = 1
the disorder smooths the ripple and cuts sigma_xx to a third of the clean value.

Caption:

> FIG. 5. Graphene (wannier2sparse examples/gen_models.py, t = -1 eV, nearest-neighbour distance 1 A, supercell 48 x 48 x 1 = 4608 orbitals, Bloch mesh k = n/48) evaluated by lsquant.apps.stage5_bases in three bases from one model file, M = 256 Chebyshev moments, Jackson kernel. Top: DOS per cell; bottom: sigma_xx (Kubo-Greenwood) in e^2/h with the true cell area 2.5981 A^2. Left column W = 0: real space (CSR, random phase, 32 states, seed 12345, solid), Bloch (block-diagonal H(k), the same 32 vectors mapped by B, dashed), eigen (exact trace in the eigenbasis of each H(k), dotted), window [-3.30, 3.30] eV from the eigenvalues padded 10 %, broadening 40.4971 meV (= pi a / M, the Jackson resolution of M = 256, a = 3.30 eV); markers: G9, tools/hr_exactdiag.py (wannier2sparse 3f553db) dos and sigma --comp xx on the 48 x 48 mesh (i + 1/2)/48 with Gaussian sigma = 40.4971 meV; residual basis - G9 below. Right column W = 1 eV Anderson disorder (uniform in [-0.5, 0.5] eV, numpy default_rng seed 7): real space (CSR + diag V) and Bloch (H(k) + B V B^+ through FFTW), window [-3.85, 3.85] eV, broadening 47.2466 meV; clean G9 markers for reference; residual real - Bloch below (eigen skipped: no per-k eigenbasis with disorder). Max |basis - G9| / max |G9| on |E| < 2.5 eV, W = 0: DOS real 1.52e-01, Bloch 1.52e-01, eigen 1.22e-01; sigma_xx real 8.08e-01, Bloch 8.08e-01, eigen 7.68e-01. W = 1: max |real - Bloch| / max |Bloch| on the curve: DOS 2.1e-14, sigma_xx 2.5e-15. Gate (figures/stage_5_table.py, T = 16; identical at T = 1), graphene: DOS: real vs Bloch moments 2.465e-14 (W = 0), 2.105e-14 (W = 1), tol 1e-9; Bloch (exact trace) vs eigen curve 1.925e-15, tol 1e-8; eigen moments vs exact moments from hr_exactdiag's H(k) on the same Gamma mesh 2.442e-14; sigma_xx: real vs Bloch moments 2.425e-15 (W = 0), 3.869e-16 (W = 1), tol 1e-9; Bloch (exact trace) vs eigen curve 2.890e-15, tol 1e-8; eigen moments vs exact moments from hr_exactdiag's H(k) on the same Gamma mesh 1.146e-14. The G9 rows fail the codex tolerances (2e-2 DOS, 5e-2 sigma_xx): shifted vs Gamma mesh and Jackson vs Gaussian kernel (tests.md section 4.3).

## 6. Defects (localised; not fixed)

- **D1 `read_hr` accepts only the dense Wannier90 layout.**
  `python/lsquant/models/wannier.py::read_hr` requires `nrpts * W^2` hopping lines and reads
  them by position (`base = idx + ir * block_size`), assigning `ndegen[ir]` by block
  position. gen_models.py (and any hand-written tight-binding `_hr.dat`, which wannier2sparse
  and hr_exactdiag both accept) lists only the non-zero entries, in any order:
  `read_hr: expected 20 hopping lines` for graphene (6 present). Even with the right count, a
  sparse file would get the wrong `ndegen` per entry. The brief names exactly these files
  (gen_models graphene) as the self-check and gate input.
- **D2 `read_xyz` requires a comment line after the orbital count** (`raw_lines[2:]`);
  wannier2sparse's `.xyz` (gen_models, `examples/02_graphene/graphene.xyz` style) has none:
  `read_xyz: expected 2 orbital lines ..., found 1`.
- **D3 cell index of the supercell and of `B`: `n1 N2 + n2` (lsquant) vs `n1 + N1 n2`
  (wannier2sparse).** Section 4.2 rows: lsquant CSR vs wannier2sparse CSR 1.2 (graphene
  6 x 4), 1.1 (Haldane 6 x 4), 0.46 / 0.72 (Haldane 48 x 48, H / VX);
  `B H_w2s B^+` vs `H_k` 0.51 and `B VX_w2s B^+` vs `v_x` 0.91 (Haldane 48 x 48; brief
  section 3 asks 1e-10). Renumbering lsquant's cells n1-fastest gives 1e-16 everywhere, and a
  numpy supercell built with the Wannier90 convention and cell index `n1 + N1 n2` equals the
  wannier2sparse CSR to 2e-16 (graphene and Haldane, 6 x 4), while with `n1 N2 + n2` it equals
  lsquant's CSR to 2e-16. So the hopping convention and conjugation are right on both sides;
  only the cell numbering differs. Brief section 3 prescribes `iR = n1 N2 N3 + n2 N3 + n3`,
  which is not wannier2sparse's numbering, while the brief also says the real-space basis is
  the CSR wannier2sparse writes. `report.md` section 2.1 states wannier2sparse is C row-major;
  it is not. Graphene 48 x 48 passes only because the mirror n1 <-> n2 is a symmetry of the
  model and of the square supercell; the self-check passes because it never reads a
  wannier2sparse CSR. Consequence: `B` cannot be applied to a wannier2sparse operator (the
  brief's `B^+ V B` with a wannier2sparse disorder or spin operator, any later stage that
  mixes the bases through `B`) without permuting the cells. Which side changes (the brief's
  index or the implementation) is Fable's decision.
- **D4 `stage5_bases` CLI differs from brief section 6.** `--model-dir` instead of `--model`;
  `--disorder-w` instead of `--disorder`; no `--disorder-seed` (the disorder uses the state
  seed: codex 2.4 cannot be run with the CLI alone, see section 3); `--window` is mandatory
  (brief: from the eigenvalues); `--seed` is `action="append"` with a heuristic that treats
  any token naming an existing `<tok>_hr.dat` as the model; the output is the new
  `.chebmom` format, not the legacy layouts (`.chebmom1D`/`2D`), the curve names gain
  `_<basis>_` and `mean_` (legacy: `mean...`), and no `<stem>.basis.json` is written.
- **D5 EigenBasis does not evaluate words exactly (brief sections 0, 3).** It builds
  `DiagonalOperator(eps)` and the rotated operators and runs the same Chebyshev engine with
  BlockLocalStates; there is no evaluation with the sympy `SpectralFunction` on the
  eigenvalues, and so no report of the difference between the table form and `f(eps; E)` (the
  check of AGENTS.md section 6 that the coefficient tables reproduce `f(eps; E)`). Its
  moments are exact (section 4.3), but the "exact concatenator" is absent.

## 7. Smells

- **A basis name in the engine:** none. `grep -rn basis lsquant/core`: empty;
  `grep -rni "bloch|eigen|wannier" lsquant/core/chebyshev`: one comment ("eigenvalue
  spectrum", `rescaled_operator.hpp:43`). BlochTransform, BlockDiagonalOperator,
  SandwichOperator, BlockLocalStates, TransformedStates live in `core/algebra`.
- **A block format file:** none. Blocks travel as one `(N_k, W, W)` complex array
  (`BlockDiagonalOperator(blocks)`), no reader or writer of blocks; no `fstream` in the new
  files.
- **Gershgorin on a Bloch run:** yes, by default. `BlochBasis.window()` and
  `EigenBasis.window()` return the Gershgorin window of the real-space CSR (docstring
  "matching RealBasis"), against brief section 4 ("never from Gershgorin; say so in the
  docstring"); `stage5_bases --window gershgorin` on `--basis bloch` uses it. There is no
  `EigenBasis.window(pad)` from the eigenvalues. The gate never used it (explicit window).
- **FFTW included from `core/chebyshev`:** no (`grep -rni fftw lsquant/core/chebyshev
  lsquant/io`: empty; FFTW only in `bloch_transform.cpp`).
- **The engine changed:** `lsquant/core/chebyshev/rescaled_operator.{hpp,cpp}` gained a
  constructor from a general `LinearOperator` that applies `scale H + shift` on the fly (85
  lines); brief section 2 says the engine does not change. `MomentEngine`, `Plan`,
  `ChebyshevRecursion` are untouched, and G1 to G8 are byte-identical, so the CSR path keeps
  the precomputed rescaled copy; the on-the-fly path is the one every composite operator
  (disorder sum, block operators) takes. `python/lsquant/moments.py`: 1 line.
- **Non-reentrant operators:** `SandwichOperator` keeps `mutable Vector t1_, t2_` and
  `BlochTransform` one FFTW workspace (in/out buffers, plan) per instance; `Multiply` is not
  thread-safe if the engine ever applies the same operator to several states concurrently.
  Fine today (results equal at 1 and 16 threads).
- **The real-space basis is not wannier2sparse's CSR:** `RealBasis` rebuilds the supercell in
  Python (`models/wannier.py::supercell_csr`); see D3.
- **Self-check scope** (section 2.1): own model writer, square 8 x 8 and 16 x 16 supercells,
  no dense-B check, no wannier2sparse CSR, disorder seed = state seed.
- **Output precision:** the curves are written with 6 significant digits (legacy format);
  any 1e-8 comparison has to go through the moments.
- **report.md** says "All checks (a)-(f) and the three-basis runs pass within machine
  precision" and "Open questions: none"; it does not mention that the gen_models files cannot
  be read, nor that the FFTW sign is the opposite of the brief's (it gives the sign but not
  the conflict).
- **Oracle (wannier2sparse, not LSQUANT):** `hr_exactdiag.py sigma` uses `|a1_x a2_y|` as the
  cell area (half the honeycomb cell); its mesh is shifted; `gen_models.py` is not
  parameterised by N (codex "gen_models.py, N = 48": N is the supercell, set by the gate).
