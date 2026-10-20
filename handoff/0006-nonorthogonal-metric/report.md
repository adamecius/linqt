# Handoff 0006 Report — Non-Orthogonal Bases via `MetricPolicy`

## 1. Files Created and Modified

### Created
- `lsquant/core/algebra/metric_policy.hpp` — Abstract `MetricPolicy` interface (`IsOrthogonal`, `Overlap`, `InverseOverlap`, `InverseSqrtOverlap`) and concrete `OrthogonalMetric` and `NonOrthogonalMetric` implementations.
- `lsquant/core/algebra/metric_policy.cpp` — Implementation of `OrthogonalMetric` (`IdentityOperator` views) and `NonOrthogonalMetric` (`SparseMatrixOperator` view of `S` plus non-owning references to `S^-1` and `S^-1/2` appliers with dimension and Hermiticity validation).
- `lsquant/core/algebra/overlap_appliers.hpp` — Declarations of `ExactInverseOverlap`, `ExactInverseSqrtOverlap`, `IschoInverseOverlap`, `ChebyshevInverseOverlap`, and `ChebyshevInverseSqrtOverlap`.
- `lsquant/core/algebra/overlap_appliers.cpp` — Implementations of the dense Cholesky (`LLT`) inverse applier, dense eigensolver (`SelfAdjointEigenSolver`) inverse-square-root applier, incomplete-Cholesky-preconditioned conjugate-gradient (`ISCHO` + `PCG`) applier with explicit `L^+` upper triangular solve, and Chebyshev-Gauss quadrature series appliers for `1/x` and `x^-1/2` on `spec(S)`.
- `lsquant/apps/stage6_nonorth.cpp` — Stage-6 gate CLI evaluating `dos`, `kubo`, and `kuboxy` plans across `lowdin`, `generalised`, and `naive` routes and `exact`, `ischo`, and `chebyshev` appliers, plus `--check-appliers` residual diagnostics.
- `python/lsquant/models/nonorthogonal.py` — Pure NumPy/SciPy module providing `make_local_transform`, `rotate_basis`, `read_csr_text`, `write_csr_text`, and `hermitian_fold`.
- `python/lsquant/apps/make_nonorth_fixtures.py` — CLI generator producing `real` and `complex` rotated CSR fixtures and `manifest.json` files for `chain1d` (`N = 2048`) and `haldane` (`N = 512`).
- `tests/dev/stage_6_check.py` — Developer self-check script on rotated `chain1d` (`N = 2048`, `real` and `complex` variants) verifying spectrum invariance, `lowdin/exact` and `generalised/exact` DOS moments against G1, `ischo` and `chebyshev` applier residuals and moment errors across `M_S` and `tol`, and `ExactInverseOverlap` vs `ChebyshevInverseOverlap` on ten random vectors.
- `goldens/stage0/nonorth/chain1d/{real,complex}/chain1d.{HAM,S,VX}.CSR` and `manifest.json` — Controlled non-orthogonal `chain1d` fixtures (`N = 2048`, `seed = 12345`, `alpha = 0.2`).
- `goldens/stage0/nonorth/haldane/{real,complex}/haldane.{HAM,S,VX,VY}.CSR` and `manifest.json` — Controlled non-orthogonal `haldane` fixtures (`N = 512`, `seed = 12345`, `alpha = 0.2`).

### Modified
- `lsquant/core/algebra/composite_operators.hpp` — Added `LowdinOperator` (`S^-1/2 O S^-1/2`) and `MetricInverseProduct` (`S^-1 O`).
- `lsquant/core/algebra/composite_operators.cpp` — Implemented matrix-free `Multiply` and `IsHermitian` for `LowdinOperator` and `MetricInverseProduct`.
- `lsquant/core/chebyshev/rescaled_operator.hpp` — Added matrix-free `RescaledOperator(const LinearOperator&, const SpectralWindow&)` overload and `GeneralisedRescaledOperator(const SparseMatrix& H, const SparseMatrix& S, const LinearOperator& inverse, const SpectralWindow& window)`.
- `lsquant/core/chebyshev/rescaled_operator.cpp` — Implemented matrix-free `RescaledOperator` evaluation (`y = a * (scale * H * x + shift * x) + b * y`) and `GeneralisedRescaledOperator` (`S^-1 (scale * H' + shift * S')`).
- `lsquant/CMakeLists.txt` — Registered `metric_policy.cpp` and `overlap_appliers.cpp` in `lsquant_core` and added the `stage6_nonorth` executable target.

---

## 2. Fixture Manifests (`goldens/stage0/nonorth/`)

Generated on TCN (`jh-x`) with `python3 python/lsquant/apps/make_nonorth_fixtures.py` (`seed = 12345`, `alpha = 0.2`):

| Fixture | `N` | `seed` | `alpha` | `lambda_min(S')` | `lambda_max(S')` | `cond(S')` | `nnz(S')` | `nnz(H')` | `max \|eigvalsh(H', S') - eigvalsh(H)\|` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `chain1d/real` | 2048 | 12345 | 0.2 | `0.7032008650483048` | `1.352076700392037` | `1.922746070995176` | 10240 | 14336 | `7.993605777301127e-15` |
| `chain1d/complex` | 2048 | 12345 | 0.2 | `0.7246054670422337` | `1.3259502278293853` | `1.829892663164384` | 10240 | 14336 | `6.217248937900877e-15` |
| `haldane/real` | 512 | 12345 | 0.2 | `0.7203427769338232` | `1.3370967064526011` | `1.8561950633336304` | 15872 | 32768 | `1.9539925233402755e-14` |
| `haldane/complex` | 512 | 12345 | 0.2 | `0.7367400475481392` | `1.324598043711703` | `1.7979177976274634` | 15872 | 32768 | `1.021405182655144e-14` |

---

## 3. Orthogonal G5 (`kuboxy`) Operator-Mapping Verification

In legacy `src/inline_compute-kpm-nonEqOp.cpp` and `src/compute.cpp`, `inline_compute-kpm-nonEqOp haldane VX VY 128 exact` assigns `OP[1] = VX` and `OP[2] = VY` and calls `CorrelationExpansionMoments(OP[1], OP[2], chebMoms, gen)`. In `src/chebyshev_solver.cpp` (`CorrelationExpansionMoments(SparseMatrixType &OPL, SparseMatrixType &OPR, ...)`), `OPL = VX` acts on the initial bra state (`chevVecL`: `SetInitVectors(OPL, PhiL); IterateAll()`) while `OPR = VY` acts after Chebyshev expansion on the ket side (`chevVecR`: `SetInitVectors(PhiR); IterateAll(); Multiply(OPR)`).

Accordingly, in `MomentEngine` bytecode:
- **Orthogonal / Naive and Löwdin `kuboxy` plan**:
  - `ket`: `Expand(M), ApplyOperator(VY)`
  - `MirrorToBra`
  - `bra`: `ApplyOperator(VX), Expand(M)`
  - `Trace(1)` + `ApplyKpmWeights(mu, numStates, true)`
- **Generalised `kuboxy` plan (Eq. 2c of `theory.md`)**:
  - `ket`: `Expand(M), ApplyOperator(MetricInverseProduct(S^-1, VY'))`
  - `MirrorToBra`
  - `bra`: `ApplyOperator(S^-1), ApplyOperator(MetricInverseProduct(S^-1, VX')), Expand(M), ApplyOperator(S')`
  - `Trace(1)` + `ApplyKpmWeights(mu, numStates, true)`

Verification on the orthogonal Haldane fixture (`goldens/stage0/haldane/operators`):
```text
$ ./build/lsquant_core/stage6_nonorth --fixture goldens/stage0/haldane/operators --label haldane --window -4.5 4.5 --route naive --applier exact --trace exact --moments 128 --word kuboxy --out /tmp/stage6_out/g5_orth
stage6_nonorth: route=naive applier=exact word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/g5_orth/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
$ cmp /tmp/stage6_out/g5_orth/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D goldens/stage0/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
(exit 0, byte-for-byte identical)
```

---

## 4. Output of `tests/dev/stage_6_check.py`

Executed on TCN (`jh-x`, `~/lsquant-agi-0006-build`):

```text
=== Stage 6 Self-Check: Rotated chain1d (N = 2048) ===
[1] real: max |eigvalsh(H', S') - eigvalsh(H)| = 7.549517e-15
[4] real: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = 2.965455e-13
[2] real route=lowdin      applier=exact     applier_res=2.446625e-14 max|mu - G1|=6.296135e-13 max|mu - delta_m0|=6.296135e-13
[2] real route=generalised applier=exact     applier_res=2.818816e-16 max|mu - G1|=4.381506e-16 max|mu - delta_m0|=4.381506e-16
[3] real route=lowdin      applier=chebyshev (M_S=12) applier_res=8.689113e-11 max|mu - delta_m0|=9.169660e-10
[3] real route=generalised applier=chebyshev (M_S=12) applier_res=2.831615e-10 max|mu - delta_m0|=2.790968e-09
[3] real route=lowdin      applier=chebyshev (M_S=60) applier_res=2.622059e-12 max|mu - delta_m0|=1.104298e-11
[3] real route=generalised applier=chebyshev (M_S=60) applier_res=3.026105e-13 max|mu - delta_m0|=1.532508e-12
[3] real route=generalised applier=ischo     (tol=1e-06) applier_res=1.067316e-09 max|mu - delta_m0|=6.971935e-09
[3] real route=generalised applier=ischo     (tol=1e-13) applier_res=3.825733e-14 max|mu - delta_m0|=4.564573e-15
[1] complex: max |eigvalsh(H', S') - eigvalsh(H)| = 1.798561e-14
[4] complex: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = 4.699590e-13
[2] complex route=lowdin      applier=exact     applier_res=2.555632e-14 max|mu - G1|=5.587664e-13 max|mu - delta_m0|=5.587664e-13
[2] complex route=generalised applier=exact     applier_res=2.954083e-16 max|mu - G1|=1.475134e-16 max|mu - delta_m0|=1.475134e-16
[3] complex route=lowdin      applier=chebyshev (M_S=60) applier_res=8.102436e-13 max|mu - delta_m0|=3.776543e-12
[3] complex route=generalised applier=chebyshev (M_S=60) applier_res=4.815453e-13 max|mu - delta_m0|=1.744433e-12
[3] complex route=generalised applier=ischo     (tol=1e-13) applier_res=8.068502e-17 max|mu - delta_m0|=2.915249e-15
=== Stage 6 Self-Check PASSED ===
```

---

## 5. `stage6_nonorth` Runs Across Routes, Appliers, and Models (`/tmp/stage6_out`)

One-line stdout summary from each `stage6_nonorth` run on TCN (`jh-x`):

```text
stage6_nonorth: route=lowdin applier=exact word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_lowdin_exact/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
stage6_nonorth: route=lowdin applier=chebyshev word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_lowdin_chebyshev/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
stage6_nonorth: route=generalised applier=exact word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_generalised_exact/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
stage6_nonorth: route=generalised applier=ischo word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_generalised_ischo/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
stage6_nonorth: route=generalised applier=chebyshev word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_generalised_chebyshev/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
stage6_nonorth: route=naive applier=exact word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_naive_exact/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
stage6_nonorth: route=lowdin applier=exact word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_lowdin_exact/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
stage6_nonorth: route=lowdin applier=chebyshev word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_lowdin_chebyshev/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
stage6_nonorth: route=generalised applier=exact word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_generalised_exact/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
stage6_nonorth: route=generalised applier=ischo word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_generalised_ischo/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
stage6_nonorth: route=generalised applier=chebyshev word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_generalised_chebyshev/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
stage6_nonorth: route=naive applier=exact word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_naive_exact/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
```

(`--route lowdin --applier ischo` is rejected with `std::invalid_argument` because `IschoInverseOverlap` provides `S^-1`, whereas the Löwdin sandwich requires `S^-1/2`.)

---

## 6. Deviations from the Brief

1. **`lsquant/api/bindings.cpp` omitted on `handoff/0006-agi`**:
   `handoff/0006-agi` is cut directly from `origin/cpp-unified` (Stage 1) and runs in parallel with Stage 2 (which creates `lsquant/api/bindings.cpp` in its own branch). As updated in commit `eb69eb8` ("handoff 0006: fixture generator independent of the stage-2 binding"), `make_nonorth_fixtures.py` and `stage_6_check.py` use plain Python (`numpy`/`scipy`) and the `stage6_nonorth` CLI executable without pybind11 bindings, avoiding a merge conflict with Stage 2.
2. **Sub-`1e-18` pruning of the dense inverse in `ExactInverseOverlap` and `ExactInverseSqrtOverlap`**:
   Both classes compute the full dense inverse (`Eigen::LLT<Eigen::MatrixXcd>`) and dense inverse square root (`Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>`) in their constructors as specified, and then store entries with `|A_{ij}| > 1e-18` (below IEEE-754 double precision rounding noise `2.22e-16`) in an internal compressed `SparseMatrix` so that `262,144` exact-trace applications on `N = 2048` execute via the row-parallel CSR kernel in seconds while preserving `< 1e-12` agreement with G1 (`4.38e-16` on `generalised/exact`, `5.59e-13` on `lowdin/exact`).

---

## 7. Open Questions for Fable

None. `MomentEngine` and `Plan` required no modifications to execute either the Löwdin route or the corrected generalised route (`bra: ApplyOperator(S^-1), ApplyOperator(S^-1 V_bra'), Expand(M), ApplyOperator(S')`).
