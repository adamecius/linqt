# Handoff 0006 Report — Non-Orthogonal Bases via `MetricPolicy`

## 1. Files Created and Modified

### Created
- `lsquant/core/algebra/metric_policy.hpp` — Abstract `MetricPolicy` interface (`IsOrthogonal`, `Overlap`, `InverseOverlap`, `InverseSqrtOverlap`) and concrete `OrthogonalMetric` and `NonOrthogonalMetric` implementations.
- `lsquant/core/algebra/metric_policy.cpp` — Implementation of `OrthogonalMetric` (`IdentityOperator` views) and `NonOrthogonalMetric` (`SparseMatrixOperator` view of `S` plus non-owning references to `S^-1` and `S^-1/2` appliers with dimension and Hermiticity validation).
- `lsquant/core/algebra/overlap_appliers.hpp` — Declarations of `ExactInverseOverlap`, `ExactInverseSqrtOverlap`, `IschoInverseOverlap`, `ChebyshevInverseOverlap`, and `ChebyshevInverseSqrtOverlap`.
- `lsquant/core/algebra/overlap_appliers.cpp` — Implementations of the dense Cholesky (`LLT`) inverse applier, dense eigensolver (`SelfAdjointEigenSolver`) inverse-square-root applier, incomplete-Cholesky-preconditioned conjugate-gradient (`ISCHO` + `PCG`) applier with explicit `L^+` upper triangular solve, and Chebyshev-Gauss quadrature series appliers for `1/x` and `x^-1/2` on `spec(S)`.
- `lsquant/apps/stage6_nonorth.cpp` — Stage-6 gate CLI evaluating `dos`, `kubo`, and `kuboxy` plans across `lowdin`, `generalised`, and `naive` routes and `exact`, `ischo`, and `chebyshev` appliers, plus `--check-appliers` residual diagnostics.
- `python/lsquant/models/__init__.py` — Package marker for `lsquant.models`.
- `python/lsquant/models/nonorthogonal.py` — Pure NumPy/SciPy module providing `make_local_transform`, `rotate_basis`, `read_csr_text`, `write_csr_text`, and `hermitian_fold`.
- `python/lsquant/apps/make_nonorth_fixtures.py` — CLI generator producing `real` and `complex` rotated CSR fixtures and `manifest.json` files for `chain1d` (`N = 2048`) and `haldane` (`N = 512`).
- `tests/dev/stage_6_check.py` — Developer self-check script on rotated `chain1d` (`N = 2048`, `real` and `complex` variants) verifying spectrum invariance, `lowdin/exact` and `generalised/exact` DOS moments against G1, `ischo` and `chebyshev` applier residuals and moment errors across `M_S` and `tol`, and `ExactInverseOverlap` vs `ChebyshevInverseOverlap` on ten random vectors.
- `goldens/stage0/nonorth/chain1d/{real,complex}/chain1d.{HAM,S,VX}.CSR` and `manifest.json` — Controlled non-orthogonal `chain1d` fixtures (`N = 2048`, `seed = 12345`, `alpha = 0.55`, `cond(S')` in `[5, 10]`, with `window` and `.CSR` `sha256` digests).
- `goldens/stage0/nonorth/haldane/{real,complex}/haldane.{HAM,S,VX,VY}.CSR` and `manifest.json` — Controlled non-orthogonal `haldane` fixtures (`N = 512`, `seed = 12345`, `alpha = 0.55`, `cond(S')` in `[5, 10]`, with `window` and `.CSR` `sha256` digests).

### Modified
- `lsquant/core/algebra/composite_operators.hpp` — Added `LowdinOperator` (`S^-1/2 O S^-1/2`) and `MetricInverseProduct` (`S^-1 O`).
- `lsquant/core/algebra/composite_operators.cpp` — Implemented matrix-free `Multiply` and `IsHermitian` for `LowdinOperator` and `MetricInverseProduct`.
- `lsquant/core/chebyshev/rescaled_operator.hpp` — Added matrix-free `RescaledOperator(const LinearOperator&, const SpectralWindow&)` overload and `GeneralisedRescaledOperator(const SparseMatrix& H, const SparseMatrix& S, const LinearOperator& inverse, const SpectralWindow& window)`.
- `lsquant/core/chebyshev/rescaled_operator.cpp` — Implemented matrix-free `RescaledOperator` evaluation (`y = a * (scale * H * x + shift * x) + b * y`) and `GeneralisedRescaledOperator` (`S^-1 (scale * H' + shift * S')`).
- `lsquant/api/bindings.cpp` — Bound `MetricPolicy`, `OrthogonalMetric`, `NonOrthogonalMetric`, `ExactInverseOverlap`, `ExactInverseSqrtOverlap`, `IschoInverseOverlap` (`last_iterations`, `last_error`), `ChebyshevInverseOverlap` (`effective_order`), `ChebyshevInverseSqrtOverlap` (`effective_order`), `LowdinOperator`, `MetricInverseProduct`, `GeneralisedRescaledOperator`, and the `LinearOperator` overload of `RescaledOperator` with `py::keep_alive` lifetime annotations.
- `python/lsquant/_core.pyi` — Added type stubs for the Stage-6 `MetricPolicy` classes, overlap appliers, and non-orthogonal composite/rescaled operators.
- `lsquant/CMakeLists.txt` — Registered `metric_policy.cpp` and `overlap_appliers.cpp` in `lsquant_core` and added the `stage6_nonorth` executable target.

---

## 2. Fixture Manifests (`goldens/stage0/nonorth/`)

Generated on TCN (`jh-x`) with `python3 python/lsquant/apps/make_nonorth_fixtures.py` (`seed = 12345`, `alpha = 0.55`):

| Fixture | `N` | `seed` | `alpha` | `window` | `lambda_min(S')` | `lambda_max(S')` | `cond(S')` | `nnz(S')` | `nnz(H')` | `max \|eigvalsh(H', S') - eigvalsh(H)\|` |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| `chain1d/real` | 2048 | 12345 | 0.55 | `[-2.0, 2.0]` | `0.3137719327555207` | `2.1035539282020235` | `6.704085702404217` | 10240 | 14336 | `7.327471962526033e-15` |
| `chain1d/complex` | 2048 | 12345 | 0.55 | `[-2.0, 2.0]` | `0.36043374446001014` | `2.020634828262849` | `5.606120013235988` | 10240 | 14336 | `9.547918011776346e-15` |
| `haldane/real` | 512 | 12345 | 0.55 | `[-4.5, 4.5]` | `0.35741617850047247` | `2.069165908085089` | `5.789234042975348` | 15872 | 32768 | `1.8207657603852567e-14` |
| `haldane/complex` | 512 | 12345 | 0.55 | `[-4.5, 4.5]` | `0.38717798606128334` | `2.020663701889391` | `5.218952974174405` | 15872 | 32768 | `1.3766765505351941e-14` |

### Fixture `.CSR` SHA-256 Digests
- `chain1d/real`:
  - `chain1d.S.CSR`: `146d1b4abd9e540b140542e7cb8846c8165a8da8171e491ef62e0cb8d1b3a871`
  - `chain1d.HAM.CSR`: `9dd033d003a26030023b7850a2dfe7549e6126eb1e5ed9db4839ec4c81a641e9`
  - `chain1d.VX.CSR`: `333a74ea8533ce898ca259ab8f8692ba88b6eefc4917be271cdbe2904ad957da`
- `chain1d/complex`:
  - `chain1d.S.CSR`: `9fdf01b77af35045515ef4184369f80426f9eeee1240cbbdd592a5e67e36caed`
  - `chain1d.HAM.CSR`: `e6a2f29b3bcecd1266d488a022b45da073e973c7cfacc9325d1597eea597db31`
  - `chain1d.VX.CSR`: `9c232ef939e11f27ca84d0411168f189c6dd000115921caeb7ce6047811b2b02`
- `haldane/real`:
  - `haldane.S.CSR`: `a8e81e919467fb6501325dbdba7e93fbdaf54c95d98b2e54dc7746d2a9e24fe2`
  - `haldane.HAM.CSR`: `9357600923caeefaeff7d458c7e0c4c4e7fd85912b8556f52d60af1f995d2fc0`
  - `haldane.VX.CSR`: `531b58a1fb00422cf3ae7fe2028f33065fdbf55f65e62977454a3d0efa5d17fd`
  - `haldane.VY.CSR`: `0a16ddb3a5b5e440a51b94c52f72fec6d4a92bc6dd8104604c94bd569494b8a7`
- `haldane/complex`:
  - `haldane.S.CSR`: `557f35d954dfc88169d941c0378c48566f705b7c38bbff37cf1bb40c43834197`
  - `haldane.HAM.CSR`: `224f23bfb1c1d2d8e30bbd5badce20d6b907cc94125566d3ed28e903adca1e6f`
  - `haldane.VX.CSR`: `d658cda1d87869b97264191fcbf191a6dc7a9af089b3dd5aa97510d24b7a084a`
  - `haldane.VY.CSR`: `7db4c32cce9bce8a71dcf41937094496f141647aaf1e4c6782c2cde82b6c6976`

---

## 3. Orthogonal G5 (`kuboxy`) Operator-Mapping Verification and Stage 1–3 Self-Checks

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

Verification on the orthogonal Haldane fixture (`goldens/stage0/haldane/operators`) after merging `origin/cpp-unified` (`552de6e`):
```text
$ ./build/lsquant_core/stage6_nonorth --fixture goldens/stage0/haldane/operators --label haldane --window -4.5 4.5 --route naive --applier exact --trace exact --moments 128 --word kuboxy --out /tmp/stage6_out/g5_orth
stage6_nonorth: route=naive applier=exact word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/g5_orth/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D order=60 tol=1.000000e-13 max_iterations=200 residual=0.000000e+00
$ cmp /tmp/stage6_out/g5_orth/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D goldens/stage0/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
(exit 0, byte-for-byte identical)
```

Stage-1, Stage-2, and Stage-3 self-checks rerun on TCN (`jh-x`) after the `origin/cpp-unified` merge and `lsquant._core` binding update:
```text
$ ./build/lsquant_core/stage_1_check
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES

$ PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_2_check.py
(a) dos_plan array_equal vs G1:  True
(a) kubo_plan array_equal vs G2: True
(b) DeltaFunction.table vs legacy_delta_table max_abs_diff = 2.220446e-16
(c) dos_curve vs G3 differing lines: 0 (total lines: 1590)

$ PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_3_check.py
(0) SparseMatrix.from_csr (scipy) dos array_equal vs ReadCsrText: True
(a) Haldane VX-VY array_equal vs G5: True
(a) Haldane VX-VX array_equal vs G5: True
(b) Graphene Gershgorin BandWidth = 15.0366847652085, array_equal vs G4: True
(c) GreenRetarded.table vs compat max_abs_diff           = 1.148853e-14
(c) GreenRetardedDerivative.table vs compat max_abs_diff = 2.702165e-13 (max_rel_diff = 7.173405e-15)
(d) G5 Bastin (10 meV)           differing lines: 0 (total lines: 3839)
(d) G5 Sea (10 meV)              differing lines: 0 (total lines: 3839)
(d) G5 Surf (10 meV)             differing lines: 0 (total lines: 3839)
(d) G5 Greenwood (10 meV)        differing lines: 0 (total lines: 127)
(d) G4 Bastin (64 meV)           differing lines: 0 (total lines: 1919)
(d) G4 Sea (10 meV)              differing lines: 0 (total lines: 1919)
(d) G4 Surf (10 meV)             differing lines: 0 (total lines: 1919)
(d) G4 Bastin FFTgrid (10 meV)   differing lines: 0 (total lines: 63)
(d) G4 Greenwood FFTgrid (10 meV) differing lines: 0 (total lines: 63)
```

---

## 4. Output of `tests/dev/stage_6_check.py`

Executed on TCN (`jh-x`, `~/lsquant-agi-0006-build`):

```text
=== Stage 6 Self-Check: Rotated chain1d (N = 2048) ===
[1] real: max |eigvalsh(H', S') - eigvalsh(H)| = 2.065015e-14
[4] real: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = 1.452447e-12
[2] real route=lowdin      applier=exact     applier_res=2.686186e-14 max|mu - G1|=7.956630e-13 max|mu - delta_m0|=7.956630e-13
[2] real route=generalised applier=exact     applier_res=3.918055e-16 max|mu - G1|=4.664319e-16 max|mu - delta_m0|=4.664319e-16
[3] real route=lowdin      applier=chebyshev (M_S=12, eff_order=12) applier_res=7.094662e-05 max|mu - delta_m0|=5.308428e-04
[3] real route=generalised applier=chebyshev (M_S=12, eff_order=12) applier_res=2.580411e-04 max|mu - delta_m0|=1.859680e-03
[3] real route=lowdin      applier=chebyshev (M_S=60, eff_order=38) applier_res=1.781045e-12 max|mu - delta_m0|=9.506879e-12
[3] real route=generalised applier=chebyshev (M_S=60, eff_order=41) applier_res=1.687868e-12 max|mu - delta_m0|=1.517385e-12
[3] real route=generalised applier=ischo     (tol=1e-06) applier_res=5.071275e-07 max|mu - delta_m0|=1.675476e-08
[3] real route=generalised applier=ischo     (tol=1e-13) applier_res=1.345871e-16 max|mu - delta_m0|=2.737152e-15
[1] complex: max |eigvalsh(H', S') - eigvalsh(H)| = 7.993606e-15
[4] complex: max ||S^-1_cheb v - S^-1_exact v|| / ||S^-1_exact v|| (M_S=60) = 9.664000e-13
[2] complex route=lowdin      applier=exact     applier_res=3.004640e-14 max|mu - G1|=6.211824e-13 max|mu - delta_m0|=6.211824e-13
[2] complex route=generalised applier=exact     applier_res=4.162287e-16 max|mu - G1|=1.022115e-15 max|mu - delta_m0|=1.022115e-15
[3] complex route=lowdin      applier=chebyshev (M_S=60, eff_order=33) applier_res=1.782907e-12 max|mu - delta_m0|=9.336730e-13
[3] complex route=generalised applier=chebyshev (M_S=60, eff_order=36) applier_res=1.138146e-12 max|mu - delta_m0|=2.915651e-12
[3] complex route=generalised applier=ischo     (tol=1e-13) applier_res=1.491737e-16 max|mu - delta_m0|=2.860823e-15
=== Stage 6 Self-Check PASSED ===
```

---

## 5. `stage6_nonorth` Runs Across Routes, Appliers, and Models (`/tmp/stage6_out`)

One-line stdout summary from each `stage6_nonorth` run on TCN (`jh-x`):

```text
stage6_nonorth: route=lowdin applier=exact word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_lowdin_exact/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D order=60 tol=1.000000e-13 max_iterations=200 residual=2.791140e-14
stage6_nonorth: route=lowdin applier=chebyshev word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_lowdin_chebyshev/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D order=60 effective_order=33 tol=1.000000e-13 max_iterations=200 residual=1.758820e-12
stage6_nonorth: route=generalised applier=exact word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_generalised_exact/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D order=60 tol=1.000000e-13 max_iterations=200 residual=4.014038e-16
stage6_nonorth: route=generalised applier=ischo word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_generalised_ischo/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D order=60 tol=1.000000e-13 max_iterations=200 iterations=0 residual=8.810763e-16
stage6_nonorth: route=generalised applier=chebyshev word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_generalised_chebyshev/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D order=60 effective_order=36 tol=1.000000e-13 max_iterations=200 residual=1.131409e-12
stage6_nonorth: route=naive applier=exact word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage6_out/chain1d_complex_dos_naive_exact/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D order=60 tol=1.000000e-13 max_iterations=200 residual=0.000000e+00
stage6_nonorth: route=lowdin applier=exact word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_lowdin_exact/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D order=60 tol=1.000000e-13 max_iterations=200 residual=1.523399e-14
stage6_nonorth: route=lowdin applier=chebyshev word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_lowdin_chebyshev/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D order=60 effective_order=31 tol=1.000000e-13 max_iterations=200 residual=2.213699e-12
stage6_nonorth: route=generalised applier=exact word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_generalised_exact/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D order=60 tol=1.000000e-13 max_iterations=200 residual=1.034914e-15
stage6_nonorth: route=generalised applier=ischo word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_generalised_ischo/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D order=60 tol=1.000000e-13 max_iterations=200 iterations=6 residual=8.559000e-14
stage6_nonorth: route=generalised applier=chebyshev word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_generalised_chebyshev/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D order=60 effective_order=34 tol=1.000000e-13 max_iterations=200 residual=1.136662e-12
stage6_nonorth: route=naive applier=exact word=kuboxy dim=512 states=512 (exact) M=128 -> /tmp/stage6_out/haldane_complex_kuboxy_naive_exact/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D order=60 tol=1.000000e-13 max_iterations=200 residual=0.000000e+00
```

(`--route lowdin --applier ischo` is rejected with `std::invalid_argument` because `IschoInverseOverlap` provides `S^-1`, whereas the Löwdin sandwich requires `S^-1/2`.)

---

## 6. Deviations from the Brief

1. **Stage-2 and Stage-3 merge (`origin/cpp-unified` at `552de6e`) and `lsquant._core` binding**:
   `handoff/0006-agi` now merges `origin/cpp-unified` (`552de6e`, closing stages 2 and 3) and binds `MetricPolicy`, `OrthogonalMetric`, `NonOrthogonalMetric`, `ExactInverseOverlap`, `ExactInverseSqrtOverlap`, `IschoInverseOverlap` (`last_iterations`, `last_error`), `ChebyshevInverseOverlap` (`effective_order`), `ChebyshevInverseSqrtOverlap` (`effective_order`), `LowdinOperator`, `MetricInverseProduct`, `GeneralisedRescaledOperator`, and the `LinearOperator` overload of `RescaledOperator` in `lsquant/api/bindings.cpp` and `python/lsquant/_core.pyi`. `make_nonorth_fixtures.py` and `stage_6_check.py` remain self-contained on `numpy`/`scipy` and `stage6_nonorth`.
2. **Sub-`1e-18` pruning of the dense inverse in `ExactInverseOverlap` and `ExactInverseSqrtOverlap`**:
   Both classes compute the full dense inverse (`Eigen::LLT<Eigen::MatrixXcd>`) and dense inverse square root (`Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>`) in their constructors as specified, and then store entries with `|A_{ij}| > 1e-18` (below IEEE-754 double precision rounding noise `2.22e-16`) in an internal compressed `SparseMatrix` while preserving `< 1e-12` agreement with G1.

---

## 7. Open Questions for Fable

None. `MomentEngine` and `Plan` required no modifications to execute either the Löwdin route or the corrected generalised route (`bra: ApplyOperator(S^-1), ApplyOperator(S^-1 V_bra'), Expand(M), ApplyOperator(S')`).
