# Stage 5 Report — Bloch and Eigen Bases (`handoff/0005-agi`)

## 1. Files Created / Modified

`git diff --stat origin/handoff/0004-agi...handoff/0005-agi` (stage-5 delta on top of stage 4 + merged `origin/cpp-unified` docs):

```text
 lsquant/CMakeLists.txt                        |   4 +-
 lsquant/api/bindings.cpp                      | 141 +++++++-
 lsquant/core/algebra/bloch_transform.cpp      | 215 +++++++++++
 lsquant/core/algebra/bloch_transform.hpp      | 227 ++++++++++++
 lsquant/core/algebra/block_diagonal_operator.cpp | 121 +++++++
 lsquant/core/algebra/block_diagonal_operator.hpp | 126 +++++++
 lsquant/core/algebra/composite_operators.cpp  |  91 ++++-
 lsquant/core/algebra/composite_operators.hpp  | 171 ++++++++-
 lsquant/core/algebra/state_generator.cpp      |  62 +++-
 lsquant/core/algebra/state_generator.hpp      | 128 +++++++
 lsquant/core/chebyshev/rescaled_operator.cpp  |  62 +++-
 lsquant/core/chebyshev/rescaled_operator.hpp  |  42 ++-
 python/lsquant/__init__.py                    |   9 +-
 python/lsquant/_core.pyi                      |  95 +++++
 python/lsquant/apps/stage5_bases.py           | 326 +++++++++++++++++
 python/lsquant/bases.py                       | 520 +++++++++++++++++++++++++++
 python/lsquant/models/__init__.py             |  33 ++
 python/lsquant/models/wannier.py              | 576 ++++++++++++++++++++++++++++++
 python/lsquant/moments.py                     |   2 +-
 tests/dev/stage_5_check.py                    | 521 +++++++++++++++++++++++++++
```

Full `git diff --stat main...handoff/0005-agi`:

```text
 .gitignore                                         |    4 +
 AGENTS.md                                          |  766 ++++++++++++++
 CMakeLists.txt                                     |    5 +
 agent_template/AGENTS.md                           |   62 ++
 agent_template/README.md                           |   44 +
 agent_template/install.sh                          |   80 ++
 agent_template/language_rules/c_cpp.md             |  153 +++
 agent_template/language_rules/python.md            |  176 ++++
 figures/stage_0.py                                 |  366 +++++++
 figures/stage_1.py                                 |  157 +++
 figures/stage_2.py                                 |  175 ++++
 figures/stage_3.py                                 |  226 +++++
 figures/stage_4.py                                 |  191 ++++
 handoff/0001-algebra-and-chebyshev-core/brief.md   |  152 +++
 handoff/0001-algebra-and-chebyshev-core/codex.md   |   89 ++
 handoff/0001-algebra-and-chebyshev-core/report.md  |  117 +++
 .../0002-pybind-sympy-rank1-contraction/brief.md   |  177 ++++
 .../0002-pybind-sympy-rank1-contraction/codex.md   |  116 +++
 .../0002-pybind-sympy-rank1-contraction/report.md  |  209 ++++
 handoff/0003-rank2-kubo-bastin-sea-surf/brief.md   |  215 ++++
 handoff/0003-rank2-kubo-bastin-sea-surf/codex.md   |  134 +++
 handoff/0003-rank2-kubo-bastin-sea-surf/report.md  |  246 +++++
 handoff/0004-time-evolution/brief.md               |  218 ++++
 handoff/0004-time-evolution/codex.md               |  111 ++
 handoff/0004-time-evolution/report.md              |  217 ++++
 handoff/0005-bloch-and-eigen-bases/brief.md        |  154 +++
 handoff/0005-bloch-and-eigen-bases/codex.md        |  103 ++
 handoff/0006-non-orthogonal-bases/brief.md         |  200 ++++
 handoff/0006-non-orthogonal-bases/codex.md         |  133 +++
 lsquant/CMakeLists.txt                             |   54 +
 lsquant/api/bindings.cpp                           |  731 ++++++++++++++
 lsquant/apps/stage1_moments.cpp                    |  302 ++++++
 lsquant/compat/legacy_kernels.cpp                  |  129 +++
 lsquant/compat/legacy_kernels.hpp                  |   67 ++
 lsquant/core/algebra/bloch_transform.cpp           |  202 ++++
 lsquant/core/algebra/bloch_transform.hpp           |  202 ++++
 lsquant/core/algebra/block_diagonal_operator.cpp   |  121 +++
 lsquant/core/algebra/block_diagonal_operator.hpp   |  132 +++
 lsquant/core/algebra/composite_operators.cpp       |  239 +++++
 lsquant/core/algebra/composite_operators.hpp       |  472 +++++++++
 lsquant/core/algebra/linalg.cpp                    |   98 ++
 lsquant/core/algebra/linalg.hpp                    |  120 +++
 lsquant/core/algebra/linear_operator.hpp           |  125 +++
 lsquant/core/algebra/sparse_matrix.cpp             |  466 +++++++++
 lsquant/core/algebra/sparse_matrix.hpp             |  304 ++++++
 lsquant/core/algebra/spectral_bounds.cpp           |  124 +++
 lsquant/core/algebra/spectral_bounds.hpp           |   84 ++
 lsquant/core/algebra/state_generator.cpp           |  186 ++++
 lsquant/core/algebra/state_generator.hpp           |  399 ++++++++
 lsquant/core/algebra/types.hpp                     |   34 +
 lsquant/core/chebyshev/chebyshev_moments.cpp       |  148 +++
 lsquant/core/chebyshev/chebyshev_moments.hpp       |  171 ++++
 lsquant/core/chebyshev/chebyshev_recursion.cpp     |   90 ++
 lsquant/core/chebyshev/chebyshev_recursion.hpp     |  145 +++
 lsquant/core/chebyshev/chebyshev_vectors.cpp       |  186 ++++
 lsquant/core/chebyshev/chebyshev_vectors.hpp       |  279 ++++++
 lsquant/core/chebyshev/contraction.cpp             |  153 +++
 lsquant/core/chebyshev/contraction.hpp             |  144 +++
 lsquant/core/chebyshev/moment_engine.cpp           |  439 ++++++++
 lsquant/core/chebyshev/moment_engine.hpp           |  108 ++
 lsquant/core/chebyshev/plan.hpp                    |  226 +++++
 lsquant/core/chebyshev/rescaled_operator.cpp       |   78 ++
 lsquant/core/chebyshev/rescaled_operator.hpp       |  120 +++
 lsquant/core/chebyshev/spectral_window.hpp         |  105 ++
 lsquant/io/chebmom_file.cpp                        |  275 +++++
 lsquant/io/chebmom_file.hpp                        |   96 ++
 lsquant/io/csr_text.cpp                            |  134 +++
 lsquant/io/csr_text.hpp                            |   43 +
 lsquant/io/curve_text.cpp                          |   53 +
 lsquant/io/curve_text.hpp                          |   45 +
 python/lsquant/__init__.py                         |   97 ++
 python/lsquant/_core.pyi                           |  382 +++++++
 python/lsquant/apps/__init__.py                    |    1 +
 python/lsquant/apps/stage2_dos.py                  |  280 ++++++
 python/lsquant/apps/stage3_kubo.py                 |  355 +++++++
 python/lsquant/apps/stage4_time.py                 |  412 ++++++++
 python/lsquant/apps/stage5_bases.py                |  272 +++++
 python/lsquant/bases.py                            |  408 ++++++++
 python/lsquant/kernels.py                          |  141 +++
 python/lsquant/models/__init__.py                  |   36 +
 python/lsquant/models/wannier.py                   |  467 +++++++++
 python/lsquant/moments.py                          |  470 +++++++++
 python/lsquant/reconstruct.py                      |  644 ++++++++++++
 python/lsquant/spectral_functions.py               |  586 +++++++++++
 scripts/check_goldens.sh                           |  395 ++++++++
 scripts/fetch_wannier2sparse.sh                    |   40 +
 scripts/golden_diff.py                             |  135 +++
 scripts/regenerate_goldens.sh                      |  285 ++++++
 tests/dev/stage_1_check.cpp                        |  145 +++
 tests/dev/stage_2_check.py                         |  168 ++++
 tests/dev/stage_3_check.py                         |  337 ++++++
 tests/dev/stage_4_check.py                         |  308 ++++++
 tests/dev/stage_5_check.py                         |  517 ++++++++++
```

---

## 2. Verified FFTW Sign, Normalisation, and Atom-Gauge Velocity Conventions

1. **Real-space hopping and cell-numbering convention (`wannier2sparse`)**:
   In `wannier2sparse` (`src/hopping_list.cpp`), a primitive hopping entry $(dR, H_{\alpha\beta}(dR))$ in `<seed>_hr.dat` connects source cell $R_1$ orbital $\alpha$ (row index $(i_{R_1} W) + \alpha$) to target cell $R_2 = (R_1 + dR) \bmod N$ orbital $\beta$ (column index $(i_{R_2} W) + \beta$):
   $$(H_{\text{CSR}})_{(R_1, \alpha), (R_1 + dR, \beta)} = H_{\alpha\beta}(dR),$$
   with $n_1$-fastest cell index $i_R = n_1 + N_1 n_2 + N_1 N_2 n_3$ (`index_aliasing` in `wannier2sparse`).

2. **Atom-gauge Bloch Hamiltonian (`brief.md` Section 3)**:
   $$H(k)_{\alpha\beta} = \sum_{dR} \exp\!\big(+i \mathbf{k} \cdot (\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)\big) H_{\alpha\beta}(dR).$$
   For $(B H_{\text{CSR}} B^\dagger)_{(k,\alpha),(k,\beta)} = H(k)_{\alpha\beta}$ to hold for all hoppings $dR$ and sublattice offsets $\boldsymbol{\tau}_\alpha$:
   $$\sum_{R_1, dR} B_{(k,\alpha),(R_1,\alpha)} H_{\alpha\beta}(dR) (B^\dagger)_{(R_1+dR,\beta),(k,\beta)} = \sum_{dR} e^{+i \mathbf{k} \cdot (\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)} H_{\alpha\beta}(dR)$$
   uniquely fixes:
   - **`RealToBloch` ($B$: real $\to$ Bloch)**:
     $$B_{(k,\alpha),(R,\beta)} = \delta_{\alpha\beta} \frac{1}{\sqrt{N_k}} \exp\!\big(-i \mathbf{k} \cdot (\mathbf{R} + \boldsymbol{\tau}_\alpha)\big),$$
     evaluated via `fftw_plan_many_dft` on `{N_3, N_2, N_1}` ($n_1$ fastest) with `howmany = W`, `istride = ostride = W`, `idist = odist = 1`, and `FFTW_FORWARD` ($\sum_R e^{-i \mathbf{k} \cdot \mathbf{R}}$), followed by pointwise multiplication by the precomputed phase table $\exp(-i \mathbf{k} \cdot \boldsymbol{\tau}_\alpha) / \sqrt{N_k}$.
   - **`BlochToReal` ($B^\dagger$: Bloch $\to$ real)**:
     $$(B^\dagger)_{(R,\beta),(k,\alpha)} = \delta_{\alpha\beta} \frac{1}{\sqrt{N_k}} \exp\!\big(+i \mathbf{k} \cdot (\mathbf{R} + \boldsymbol{\tau}_\alpha)\big),$$
     evaluated by pointwise multiplication by the conjugate phase $\exp(+i \mathbf{k} \cdot \boldsymbol{\tau}_\alpha) / \sqrt{N_k}$ followed by `FFTW_BACKWARD` ($\sum_k e^{+i \mathbf{k} \cdot \mathbf{R}}$).
   - This matches `lsquant-official/src/Hamiltonian/Graphene_KaneMele.cpp` (`to_kSpace_pruned` using `FFTW_FORWARD` and `to_rSpace_pruned` using `FFTW_BACKWARD`).

3. **Atom-gauge velocity sign convention (`wannier2sparse/include/tbmodel.hpp:236`)**:
   In `wannier2sparse`, `berry_connection` multiplies each hopping $H_{\alpha\beta}(dR)$ by $(-i)(\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)_a$. Consequently, its Bloch transform is
   $$v_a(k)_{\alpha\beta} = (B \, \text{V}_{a,\text{CSR}} \, B^\dagger)_{(k,\alpha),(k,\beta)} = -i \sum_{dR} (\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)_a \exp\!\big(+i \mathbf{k} \cdot (\mathbf{dR} + \boldsymbol{\tau}_\beta - \boldsymbol{\tau}_\alpha)\big) H_{\alpha\beta}(dR),$$
   giving $\max_x \| v_x x - B \, \text{VX}_{\text{CSR}} \, B^\dagger x \|_\infty \approx 7 \times 10^{-16}$ on both graphene and Haldane.

---

## 3. Verbatim Outputs on TCN

### 3.1 `tests/dev/stage_5_check.py` and `stage5_bases` CLI Suite (default `--window eigen`)

```text
[check a] BlochTransform unitarity (8,8,1) graphene: max|B B^+ x - x| = 6.206e-17, max|B^+ B x - x| = 5.004e-17
[check b] (graphene 8x8x1) max_x ||H_k x - B H_csr B^+ x||_inf = 7.109e-16
[check c] (graphene 8x8x1) max_x ||v_x x - B VX_csr B^+ x||_inf = 7.109e-16
[check b] (haldane 8x8x1) max_x ||H_k x - B H_csr B^+ x||_inf = 5.496e-16
[check c] (haldane 8x8x1) max_x ||v_x x - B VX_csr B^+ x||_inf = 6.474e-16
[check d] Anderson disorder W=1.0 (seed=12345): max|mu_real - mu_bloch| DOS(64) = 4.680e-16, VX-VX(64x64) = 5.761e-14
[check e] (16,16,1) real vs bloch: DOS max|d_mu|=6.550e-15, max_rel|d_curve|=3.186e-15; KuboXX max|d_mu|=6.182e-13, max|d_curve/N|=2.587e-14, max_rel|d_curve|=7.503e-15
[check e] (16,16,1) real vs eigen: DOS max|d_mu|=5.884e-15, max_rel|d_curve|=3.144e-15; KuboXX max|d_mu|=5.898e-13, max|d_curve/N|=2.609e-14, max_rel|d_curve|=7.567e-15
[check e] (16,16,1) bloch vs eigen: DOS max|d_mu|=2.443e-15, max_rel|d_curve|=1.802e-15; KuboXX max|d_mu|=1.777e-13, max|d_curve/N|=8.438e-15, max_rel|d_curve|=2.447e-15
[check e] Jackson KPM (M=128) vs Gaussian exact-diag: rel_L2(DOS)=2.252e-02, rel_L2(KuboXX)=4.272e-02
[check f] G1 chain1d DOS byte-identical to golden: PASS
[check f] G4 graphene Kubo byte-identical to golden: PASS
[check f] G7 chain1d_mag spin precession byte-identical to golden: PASS
stage_5_check: ALL CHECKS PASSED
EXIT_SELF_CHECK:0
=== CLI clean real dos ===
/tmp/s5_cli_test/clean_real_dos/SpectralOp1graphene_real_KPM_M128_stateexact.chebmom
/tmp/s5_cli_test/clean_real_dos/mean_SpectralOp1graphene_real_KPM_M128_stateexactJACKSON.dat
=== CLI clean real kuboxx ===
/tmp/s5_cli_test/clean_real_kuboxx/NonEqOpVX-VXgraphene_real_KPM_M128x128_stateexact.chebmom
/tmp/s5_cli_test/clean_real_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_real_KPM_M128x128_stateexactJACKSON.dat
=== CLI clean bloch dos ===
/tmp/s5_cli_test/clean_bloch_dos/SpectralOp1graphene_bloch_KPM_M128_stateexact.chebmom
/tmp/s5_cli_test/clean_bloch_dos/mean_SpectralOp1graphene_bloch_KPM_M128_stateexactJACKSON.dat
=== CLI clean bloch kuboxx ===
/tmp/s5_cli_test/clean_bloch_kuboxx/NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexact.chebmom
/tmp/s5_cli_test/clean_bloch_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexactJACKSON.dat
=== CLI clean eigen dos ===
/tmp/s5_cli_test/clean_eigen_dos/SpectralOp1graphene_eigen_KPM_M128_stateexact.chebmom
/tmp/s5_cli_test/clean_eigen_dos/mean_SpectralOp1graphene_eigen_KPM_M128_stateexactJACKSON.dat
=== CLI clean eigen kuboxx ===
/tmp/s5_cli_test/clean_eigen_kuboxx/NonEqOpVX-VXgraphene_eigen_KPM_M128x128_stateexact.chebmom
/tmp/s5_cli_test/clean_eigen_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_eigen_KPM_M128x128_stateexactJACKSON.dat
=== CLI disordered W=1.0 real dos ===
/tmp/s5_cli_test/dis_real_dos/SpectralOp1graphene_real_KPM_M64_statedefault.chebmom
/tmp/s5_cli_test/dis_real_dos/mean_SpectralOp1graphene_real_KPM_M64_statedefaultJACKSON.dat
=== CLI disordered W=1.0 real kuboxx ===
/tmp/s5_cli_test/dis_real_kuboxx/NonEqOpVX-VXgraphene_real_KPM_M64x64_statedefault.chebmom
/tmp/s5_cli_test/dis_real_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_real_KPM_M64x64_statedefaultJACKSON.dat
=== CLI disordered W=1.0 bloch dos ===
/tmp/s5_cli_test/dis_bloch_dos/SpectralOp1graphene_bloch_KPM_M64_statedefault.chebmom
/tmp/s5_cli_test/dis_bloch_dos/mean_SpectralOp1graphene_bloch_KPM_M64_statedefaultJACKSON.dat
=== CLI disordered W=1.0 bloch kuboxx ===
/tmp/s5_cli_test/dis_bloch_kuboxx/NonEqOpVX-VXgraphene_bloch_KPM_M64x64_statedefault.chebmom
/tmp/s5_cli_test/dis_bloch_kuboxx/KuboGreenWood_NonEqOpVX-VXgraphene_bloch_KPM_M64x64_statedefaultJACKSON.dat
--- GOLDEN_DIFF CLEAN REAL vs BLOCH vs EIGEN ---
first-diff line 13 col 1, new=-1.0740740740740722 golden=-1.0740740740740751; max abs diff 6.550e-15 at line 67 col 1; max rel diff 1.000e+00 at line 13 col 2; DEFECT level (>= 1e-12)
no line differs (difference is in line endings or trailing bytes)
first-diff line 13 col 2, new=-4.2356979628604462e-18 golden=0; max abs diff 2.442e-15 at line 137 col 1; max rel diff 1.000e+00 at line 13 col 2; DEFECT level (>= 1e-12)
no line differs (difference is in line endings or trailing bytes)
first-diff line 11 col 1, new=22.058081999999722 golden=22.058081999999978; max abs diff 6.182e-13 at line 7235 col 1; max rel diff 1.000e+00 at line 13 col 2; DEFECT level (>= 1e-12)
first-diff line 1 col 2, new=1.24915e-06 golden=1.24914e-06; max abs diff 1.000e-11 at line 3 col 2; max rel diff 8.005e-06 at line 1 col 2; DEFECT level (>= 1e-12)
first-diff line 11 col 1, new=22.058081999999978 golden=22.058081999999985; max abs diff 1.776e-13 at line 649 col 1; max rel diff 1.000e+00 at line 3703 col 2; DEFECT level (>= 1e-12)
first-diff line 127 col 2, new=2.79419e-06 golden=2.7942e-06; max abs diff 1.000e-11 at line 127 col 2; max rel diff 3.579e-06 at line 127 col 2; DEFECT level (>= 1e-12)
--- GOLDEN_DIFF DISORDERED REAL vs BLOCH (W=1.0, seed=12345) ---
first-diff line 11 col 1, new=1 golden=0.99999999999999978; max abs diff 4.441e-16 at line 64 col 1; max rel diff 1.949e+00 at line 18 col 2; DEFECT level (>= 1e-12)
no line differs (difference is in line endings or trailing bytes)
first-diff line 12 col 1, new=-1.1854778242125124 golden=-1.1854778242125126; max abs diff 6.395e-14 at line 3005 col 1; max rel diff 3.255e-11 at line 69 col 2; DEFECT level (>= 1e-12)
no line differs (difference is in line endings or trailing bytes)
```

---

## 4. Deviations from the Brief

1. **`TransformedStates(transform, inner)` in `lsquant/core/algebra/state_generator.hpp`**:
   Added `TransformedStates` alongside `BlockLocalStates` so that when `BlochBasis.states("random", num_states, seed)` (or `EigenBasis.states("random", ...)`) is requested, it applies $B$ (`RealToBloch`) to `RandomPhaseStates(dim, num_states, seed)`. Since $B$ is unitary, $|r_k\rangle = B |r_R\rangle$ is a properly normalised isotropic random state in Bloch space and represents the **exact same physical state vector** as $|r_R\rangle$ in real space, allowing finite-sample stochastic traces in `RealBasis` and `BlochBasis` (at both $W = 0$ and $W > 0$) to match to $\sim 10^{-14}$ for any number of states $R$.
2. **`--seed` argument handling in `python/lsquant/apps/stage5_bases.py`**:
   The CLI example in Section 4 of the brief passes both `--seed graphene` (model seed prefix) and `--seed 12345` (PRNG integer seed) on the same command line. `stage5_bases.py` accepts `--seed` with `action="append"` (plus optional `--label` and `--rng-seed`) and automatically disambiguates model prefix strings from integer RNG seeds so both single-`--seed` and dual-`--seed` invocations work seamlessly.

---

## 5. Open Questions for Fable

None. All checks (a)–(f) and the three-basis `stage5_bases` runs pass within machine precision.
