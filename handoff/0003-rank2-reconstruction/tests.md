# Handoff 0003 — tester report (stage 3 gate: G4 and G5 from Python; tutorial 03)

Tester: Codex role (Opus subagent driven by Fable). Machine: TCN (`fatnode-tacn`),
checkout `~/lsquant-unified`, branch `handoff/0003-codex` cut from `cpp-unified`
(`e1e6b28`). Nothing under `lsquant/`, `python/lsquant/`, `tests/dev/`, the legacy tree,
`goldens/`, a seed or a tolerance was edited. Every run used an explicit
`OMP_NUM_THREADS` (16, and 1 where the gate says so; the stage-1 part of the harness
keeps its own 1, 8 and 64).

The developer branch was tested twice:

| Run | `origin/handoff/0003-agi` | merge into `handoff/0003-codex` | logs |
|---|---|---|---|
| 1 | `c01b026` | fast-forward, no conflict | `figures/out/stage_3/*_run1_c01b026.log` |
| 2 (final) | `0774db7` (review-1 change: keyword-only `kubo_plan(M, *, left_index, right_index, word)`, `dos_plan(M, N, word)`) | fast-forward, no conflict (`merge_run2.log`) | `figures/out/stage_3/{build,stage_3_check,check_goldens_stage3,figure_run}.log` |

**Result (final run, `0774db7`): GATE PASS (stage 3), 41 of 41 files IDENTICAL at
OMP_NUM_THREADS 16 and 1; both greps empty. Run 1 (`c01b026`) gave the same 41
IDENTICAL verdicts. Tutorial 03 from Python: 1 of 4 files IDENTICAL, 3 DIFFER
(Bastin, Sea, Greenwood); C = +1.017016743089, equal to the golden's C to all printed
digits. The cause is localised (section 6): `tutorial03.py` reconstructs from the
in-memory moments, while the legacy tutorial reconstructed from the `.chebmom2D` file
it had written at legacy text precision. Through that write/read round trip all four
files are IDENTICAL.**

## 1. Build (fresh, both runs)

```
$ rm -rf build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-unified/build
cmake exit: 0
$ cmake --build build -j 32
[100%] Built target stage_1_check
[100%] Built target _core
build exit: 0
```

Full logs: `figures/out/stage_3/build.log` (final) and `build_run1_c01b026.log`.
`compat/legacy_kernels.cpp` compiles into `_core` with the same `-O3 -Wall -Wextra`
flags (`-flto=auto` on `_core`). Warnings: **0 under `lsquant/`** (core, io, api,
compat, apps). The 101 warnings in each log are all in the legacy tree, as in stages 1
and 2: 62 in `include/chebyshev_moments.hpp`, 32 in `include/units.hpp`, 7 in
`src/chebyshev_momentsLocal.cpp`.

## 2. Developer self-check (final run, `0774db7`)

```
$ PYTHONPATH=build/lsquant_core:python OMP_NUM_THREADS=16 python3 tests/dev/stage_3_check.py
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
wall 27.12 s
exit: 0
```

Exit 0. Run 1 (`c01b026`) printed the same lines (only the wall time differs, 27.42 s).
It matches the developer's report.

## 3. Item 0 of the brief (`SparseMatrix.from_csr`), tester check

Done with a scratch script that is not committed (the brief asks for this check by name;
it is not a test suite). It parses `goldens/stage0/chain1d/operators/chain1d.HAM.CSR`
into a `scipy.sparse.csr_matrix` (scipy gives `int32` indices), passes the arrays
converted to `int64`/`complex128` to `from_csr`, and then tries the wrong buffers.
Output at `0774db7`:

```
scipy csr: shape (2048, 2048) nnz 4096 dtypes int32 int32 complex128
int64/complex128 from_csr + hermitian_fold: dos moments array_equal vs ReadCsrText: True
indptr int32                              : TypeError: SparseMatrix.from_csr: argument 'row_pointers' must be a 1D C-contiguous int64 numpy.ndarray
indices int32                             : TypeError: SparseMatrix.from_csr: argument 'columns' must be a 1D C-contiguous int64 numpy.ndarray
scipy native int32 indptr/indices         : TypeError: SparseMatrix.from_csr: argument 'row_pointers' must be a 1D C-contiguous int64 numpy.ndarray
indices non-contiguous int64 (stride 2)   : TypeError: SparseMatrix.from_csr: argument 'columns' must be a 1D C-contiguous int64 numpy.ndarray
values complex64                          : TypeError: SparseMatrix.from_csr: argument 'values' must be a 1D C-contiguous complex128 numpy.ndarray
values float64                            : TypeError: SparseMatrix.from_csr: argument 'values' must be a 1D C-contiguous complex128 numpy.ndarray
indptr python list                        : TypeError: SparseMatrix.from_csr: argument 'row_pointers' must be a 1D C-contiguous int64 numpy.ndarray
```

Source check: the binding validates each buffer (`ValidateCsrBuffer1D`, forcecast
off) and passes the raw pointers to the new pointer overload
`SparseMatrix::FromCsr(dim, rpSize, rpPtr, colSize, colPtr, valSize, valPtr)`. That
overload assigns an `Eigen::Map<const SparseMatrix>` into the owned matrix, which is
the single accepted copy. There is no `std::vector` staging on the Python path. The
by-value `std::vector` overload remains for C++ callers and delegates to the pointer
overload.

## 4. Gate (`scripts/check_goldens.sh 3 build`)

Harness extension (this handoff, `scripts/check_goldens.sh`): stage 3 runs stage 1 and
stage 2 in full (their comparisons and greps). It then runs
`python3 -m lsquant.apps.stage3_kubo` at `OMP_NUM_THREADS=16` and again at `1`, with
`PYTHONPATH=build/lsquant_core:python`:

1. G5 moments: Haldane `--left VX --right VY` and `--right VX`, `--window -4.5 4.5
   --trace exact --states 512 --moments 128`; `cmp` both `.chebmom2D` against
   `goldens/stage0/haldane/`.
2. G5 curves from the golden moment files at 10 meV: `bastin`, `sea` and `surf` on
   `VX-VY`, `greenwood` on `VX-VX`.
3. G4 moments: graphene `VX VX`, `--window gershgorin --trace random --seed 12345
   --states 1 --moments 64`.
4. G4 curves from the golden moment file: `bastin` at 64 meV, `sea` and `surf` at
   10 meV (`sea_surf/`), `greenwood_fft` and `bastin_fft` at 10 meV (`fft_grid/`).
5. `grep -rn "compat/" lsquant/core lsquant/io` must be empty. The stage-2 code-line
   grep for `delta|DOS|Kubo|basis|std::cout|getenv` on `lsquant/core` and `lsquant/io`
   is part of the stage-2 rerun.

Command: `OMP_NUM_THREADS=16 bash scripts/check_goldens.sh 3 build 2>&1 | tee
figures/out/stage_3/check_goldens_stage3.log`.

### 4.1 Final run (`0774db7`), verbatim — exit 0

```
== stage 1: /home/ICN2/jgarcia/lsquant-unified/build/lsquant_core/stage1_moments vs goldens/stage0/chain1d (KPM_SEED=12345)
IDENTICAL [T=1] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=1] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=8] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=8] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=64] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=64] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=8, batch 16] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2, batched)
IDENTICAL [random, run2 vs run1] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom
IDENTICAL [random, run2 vs run1] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D
IDENTICAL [random, run2 vs run1] SpectralOp1chain1dKPM_M64_statedefault.chebmom
IDENTICAL [random, run2 vs run1] SpectralOp1chain1dKPM_M64_statedefault.chebmom1D
GREP      grep -rn "basis\|std::cout\|getenv" lsquant/core: empty
== stage 2: python3 -m lsquant.apps.stage2_dos, PYTHONPATH=/home/ICN2/jgarcia/lsquant-unified/build/lsquant_core:python, vs goldens/stage0/chain1d (KPM_SEED=12345)
          [T=1] (a) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=1] (a) --moments G1 -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=1] (b) printed: SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
          [T=1] (b) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=1] (b) --operators -> chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=1] (b) --operators -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=16] (a) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=16] (a) --moments G1 -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=16] (b) printed: SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
          [T=16] (b) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=16] (b) --operators -> chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=16] (b) --operators -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
GREP      grep -rn "delta\|DOS\|Kubo\|basis\|std::cout\|getenv" lsquant/core lsquant/io, code lines: empty
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:27: * recording whether the `(2 - delta_{m0})` KPM weights (`kpmWeighted`) and Hermitian
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:32: *      Never apply KPM `(2 - delta_{m0})` weights or symmetrisation twice; check
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:158: * @brief Apply the legacy rank-2 KPM `(2 - delta_{m0})` weights and optional Hermitian symmetrisation.
== stage 3: python3 -m lsquant.apps.stage3_kubo, PYTHONPATH=/home/ICN2/jgarcia/lsquant-unified/build/lsquant_core:python, vs goldens/stage0/{haldane,graphene}
          [T=16] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=16] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
          [T=16] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=16] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL [T=16] haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
IDENTICAL [T=16] haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
          [T=16] G5 bastin printed: KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 sea printed: KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 surf printed: KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 greenwood printed: KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL [T=16] haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
          [T=16] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
          [T=16] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom
IDENTICAL [T=16] graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D (G4 moments)
          [T=16] G4 bastin 64 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 sea 10 meV printed: KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 surf 10 meV printed: KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 greenwood_fft 10 meV printed: KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
          [T=16] G4 bastin_fft 10 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
IDENTICAL [T=16] graphene/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 64 meV)
IDENTICAL [T=16] graphene/sea_surf/KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/sea_surf/KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/fft_grid/KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/fft_grid/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
          [T=1] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=1] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
          [T=1] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=1] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL [T=1] haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
IDENTICAL [T=1] haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
          [T=1] G5 bastin printed: KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 sea printed: KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 surf printed: KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 greenwood printed: KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL [T=1] haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
          [T=1] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
          [T=1] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom
IDENTICAL [T=1] graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D (G4 moments)
          [T=1] G4 bastin 64 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 sea 10 meV printed: KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 surf 10 meV printed: KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 greenwood_fft 10 meV printed: KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
          [T=1] G4 bastin_fft 10 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
IDENTICAL [T=1] graphene/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 64 meV)
IDENTICAL [T=1] graphene/sea_surf/KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/sea_surf/KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/fft_grid/KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/fft_grid/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
GREP      grep -rn "compat/" lsquant/core lsquant/io: empty
== 41 files compared; failures: 0
GATE PASS (stage 3)
```

### 4.2 Run 1 (`c01b026`, before the review-1 change), verbatim — exit 0

```
== stage 1: /home/ICN2/jgarcia/lsquant-unified/build/lsquant_core/stage1_moments vs goldens/stage0/chain1d (KPM_SEED=12345)
IDENTICAL [T=1] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=1] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=8] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=8] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=64] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=64] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=8, batch 16] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2, batched)
IDENTICAL [random, run2 vs run1] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom
IDENTICAL [random, run2 vs run1] NonEqOpVX-VXchain1dKPM_M64x64_statedefault.chebmom2D
IDENTICAL [random, run2 vs run1] SpectralOp1chain1dKPM_M64_statedefault.chebmom
IDENTICAL [random, run2 vs run1] SpectralOp1chain1dKPM_M64_statedefault.chebmom1D
GREP      grep -rn "basis\|std::cout\|getenv" lsquant/core: empty
== stage 2: python3 -m lsquant.apps.stage2_dos, PYTHONPATH=/home/ICN2/jgarcia/lsquant-unified/build/lsquant_core:python, vs goldens/stage0/chain1d (KPM_SEED=12345)
          [T=1] (a) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=1] (a) --moments G1 -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=1] (b) printed: SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
          [T=1] (b) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=1] (b) --operators -> chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=1] (b) --operators -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=16] (a) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=16] (a) --moments G1 -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
          [T=16] (b) printed: SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
          [T=16] (b) printed: meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat
IDENTICAL [T=16] (b) --operators -> chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=16] (b) --operators -> chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat (G3)
GREP      grep -rn "delta\|DOS\|Kubo\|basis\|std::cout\|getenv" lsquant/core lsquant/io, code lines: empty
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:27: * recording whether the `(2 - delta_{m0})` KPM weights (`kpmWeighted`) and Hermitian
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:32: *      Never apply KPM `(2 - delta_{m0})` weights or symmetrisation twice; check
          excluded comment line: lsquant/core/chebyshev/chebyshev_moments.hpp:158: * @brief Apply the legacy rank-2 KPM `(2 - delta_{m0})` weights and optional Hermitian symmetrisation.
== stage 3: python3 -m lsquant.apps.stage3_kubo, PYTHONPATH=/home/ICN2/jgarcia/lsquant-unified/build/lsquant_core:python, vs goldens/stage0/{haldane,graphene}
          [T=16] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=16] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
          [T=16] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=16] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL [T=16] haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
IDENTICAL [T=16] haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
          [T=16] G5 bastin printed: KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 sea printed: KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 surf printed: KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=16] G5 greenwood printed: KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL [T=16] haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=16] haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
          [T=16] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
          [T=16] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom
IDENTICAL [T=16] graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D (G4 moments)
          [T=16] G4 bastin 64 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 sea 10 meV printed: KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 surf 10 meV printed: KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=16] G4 greenwood_fft 10 meV printed: KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
          [T=16] G4 bastin_fft 10 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
IDENTICAL [T=16] graphene/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 64 meV)
IDENTICAL [T=16] graphene/sea_surf/KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/sea_surf/KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/fft_grid/KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
IDENTICAL [T=16] graphene/fft_grid/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
          [T=1] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=1] G5 VX-VY moments printed: NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom
          [T=1] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
          [T=1] G5 VX-VX moments printed: NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom
IDENTICAL [T=1] haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
IDENTICAL [T=1] haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D (G5 moments)
          [T=1] G5 bastin printed: KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 sea printed: KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 surf printed: KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
          [T=1] G5 greenwood printed: KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
IDENTICAL [T=1] haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
IDENTICAL [T=1] haldane/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat (G5 curve)
          [T=1] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D
          [T=1] G4 VX-VX moments printed: NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom
IDENTICAL [T=1] graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D (G4 moments)
          [T=1] G4 bastin 64 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 sea 10 meV printed: KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 surf 10 meV printed: KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat
          [T=1] G4 greenwood_fft 10 meV printed: KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
          [T=1] G4 bastin_fft 10 meV printed: KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat
IDENTICAL [T=1] graphene/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 64 meV)
IDENTICAL [T=1] graphene/sea_surf/KuboBastinSea_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/sea_surf/KuboBastinSurf_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/fft_grid/KuboGreenWood_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
IDENTICAL [T=1] graphene/fft_grid/KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON_FFTgrid.dat (G4 curve, 10 meV)
GREP      grep -rn "compat/" lsquant/core lsquant/io: empty
== 41 files compared; failures: 0
GATE PASS (stage 3)
```

The two runs give identical verdicts, line by line.

## 5. Figure (`figures/stage_3.py`)

`MPLBACKEND=Agg python3 figures/stage_3.py --build build --threads 16` runs
`python3 -m lsquant.apps.tutorial03 --out <tmp>`. It compares the four `.dat` files
by name with `goldens/stage0/tutorials/03/`, which are byte-identical to the G5
`.dat` files in `goldens/stage0/haldane/` (md5 checked). It then draws the two panels.
Output: `figures/out/stage_3/fig_stage_3.png` and `fig_stage_3.caption.txt`; run log
in `figure_run.log`. The script exits 1 because files differ, and says so in the
caption.

tutorial03 output and `cmp` list (final run):

```
tutorial03 printed (5 lines):
  <tmp>/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  <tmp>/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  <tmp>/KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
  <tmp>/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
  Chern from gap plateau: C = +1.017 (+1.017016743089)
cmp list:
  KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS
  KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS
  KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
  KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS
C new = +1.017016743089, C golden = +1.017016743089
wrote /home/ICN2/jgarcia/lsquant-unified/figures/out/stage_3/fig_stage_3.png
exit: 1
```

**C = +1.017016743089** from the new Bastin curve (tutorial03 prints `C = +1.017`),
and **C = +1.017016743089** from the golden curve.

Caption (verbatim, `fig_stage_3.caption.txt`):

> FIG. 3. Stage 3 (handoff 0003): tutorial 03, the Hall and longitudinal conductivity of the Haldane model (t1 = -1, t2 = 0.15, phi = pi/2), reproduced from Python by python3 -m lsquant.apps.tutorial03 (lsquant package over the pybind11 module lsquant._core, OMP_NUM_THREADS = 16): N = 512 orbitals, window [-4.5, 4.5] eV, exact trace over all N site states, M = 128 x 128 Chebyshev moments (VX-VY for sigma_xy, VX-VX for sigma_xx), Jackson kernel at eta = 10 meV; conductivities in units of e^2/h as lsqhall.py normalises them (divided by (N/2) A_cell e^2/h, A_cell = 3 sqrt(3)/2, e = hbar = 1). (a) sigma_xy: Kubo-Bastin total, Fermi-sea and Fermi-surface parts; new as lines, the stage-0 golden (goldens/stage0/tutorials/03) as open markers (48 evenly spaced points per curve), dotted line at 1; residual new - golden below. (b) sigma_xx (Kubo-Greenwood), same encoding. Chern number from the gap plateau, C = (mean over |E| < 0.5 eV of sigma_xy / A) 2 pi with A = (N/2) A_cell: C = +1.017016743089 from the new curve (tutorial03 prints +1.017016743089), C = +1.017016743089 from the golden curve. Max |new - golden| over the four curves: 9.4e-13 e^2/h, max energy-axis difference 0.0e+00 eV. Byte comparison with goldens/stage0/tutorials/03 by name (cmp criterion), 1 of 4 IDENTICAL: KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS; KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS; KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL; KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS. AT LEAST ONE FILE DIFFERS.

## 6. Localisation of the tutorial-03 DIFFERS (not fixed)

`cmp` then `golden_diff.py` on the four files (final run). Surf: `cmp` prints nothing because the file is IDENTICAL; the "no line differs" line is `golden_diff.py` run on an identical pair.

```
== KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
<out>/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat goldens/stage0/tutorials/03/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat differ: byte 14, line 1
first-diff line 1 col 2, new=2.18661e-08 golden=2.1866e-08; max abs diff 1.000e-10 at line 186 col 2; max rel diff 9.711e-06 at line 186 col 2; DEFECT level (>= 1e-12)
== KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
<out>/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat goldens/stage0/tutorials/03/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat differ: byte 35, line 2
first-diff line 2 col 2, new=4.60745e-08 golden=4.60743e-08; max abs diff 1.000e-10 at line 411 col 2; max rel diff 3.567e-01 at line 3839 col 2; DEFECT level (>= 1e-12)
== KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat
no line differs (difference is in line endings or trailing bytes)
== KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat
<out>/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat goldens/stage0/tutorials/03/KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat differ: byte 13, line 1
first-diff line 1 col 2, new=5.06063e-08 golden=5.06082e-08; max abs diff 1.900e-12 at line 1 col 2; max rel diff 3.754e-05 at line 1 col 2; DEFECT level (>= 1e-12)
```

(Abs diffs are in raw units: 1e-10 raw is 9.4e-13 e^2/h. The large relative differences are on values near zero printed with 6 significant
digits. Every residual is at most 9.4e-13 e^2/h, and the energy axis is
identical.)

Diagnosis. The operators that `make_haldane` stages from `test/golden/haldane` are
byte-identical to `goldens/stage0/haldane/operators`. The `.desc` window is exactly
[-4.5, 4.5]. The gate shows the moments (`.chebmom2D`) and each formula applied to
the golden moment file are IDENTICAL. So moments and formulas are right. The
in-memory moments that `run_plan` returns differ from the golden file at the last
bits:

- `VX-VY`: max |dIm| 4.9e-15, real parts at 1e-31.
- `VX-VX`: max |dRe| 5.3e-15.
- About 15,400 of 16,384 elements are not bitwise equal.

The golden file holds them at the legacy text precision. The legacy `lsqhall.py`
wrote the `.chebmom2D` and reconstructed by running the `*FromChebmom` drivers on that
file. `tutorial03.py` instead passes the unrounded `ChebyshevMoments` straight to
`kubo_bastin` and the other formulas. Surf is IDENTICAL because it uses only
`real(mu)` of `VX-VY`, which sits at 1e-31.

Direct check (scratch script, same build):

```
in-memory moments (as tutorial03 does):
  KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS
  KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS
  KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
  KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat: DIFFERS
moments after WriteLegacyChebmom2D + ReadChebmom:
  KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
  KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
  KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
  KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat: IDENTICAL
```

Where: `python/lsquant/apps/tutorial03.py`, `run_tutorial03`. `mu_xy` and `mu_xx` go
to the reconstruction functions without the legacy file round trip. The gate
program `stage3_kubo` form 2 + `--formula` has the same property: it reconstructs from
the in-memory `mu` after writing the file. The gate does not exercise that path,
because it reconstructs from the golden moment files. Self-check (a) compares
moments only after `_roundtrip_legacy_2d`, so it cannot show this. The developer's
report and review 1 ("tutorial 03 gives C = +1.017, the legacy value") are true for C,
but the four files are not byte-identical. The decision on whether tutorial 03
should reconstruct from the written file (the legacy pipeline) belongs to Fable.

## 7. Smells

- **Formula or grid in `core/`: none found.** `contraction.cpp` gains only the three
  rank-2 kernels (`m0` outer, `m1` inner, OpenMP over `i`, tables passed in). A grep of
  `lsquant/core` and `lsquant/io` for `pow|acos|0.95|M_PI|cos(|sqrt(` on code lines
  finds only the stage-1 random-phase generator (`state_generator.cpp`) and the
  Lanczos norms (`spectral_bounds.cpp`). The stage-2 code-line grep is empty (three
  comment lines excluded, listed in the gate output).
- **`compat` included from `core/`: no.** `grep -rn "compat/" lsquant/core lsquant/io`
  is empty. `compat/legacy_kernels.cpp` is compiled only into `_core`.
- **numpy where `math` was required: none found.** The prefactors, `ceil(pi/eta)`,
  the uniform and FFT-node grids, the 2D Jackson weights and the cumulative integrals
  are scalar `math` loops in `reconstruct.py` and `kernels.py`. numpy is used only for
  containers and for the elementwise table algebra of Sea (`np.conj(Gr)`, `Gr - Ga`,
  `DGr + DGa`). Those operations are exact IEEE operations, the same ones the legacy
  performs, and the G4/G5 Sea curves are byte-identical.
- **Table built with Python's complex power: no.** `greenR`/`DgreenR` come from
  `_core.compat.legacy_green_tables` (libstdc++ `std::pow(complex<double>, double)`).
  The sympy `GreenRetarded`/`GreenRetardedDerivative` are only compared (1.1e-14 abs;
  7.2e-15 rel).
- **Copy in `from_csr`:** only the accepted Eigen `Map` copy remains; strict dtypes and
  contiguity are enforced with a `TypeError` naming the argument (section 3).
- **tutorial03 not byte-identical** (section 6): the reconstruction skips the legacy
  write/read of the moments. The self-check's round trip hides this.
- `tutorial03.py` stages from `test/golden/haldane/operators` in the legacy tree,
  through `make_haldane.GOLDEN`, rather than from `goldens/stage0/haldane/operators`.
  The two are byte-identical today, but stage 7 deletes the legacy tree.
- No golden, seed or tolerance was touched; `goldens/` is unchanged.
