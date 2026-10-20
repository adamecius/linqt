# Handoff 0004 Report — Stage 4: Apply Mode, Time Evolution, MSD, and Projected Words (G6, G7, G8)

## 1. Files Added and Modified

- `lsquant/core/chebyshev/plan.hpp` (modified: added `Register = {Ket, Bra, Aux}`, `ApplyDerivedSeries`, `EvolveRegister`, `EvolveDerived`, `CopyRegister`, `TimeEvolutionStep`, `TimeAxis`, and optional `scale` field on `ApplySeries` and `ApplyDerivedSeries` for the post-recursion band-center phase `scal`)
- `lsquant/core/chebyshev/chebyshev_moments.hpp` (modified: added `timeStep` and `numberOfTimeSteps` metadata fields to `ChebyshevMoments`)
- `lsquant/core/chebyshev/moment_engine.hpp`, `lsquant/core/chebyshev/moment_engine.cpp` (modified: implemented `ApplySeries`, `ApplyDerivedSeries`, generalised `MirrorToBra`, and `TimeAxis` with `{Ket, Bra, Aux}` registers in `MomentEngine::Run`)
- `lsquant/io/chebmom_file.hpp`, `lsquant/io/chebmom_file.cpp` (modified: added `WriteLegacyChebmomTD` and `.chebmomTD` auto-detection in `ReadChebmom`)
- `lsquant/compat/legacy_kernels.hpp`, `lsquant/compat/legacy_kernels.cpp` (modified: added `ComputeLegacyBesselCoefficients` using `gsl_sf_bessel_Jn`)
- `lsquant/CMakeLists.txt` (modified: linked `GSL::gsl` to `_core`; `lsquant_core` remains free of GSL)
- `lsquant/api/bindings.cpp` (modified: bound `Register`, `ApplySeries`, `ApplyDerivedSeries`, `EvolveRegister`, `EvolveDerived`, `CopyRegister`, `TimeAxis`, `ChebyshevMoments.time_step`, `ChebyshevMoments.number_of_time_steps`, `WriteLegacyChebmomTD`, and `_core.compat.bessel_series_coefficients`)
- `python/lsquant/_core.pyi` (modified: added type stubs for stage-4 bindings)
- `python/lsquant/spectral_functions.py` (modified: added `TimeEvolutionFunction`, `TimeEvolution`, `bessel_series_coefficients`, and `legacy_bessel_coefficients`)
- `python/lsquant/moments.py` (modified: added `HBAR`, `time_step_from_tmax`, `time_evolution_series`, `msd_plan`, `projected_evolution_plan`, `time_correlation_plan`, and updated `run_plan` for `TimeAxis` tensors)
- `python/lsquant/reconstruct.py` (modified: added `time_correlation_curve`)
- `python/lsquant/__init__.py` (modified: exported stage-4 symbols)
- `python/lsquant/apps/stage4_time.py` (added: stage-4 gate CLI)
- `python/lsquant/apps/tutorial02.py` (added: Tutorial 02 1D Anderson localisation application)
- `python/lsquant/apps/tutorial04.py` (added: Tutorial 04 Larmor spin precession application)
- `python/lsquant/apps/tutorial05.py` (added: Tutorial 05 2D square-lattice charge transport application)
- `tests/dev/stage_4_check.py` (added: stage-4 developer self-check)

---

## 2. CMake Configure and Build on TCN (`~/lsquant-agi-build`)

```
$ cd ~/lsquant-agi-build
$ cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -DLSQUANT_BUILD_INPUT_TOOLS=OFF
$ cmake --build build --target lsquant_core stage1_moments stage_1_check _core -j 16
```

Verbatim configure and build output:

```
-- Found pybind11: /home/ICN2/jgarcia/.local/lib/python3.10/site-packages/pybind11/include (found version "3.1.0")
-- Configuring done
-- Generating done
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build/build
Consolidate compiler generated dependencies of target lsquant_core
[  0%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[  8%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[ 16%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[ 25%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[ 33%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 41%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 50%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 66%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/core/chebyshev/contraction.cpp.o
[ 66%] Building CXX object lsquant_core/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 66%] Linking CXX static library ../lib/liblsquant_core.a
[100%] Built target lsquant_core
Consolidate compiler generated dependencies of target lsquant_core
[ 92%] Built target lsquant_core
Consolidate compiler generated dependencies of target stage1_moments
[ 92%] Building CXX object lsquant_core/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[100%] Linking CXX executable stage1_moments
[100%] Built target stage1_moments
[ 85%] Built target lsquant_core
Consolidate compiler generated dependencies of target stage_1_check
[ 92%] Building CXX object lsquant_core/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[100%] Linking CXX executable stage_1_check
[100%] Built target stage_1_check
[ 85%] Built target lsquant_core
Consolidate compiler generated dependencies of target _core
[ 85%] Building CXX object lsquant_core/CMakeFiles/_core.dir/api/bindings.cpp.o
[ 92%] Building CXX object lsquant_core/CMakeFiles/_core.dir/compat/legacy_kernels.cpp.o
[100%] Linking CXX shared module _core.cpython-310-x86_64-linux-gnu.so
[100%] Built target _core
```

Zero compiler warnings were emitted under `lsquant/` with `-Wall -Wextra`.

---

## 3. Bounds Source of G7 (`chain1d_mag`) and Legacy `RescaleCon2ChebyshevDomain` Defect Note

### Bounds Source of G7 (`chain1d_mag`)

In `src/inline_compute-kpm-TimeEvProjetedOp.cpp` (line 64), the legacy driver calls:

```cpp
spectral_bounds = chebyshev::utility::SpectralBounds(OP[0]);
```

without passing a `descriptor_path` argument (which defaults to `""`), and `test/spin_precession.sh` stages only `operators/chain1d_mag.{HAM,SX,SY,SZ,PX}.CSR` without a `BOUNDS` file. Consequently, `chebyshev::utility::SpectralBounds` takes branch 3 (the Gershgorin row-sum enclosure with `PAD = 0.10`): for `chain1d_mag.HAM.CSR` (`t = -1.0`, Zeeman half-splitting `0.1`), the maximum absolute row sum is `rho = 2 * 1.0 + 0.1 = 2.1`, so `[-rho * (1 + 0.10), +rho * (1 + 0.10)] = [-2.31, 2.31]` and `BandWidth = 4.62`, `BandCenter = 0`, matching the G7 golden header `128 4.62 0 96 1 ` (`window_from_gershgorin(ham)` reproduces this window directly).

### Legacy `RescaleCon2ChebyshevDomain` Defect (for `AGENTS.md` section 10)

In `include/chebyshev_moments.hpp` (lines 422-428), `MomentsTD::RescaleCon2ChebyshevDomain()` rescales the commutator / velocity matrix `Conmuter()` via:

```cpp
this->Conmuter().Rescale(this->ScaleFactor(), this->ShiftFactor());
```

which applies `V <- ScaleFactor * V + ShiftFactor * I` (`RescaledOperator(V, window)`), adding `ShiftFactor * I = -(BandCenter / HalfWidth) * I` to the diagonal of `[X, H~]`. Because `[X, a H + b I] = a [X, H]`, only the multiplicative `ScaleFactor` should rescale the commutator; adding `ShiftFactor * I` is only harmless when `BandCenter == 0` (`ShiftFactor == 0`, as in symmetric windows `[-W, W]`). For byte-for-byte compatibility with legacy MSD goldens, `stage4_time.py` and the tutorial apps wrap the velocity CSR matrix in `_core.RescaledOperator(vel_mat, window)`.

---

## 4. Developer Self-Check Output (`stage_1_check` and `tests/dev/stage_4_check.py`)

```
$ ./build/lsquant_core/stage_1_check
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES

$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -u tests/dev/stage_4_check.py
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
```

### Exact Legacy Command Used in Self-Check (d)

In a temporary working directory staging `goldens/stage0/chain1d_ballistic/operators/chain1d.{HAM,VX}.CSR`, `BOUNDS` (`-2 2`), and an `exact` local-state file over all `256` sites (`local\n256\n0\n...\n255\n`):

```
KPM_SEED=12345 OMP_NUM_THREADS=1 /home/ICN2/jgarcia/lsquant-unified/build/inline_compute-kpm-CorrOp chain1d VX VX 64 16 5 exact
```

which writes `TimeCorrVX-VXchain1dKPM_M64_stateexact.chebmomTD`.

### Review 1 Self-Check Rerun (`tests/dev/stage_4_check.py`)

After updating `time_evolution_series` to default to `use_legacy_bessel=False` (with `stage4_time.py`, `tutorial02.py`, `tutorial04.py`, `tutorial05.py`, and the legacy-golden checks in `stage_4_check.py` passing `use_legacy_bessel=True` explicitly) and documenting the empty-`coeffs` zero-map semantics on `ApplySeries`, `ApplyDerivedSeries`, and `msd_plan`:

```
$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -u tests/dev/stage_4_check.py
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
```

---

## 5. Gate CLI (`stage4_time.py`) and Tutorial 04 Smoke Run on TCN

```
$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage4_time \
    --operators goldens/stage0/chain1d_ballistic/operators --label chain1d --velocity VX \
    --window -2 2 --trace exact --states 256 --moments 256 --steps 64 --tmax 20 --word msd \
    --out /tmp/agi_s4_cli
/tmp/agi_s4_cli/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD
/tmp/agi_s4_cli/CorrelationVXchain1dKPM_M256_state20exact.chebmom

$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.stage4_time \
    --moments /tmp/agi_s4_cli/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD \
    --broadening-meV 10 --fermi-energy 0.0 --out /tmp/agi_s4_cli
/tmp/agi_s4_cli/meanCorrelationVXchain1dKPM_M256_state20exactEF0.000000JACKSON.dat

$ OMP_NUM_THREADS=1 PYTHONPATH=build/lsquant_core:python python3 -m lsquant.apps.tutorial04 \
    --out /tmp/agi_s4_cli/t04
/tmp/agi_s4_cli/t04/TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD
/tmp/agi_s4_cli/t04/TimeProjSY-PXchain1d_magKPM_M32_stateexact.chebmomTD
/tmp/agi_s4_cli/t04/TimeProjSZ-PXchain1d_magKPM_M32_stateexact.chebmomTD
```

---

## 6. Deviations from the Brief

1. **No `ZeroRegister` instruction needed for MSD seed initialisation.**
   In `ApplySeries`, an empty coefficient vector (`ApplySeries{{}}`) scales the active seed vector by zero (`linalg::scal(0.0, phi)`), matching the `linalg::scal(0, Phi)` pre-step of `MomentsTD::Evolve`. Because `MomentEngine::Run` initialises `seed[Aux] = |r>` at the start of each state, `msd_plan` prepares the initial MSD seeds (`seed[Ket] = 0`, `seed[Aux] = |r>`) with `ApplySeries([])` followed by `MirrorToBra()` without adding a separate `ZeroRegister` instruction.
2. **Optional `scale` scalar on `ApplySeries`, `ApplyDerivedSeries`, `EvolveRegister`, and `EvolveDerived`.**
   Legacy `MomentsTD::Evolve` and `MomentsTD::Conv_Evolve` accumulate `sum_n c_n T_n(H~) v` with pure real/imaginary Bessel coefficients `c_n` and apply the global band-center phase `exp0 = exp(-i BandCenter dt / hbar)` once at the end via `linalg::scal(exp0, Phi)`. Carrying `scale` (default `(1.0, 0.0)`) on the series instructions allows Python to pass `exp0` so the final `linalg::scal(scale, phi)` happens after the Chebyshev sum without introducing physical units or time parameters into `lsquant/core/`.
3. **G6 `EF` energies in `tests/dev/stage_4_check.py`.**
   The brief's text in section 2 listed `0, 0.5, -1.5, 1.9` while `goldens/stage0/PROVENANCE.md` (line 77) and the committed files in `goldens/stage0/chain1d_ballistic/` have `EF = 0.0, 0.5, 1.0, 1.5` (`meanCorrelationVXchain1dKPM_M256_state20exactEF{0.000000,0.500000,1.000000,1.500000}JACKSON.dat`). `tests/dev/stage_4_check.py` checks against the four committed G6 golden files `(0.0, 0.5, 1.0, 1.5)`.

---

## 7. Open Questions

None.
