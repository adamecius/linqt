# Handoff 0001 Report (Agi) — Stage 1: Core Algebra and Chebyshev Engine

## 1. Files Added and Modified

### Added Files
- `lsquant/CMakeLists.txt`
- `lsquant/core/algebra/types.hpp`
- `lsquant/core/algebra/linear_operator.hpp`
- `lsquant/core/algebra/sparse_matrix.hpp`
- `lsquant/core/algebra/sparse_matrix.cpp`
- `lsquant/core/algebra/composite_operators.hpp`
- `lsquant/core/algebra/composite_operators.cpp`
- `lsquant/core/algebra/linalg.hpp`
- `lsquant/core/algebra/linalg.cpp`
- `lsquant/core/algebra/state_generator.hpp`
- `lsquant/core/algebra/state_generator.cpp`
- `lsquant/core/algebra/spectral_bounds.hpp`
- `lsquant/core/algebra/spectral_bounds.cpp`
- `lsquant/core/chebyshev/spectral_window.hpp`
- `lsquant/core/chebyshev/rescaled_operator.hpp`
- `lsquant/core/chebyshev/rescaled_operator.cpp`
- `lsquant/core/chebyshev/chebyshev_recursion.hpp`
- `lsquant/core/chebyshev/chebyshev_recursion.cpp`
- `lsquant/core/chebyshev/chebyshev_vectors.hpp`
- `lsquant/core/chebyshev/chebyshev_vectors.cpp`
- `lsquant/core/chebyshev/chebyshev_moments.hpp`
- `lsquant/core/chebyshev/chebyshev_moments.cpp`
- `lsquant/core/chebyshev/plan.hpp`
- `lsquant/core/chebyshev/moment_engine.hpp`
- `lsquant/core/chebyshev/moment_engine.cpp`
- `lsquant/io/csr_text.hpp`
- `lsquant/io/csr_text.cpp`
- `lsquant/io/chebmom_file.hpp`
- `lsquant/io/chebmom_file.cpp`
- `lsquant/apps/stage1_moments.cpp`
- `tests/dev/stage_1_check.cpp`

### Modified Files
- `CMakeLists.txt` (`option(LSQUANT_BUILD_NEW_CORE "" ON)` and guarded `add_subdirectory(lsquant)`)

---

## 2. Build Commands and Output on TCN (`fatnode-tacn`)

### Commands
```bash
rsync -a --delete --exclude .git --exclude build /root/lsquant-agi/ jh-x:~/lsquant-agi-build/
ssh jh-x "cd ~/lsquant-agi-build && cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=ON && cmake --build build -j 16"
```

### Output
```text
-- The CXX compiler identification is GNU 11.4.0
-- Found OpenMP_CXX: -fopenmp (found version "4.5")
-- Found OpenMP: TRUE (found version "4.5")
-- Found PkgConfig: /usr/bin/pkg-config (found version "0.29.2")
-- Found GSL: /usr/include (found version "2.7.1")
-- Found Threads: TRUE
-- Checking for module 'fftw3'
--   Found fftw3, version 3.3.8
-- wannier2sparse: /home/ICN2/jgarcia/lsquant-agi-build/.cache/wannier2sparse/build/wannier2sparse
-- Configuring finished
-- Generating finished
-- Build files have been written to: /home/ICN2/jgarcia/lsquant-agi-build/build
[  1%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/algebra/sparse_matrix.cpp.o
[  2%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/algebra/linalg.cpp.o
[  3%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_vectors.cpp.o
[  5%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/algebra/state_generator.cpp.o
[  7%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/algebra/composite_operators.cpp.o
[  7%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/algebra/spectral_bounds.cpp.o
[  7%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/chebyshev/rescaled_operator.cpp.o
[  8%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_recursion.cpp.o
[ 18%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/chebyshev/moment_engine.cpp.o
[ 19%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/core/chebyshev/chebyshev_moments.cpp.o
[ 28%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/io/csr_text.cpp.o
[ 30%] Building CXX object lsquant/CMakeFiles/lsquant_core.dir/io/chebmom_file.cpp.o
[ 52%] Building CXX object lsquant/CMakeFiles/stage1_moments.dir/apps/stage1_moments.cpp.o
[ 57%] Building CXX object lsquant/CMakeFiles/stage_1_check.dir/__/tests/dev/stage_1_check.cpp.o
[ 58%] Built target stage1_moments
[ 64%] Built target stage_1_check
[100%] Built target lsquant
```
Compiler warnings under `lsquant/` and `tests/dev/stage_1_check.cpp` with `-Wall -Wextra`: 0.

---

## 3. `stage_1_check` Output at 1 and 8 Threads on TCN

### Commands
```bash
OMP_NUM_THREADS=1 build/lsquant/stage_1_check
OMP_NUM_THREADS=8 build/lsquant/stage_1_check
```

### Output (`OMP_NUM_THREADS=1`)
```text
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
```

### Output (`OMP_NUM_THREADS=8`)
```text
dos max_abs_error: 2.47025e-15
kubo max_abs_error: 1.13687e-13
OMP_NUM_THREADS=1 vs 8 bit-identical: YES
```

---

## 4. `stage1_moments` Output on `test/golden/chain1d/operators`

(`goldens/stage0` is not present in this branch yet; executed against `test/golden/chain1d/operators` with `--window -2 2 --trace exact --states 2048 --moments 64`.)

### Commands
```bash
OMP_NUM_THREADS=8 KPM_SEED=12345 build/lsquant/stage1_moments \
  --operators test/golden/chain1d/operators --label chain1d --window -2 2 \
  --trace exact --states 2048 --moments 64 --word dos --out /tmp/stage1_out

OMP_NUM_THREADS=8 KPM_SEED=12345 build/lsquant/stage1_moments \
  --operators test/golden/chain1d/operators --label chain1d --window -2 2 \
  --trace exact --states 2048 --moments 64 --word kubo --out /tmp/stage1_out
```

### Standard Output
```text
stage1_moments: word=dos dim=2048 states=2048 (exact) M=64 -> /tmp/stage1_out/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
stage1_moments: word=kubo dim=2048 states=2048 (exact) M=64x64 -> /tmp/stage1_out/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
```

### Leading Lines of Produced Files
```text
==> /tmp/stage1_out/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D <==
2048 4 0
64
1 0
0 0
0 0
0 0
0 0
0 0
0 0
0 0

==> /tmp/stage1_out/SpectralOp1chain1dKPM_M64_stateexact.chebmom <==
chebmom 1
rank 1
dims 64
window -2 2
system_size 2048
states exact 2048
word 1 delta
kpm_weighted 1
symmetrised 0
data

==> /tmp/stage1_out/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D <==
2048 4  0
64 64 
2 -0
0 0
-2 0
0 0
0 0
0 0
0 0
0 0

==> /tmp/stage1_out/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom <==
chebmom 1
rank 2
dims 64 64
window -2 2
system_size 2048
states exact 2048
word VX delta VX delta
kpm_weighted 1
symmetrised 1
data
```

---

## 5. Legacy Numerical Conventions Reproduced and Deviations

### Legacy Numerical Conventions Reproduced for Byte Identity
1. **Hermiticity Fold in `SparseMatrix::FromCsr` (`0.5 * (transp + matrix_)`):** Legacy `SparseMatrixType::ConvertFromCSR` folds Hermitian inputs via `matrix_ = 0.5 * (transp + matrix_)`. In `.CSR` files such as `chain1d.VX.CSR`, off-diagonal entries contain `-0` real components (`-0 -1`) whose conjugate sum with `+0 1` converts `-0.0` into `+0.0`. `SparseMatrix::FromCsr` and `SparseMatrix::FromTriplets` apply this exact Eigen expression whenever `||M - M^+|| / ||M|| <= 1e-8`, preserving non-Hermitian matrices untouched when `||M - M^+|| / ||M|| > 1e-8`.
2. **Real Scalar Weight Multiplication in `MomentEngine::Run` (`Trace`):** Legacy `SpectralMoments` evaluates `scal * linalg::vdot(...)` where `scal` is a `double`. Multiplying `double * std::complex<double>` scales the real and imaginary parts independently (`(scal * re, scal * im)`), whereas multiplying `std::complex<double>(scal, 0.0) * std::complex<double>(re, im)` evaluates `scal * im + 0.0 * re`, turning `-0.0` imaginary parts into `+0.0`. `MomentEngine::Run` checks `weight.imag() == 0.0` and uses `weight.real() * d` so IEEE-754 signed zeros match legacy `vdot` scaling and raw accumulation (`2 -0` in `kubo`).
3. **Legacy `.chebmom2D` Header Whitespace:** `WriteLegacyChebmom2D` emits two spaces between `BandWidth` and `BandCenter` (`"2048 4  0"`) and a trailing space after each dimension on line 2 (`"64 64 "`), matching `chebyshev::Moments2D::saveIn`.

### Deviations from the Brief
1. **Legacy `lsquant` Target Output Directory in `lsquant/CMakeLists.txt`:** The root `CMakeLists.txt` sets `CMAKE_RUNTIME_OUTPUT_DIRECTORY` to `${CMAKE_BINARY_DIR}`, and `src/CMakeLists.txt` defines an executable target named `lsquant` whose default output path is `${CMAKE_BINARY_DIR}/lsquant`. Adding `add_subdirectory(lsquant)` creates the directory `${CMAKE_BINARY_DIR}/lsquant/`, causing `/usr/bin/ld: cannot open output file ../lsquant: Is a directory` when linking the legacy `lsquant` executable. Without touching `src/CMakeLists.txt`, `lsquant/CMakeLists.txt` sets `RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/legacy_bin"` on `TARGET lsquant` when `LSQUANT_BUILD_NEW_CORE` is enabled, while placing `stage1_moments` and `stage_1_check` in `${CMAKE_BINARY_DIR}/lsquant/`.
2. **Const Accessors on `ChebyshevMoments`:** In addition to `std::vector<Scalar>& Data()`, `const std::vector<Scalar>& Data() const` was added alongside the `const` overloads of `At(...)` so `WriteChebmom`, `WriteLegacyChebmom1D`, and `WriteLegacyChebmom2D` can take `const ChebyshevMoments&`.

---

## 6. Open Questions
- None for Stage 1.
