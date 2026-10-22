# Stage 8 developer report

Editing environment only: no builds, executions, profiling or gates were run,
as explicitly required by the task. No performance numbers have been measured.
Tester measurements and Fable review remain outstanding. No push is requested.

## Item 1 — profiling

Added `bench/profile_run.py` and `bench/README.md`. The driver profiles execution
of one validated run file in a child process, prints the top 20 cumulative entries,
records attributed `_core`/pybind11 self time, samples Linux thread count and ps
CPU percent, and optionally replays under perf without cProfile. It records only
its own child PIDs. Expected effect: identify Python, engine and occupancy costs;
the profiler itself is not a performance improvement.

The tester runs the commands in `bench/README.md` for G5 and Bi at M=128, R=2
before and after every item. Profile artifacts and measurements: **pending**.
The C++ fraction is explicitly an attribution from CPython builtin labels;
inspect the recorded labels when interpreting it.

## Scope and validation

Item 2 and the stage-6b files are owned by the other developer run; `bases.py`
is excluded here. No golden, gate harness, figure, legacy file or dependency is
changed. Build commands and numerical comparisons will be supplied to the tester;
there is no build or check output to quote from this editing environment.

## Item 3 — build profiles

`CMakeLists.txt`, `lsquant/CMakeLists.txt` and `build_profile.hpp` select golden
(default, no added flags/defines) or native (`-O3 -march=native -mtune=native
-fno-math-errno`, native-only macro). Bindings expose the compiled profile and
moment provenance. `runfile.py` writes it in receipts and rejects mismatched
resolved receipts. `chebmom_file.cpp` tags native on the existing first header
line of all four layouts and reads the optional tag; historical untagged headers
mean golden. Golden output bytes, precision and numerical kernels are retained.

Expected effect: native instruction selection and vectorization on Aron. Tester:
archive actual compile commands for gcc 9.4 and 11.4, build both profiles, compare
golden flag lists against the base, rerun gate 7 on the fat node, round-trip header
provenance, and repeat item-1 profiles. All results: **pending**.

Deviation requiring review: the brief simultaneously asks for provenance in every
header, no format change and unchanged golden bytes. Native uses a trailing first
line comment; golden provenance is implicit in the unchanged historical header.
This is a backward-readable header extension for native only, not a new numeric
layout. Older LSQUANT readers need the new reader for the generic native header.

## Item 4 — threads and memory placement

`runtime.py` sets absent OpenMP environment values before `_core` and NumPy import.
It counts physical cores within Linux process affinity. Explicit caller values win.
Native `memory_placement.hpp` initializes owned large vectors with static first
touch; on Linux it discards full interior pages of fresh allocations before the
fill, avoiding the serial value-initialization trap of `std::vector<complex>`.
Failure/unavailability of madvise affects placement only. Native CSR builders and
transformations copy final storage by static rows; block storage by static blocks.
Chebyshev recursion, vector batches and engine work vectors use the same row fill.
Golden numerical paths retain their original allocations and schedules.

Expected effect: fewer remote-memory accesses and sensible default occupancy.
Tester: repeat G5 and Bi profiles, inspect `/proc/PID/numa_maps` and `numastat -p`,
compare 24/48 threads on one socket and 96 on two, and two independent processes
bound one per socket. Exact commands and throughput interpretation are in
`bench/README.md`. All measurements: **pending**. Python defaults cannot change an
OpenMP runtime that another imported extension has already initialized; launch
with explicit environment variables in that case. A caller that changes thread
counts after construction forfeits placement, but not deterministic arithmetic.

Scope limitation: native placement targets final owned operator and engine buffers;
temporary Eigen assembly storage and Python-owned input arrays retain their library
allocation policy. This avoids changing the excluded 6b adapters or NumPy itself.

## Item 5 — native kernels and optional FFT contraction

- CSR: native whole-row nnz-balanced static partitions, tunable parallel crossover
  and multi-vector sweeps. Stored column reduction order is retained. First-touch
  CSR placement now follows that same partition. `work_partition.hpp` documents
  the tuning parameters; **the default crossover is a candidate, not measured**.
- Live vectors: `LinearOperator::MultiplyBatch` provides an adapter fallback;
  CSR and block operators override it. Native `ChebyshevVectors::ApplyToAll`
  dispatches a whole window, using a reusable output batch. Golden keeps its
  original per-vector applications. Existing composite adapters use the fallback.
- Dense blocks: native single-vector Eigen GEMV and per-block GEMM for live
  vectors, with static contiguous block ownership and an L2-sized microbatch.
- Dot table: native rank-2 engine uses `linalg::MomentTable`, packing inputs and
  computing fixed 16x16 result tiles by GEMM. Full inner products stay with one
  tile owner. Packing and scratch allocation are included in benchmark timing.
- Rank-2 reconstruction: native coefficient contraction computes `mu B` by GEMM
  in fixed 32-point groups, then the diagonal contraction with A. Golden keeps
  the existing scalar reduction and static schedule.
- Optional `FftMomentContraction`: collapses the cosine/exponential double series
  into frequency sums in fixed m,n order and evaluates four single-thread FFTW
  transforms. It returns polynomial sums only; Python owns damping, Green factors,
  units and integration. `kubo_fft_grid` selects it explicitly by argument or
  `LSQUANT_FFT_BACKEND=fft`; default remains table contraction. FFTW planner access
  is serialized with the existing transform operator through one shared mutex.
  This generic Fourier polynomial evaluator lives in `core/algebra/`; no FFTW
  dependency enters `core/chebyshev/` or I/O. Tuning environment variables are
  translated to numeric settings by Python (or the benchmark), never read by core.
- `bench/kernel_bench.cpp` and optional CMake target `stage8_kernels` measure
  CSR/vector batches, W=72 blocks, dot tables and coefficient contraction, with
  output hashes and explicitly modeled FLOP/byte rates.

Expected effects: less OpenMP launch overhead on small CSR matrices, balanced row
work, matrix reuse across vectors, and vectorized dense products. FFT reduces
polynomial evaluation complexity to quadratic coefficient folding plus FFTs.
No speedup or numerical result has been measured. Tester commands, crossover and
batch sweeps, byte comparisons across thread counts and stream models are in
`bench/README.md`. All before/after profiles, residuals and rates: **pending**.

Deviations/open questions for review: (1) no measured crossover is possible in
this environment; the tester must choose it from the supplied sweep. (2) the L2
formula selects kernel microbatches, not the storage window: shrinking that window
would multiply right-side recursions in the current sectioned engine. If even one
row shard exceeds L2, streaming is unavoidable. (3) fixed GEMM tiles and FFTs
reassociate sums relative to scalar goldens; native data is profile-tagged, and
optional FFT equality to historical golden curves is an outstanding measurement,
not an asserted gate. (4) composite batch forwarding beyond the generic fallback
is left to a follow-up after profiling; no excluded stage-6b adapter is edited.
