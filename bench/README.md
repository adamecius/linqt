# Stage 8 measurements on Aron

These drivers are measurement tools, not gates. Run from the checkout root in
`~/lsq-env`, with the selected build's `lsquant_core` directory and `python/` on
`PYTHONPATH`. Keep outputs outside the checkout. Never overwrite a baseline.

Item 1, before any kernel change (repeat after each numbered item):

```bash
export PYTHONPATH="$PWD/python:$PWD/build-golden/lsquant_core"
export KPM_SEED=12345 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_NUM_THREADS=48 OMP_PROC_BIND=close OMP_PLACES=cores
python bench/profile_run.py runs/goldens/g5.json --out /tmp/stage8/item1-g5 --perf
python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item1-bi --perf
```

`BI_RUN` must name the tester's 6b Bi block spin Hall run file, including its
local SIESTA input paths. Keep its trace mode, word, metric and reconstruction
unchanged; record its path, contents and input checksums with the artifacts.
The short profile changes only M and vector count in memory.

Record `profile.txt` (top 20 cumulative entries), `profile.json` (C++ attributed
self time and fraction), and `occupancy.jsonl` (thread counts and ps CPU percent,
which is the process lifetime average, not an instantaneous sample). Inspect the
attributed entry names: CPython labels pybind11 calls as `PyCapsule` builtins;
unlabelled extension calls are excluded, so this fraction is an attribution,
not a sampling profiler's measurement of all native work. Imports are excluded.
`perf.csv` is an additional unprofiled replay, with its own exit status; missing
perf or denied counters must be recorded, never treated as zero measurements.
Record compiler, dependency versions, affinity, host load and governor as well.

Item 3: configure separate clean directories with the same dependencies and build
type. These commands require only the packages already used by LSQUANT:

```bash
cmake -S . -B build-golden -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DLSQUANT_BUILD_PROFILE=golden -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-golden -j 24
cmake -S . -B build-native -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" -DLSQUANT_BUILD_PROFILE=native -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-native -j 24
export PYTHONPATH="$PWD/python:$PWD/build-native/lsquant_core"
```

Golden adds no flags/defines to the historical command lines. Native appends
`-O3 -march=native -mtune=native -fno-math-errno` and `LSQUANT_NATIVE_BUILD=1`.
Neither profile enables fast-math or unsafe-math optimizations. Record the actual
compile commands (including inherited Release and pybind11 LTO flags). Native
headers carry `# build_profile=native` at the end of the first line; the existing
golden headers are untagged and mean golden. No numeric layout or precision changes.
The reader retains source provenance when a moment file is rewritten by another
profile. Receipts always name the loaded engine profile and reject replay under a
different one. Compare moment data only within one profile and one machine.

Item 4: import LSQUANT before another library initializes OpenMP. Defaults fill
only absent variables: close binding, core places, physical cores in the current
CPU affinity (so `numactl`/scheduler cpusets constrain the count). Set BLAS threads
to one to avoid nested teams. A static contiguous block group stays with the same
thread for initialization and multiplication; keep the team size fixed per process.
First touch applies to owned native CSR, block and recursion/batch buffers; on
Linux complete interior pages of newly allocated buffers are discarded before
the parallel fill. It never discards borrowed data. Allocation reuse and policy
still matter: record `numastat -p PID` and `/proc/PID/numa_maps`, not just affinity.

```bash
for threads in 24 48; do
  OMP_NUM_THREADS=$threads numactl --cpunodebind=0 --membind=0 python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item4-bi-socket0-t$threads
done
OMP_NUM_THREADS=96 numactl --cpunodebind=0,1 --membind=0,1 python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item4-bi-dual-t96
OMP_NUM_THREADS=48 numactl --cpunodebind=0 --membind=0 python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item4-bi-p0 > /tmp/stage8/item4-p0.log 2>&1 &
echo $! > /tmp/stage8/item4-p0.pid
OMP_NUM_THREADS=48 numactl --cpunodebind=1 --membind=1 python bench/profile_run.py "$BI_RUN" --moments 128 --states 2 --out /tmp/stage8/item4-bi-p1 > /tmp/stage8/item4-p1.log 2>&1 &
echo $! > /tmp/stage8/item4-p1.pid
wait
```

The last recipe measures two independent jobs and aggregate throughput, not the
latency of one distributed trace. Keep their outputs separate; do not average
already normalized moment files to claim one exact trace. Compare the same total
work (two serial jobs versus these two concurrent jobs). Use the recorded child
PIDs in each output directory for numastat or stopping work, never a process name.

Item 5: enable the optional measurement executable in both build directories:

```bash
cmake -S . -B build-golden -DLSQUANT_BUILD_BENCHMARKS=ON
cmake --build build-golden -j 24
cmake -S . -B build-native -DLSQUANT_BUILD_BENCHMARKS=ON
cmake --build build-native -j 24
mkdir -p /tmp/stage8/kernels
for profile in golden native; do
  for threads in 1 24 48 96; do
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels spmv 1000000 1 100 > /tmp/stage8/kernels/$profile-spmv-t$threads.json
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels spmm 1000000 16 30 > /tmp/stage8/kernels/$profile-spmm-t$threads.json
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels block 2304 32 30 > /tmp/stage8/kernels/$profile-block-t$threads.json
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels dot 165888 128 3 > /tmp/stage8/kernels/$profile-dot-t$threads.json
    OMP_NUM_THREADS=$threads build-$profile/lsquant_core/stage8_kernels contract 256 7680 3 > /tmp/stage8/kernels/$profile-contract-t$threads.json
  done
done
for cutoff in 0 4096 16384 65536 262144; do
  for rows in 64 256 1024 4096 16384 65536; do
    OMP_NUM_THREADS=48 LSQUANT_SPMV_MIN_NNZ=$cutoff numactl --cpunodebind=0 --membind=0 build-native/lsquant_core/stage8_kernels spmv $rows 1 1000 > /tmp/stage8/kernels/cutoff-$cutoff-rows-$rows.json
  done
done
for batch in 1 2 4 8 16 32 64 128; do
  OMP_NUM_THREADS=48 LSQUANT_KERNEL_BATCH=$batch numactl --cpunodebind=0 --membind=0 build-native/lsquant_core/stage8_kernels block 2304 128 30 > /tmp/stage8/kernels/block-batch-$batch.json
  OMP_NUM_THREADS=48 LSQUANT_KERNEL_BATCH=$batch numactl --cpunodebind=0 --membind=0 build-native/lsquant_core/stage8_kernels spmm 1000000 128 10 > /tmp/stage8/kernels/csr-batch-$batch.json
done
```

Repeat at least three times with distinct output paths; report medians and range.
The cutoff default is an **unmeasured candidate** until this sweep establishes a
crossover on Aron. Hashes must agree between thread counts within a profile and
configuration; report any discrepancy, never bless a tolerance automatically.
Also compare the full production moment and curve files between thread counts.

Kernel microbatches and storage windows serve different purposes. For a CSR row
shard `q=ceil(N/T)`, the candidate microbatch is
`max(1,floor((512 KiB/(2*16)-3*q)/(2*q)))`, capped at the live vector count. If a
single stream exceeds L2, this returns one and the kernel streams from memory.
For block size W, the microbatch is
`max(1,floor((512 KiB/(2*16)-W*W)/(2*W)))`, likewise capped. The half-L2 budget
leaves space for Eigen packing, stack and other streams; it is a tuning starting
point. `LSQUANT_KERNEL_BATCH` overrides either for measurement. Chebyshev storage
windows keep the requested `expansion.batch` (default M), preventing repeated
right-side recursions merely to satisfy a cache estimate. The moment table packs
the stored window and evaluates fixed 16x16 output tiles by GEMM; contraction
uses fixed groups of 32 energy points. Geometry and reduction trees of each dense
product do not depend on the team size.

FFT is optional via `kubo_fft_grid(..., backend="fft")` or
`LSQUANT_FFT_BACKEND=fft` for `bastin_fft`/`greenwood_fft` reconstruction. The default
`tables` uses the historical scalar contraction in golden, GEMM in native. Compare
FFT and tables on **identical FFT nodes**; uniform-grid curves have a different
grid and cannot be directly subtracted. Frequency folding uses FFTW already linked
by the project, with single-thread transforms and deterministic coefficient sums.
It changes floating-point association, so equality to scalar goldens needs tester
evidence and reviewer acceptance under AGENTS.md, never a silent tolerance change.

Kernel `modeled_gflops` and `modeled_gbs` use the FLOP/byte counts printed by the
driver: complex multiply-add = 8 FLOP, CSR value/index = 24 bytes (64-bit index),
row pointer = 8 bytes, vector scalar = 16 bytes. They exclude setup and include
packing in timed GEMMs. Bytes are an algorithmic stream model, not measured DRAM
traffic (cache reuse/write allocation can change hardware traffic). For contraction
and dots they are lower-bound unique input/output streams; label them accordingly.
