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
