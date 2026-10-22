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
