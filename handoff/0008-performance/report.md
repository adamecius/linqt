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
