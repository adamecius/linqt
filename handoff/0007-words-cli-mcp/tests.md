# Handoff 0007 — tester report (stage 7: words, run files, CLI, MCP)

Tester: Opus subagent driven by Fable, on TCN (`fatnode-tacn`), worktree `~/lsquant-unified-0007`,
build `~/lsquant-unified-build7`, branch `handoff/0007-codex`.

**Head tested: `9f8ccc8`** (`origin/handoff/0007-agi`, review-1 items 1 to 3), fast-forwarded from
`origin/cpp-unified` 7510ef0. A first gate on `6ec46f1` was stopped part-way (section 9) and is not
reported except where it localised a defect. No C++ changed between 6ec46f1 and 9f8ccc8
(`git diff --stat 6ec46f1 9f8ccc8 -- lsquant CMakeLists.txt` is empty), so the build from 6ec46f1 is
the build of 9f8ccc8.

## 0. Summary

- **Quick closure gate (`S7_SET=quick`, 16 threads, the subset Fable asked for): 92 rows, 13 failures.**
  Every golden that a run file produced is byte-identical: G1 1/1, G2 1/1, G3 3/3, G4 6/6, G4c 2/2,
  G5 6/6, G7 3/3 (g7_sx, g7_sy, g7_sz), tutorial 01 24/24, tutorial 03 9/9, tutorial 04 4/4,
  G9 graphene real 8/8 (against `stage5_bases` with the same arguments), MCP 7/7.
  The 13 failures come from two run-time defects:
  - **D1**: every G10 run file fails (`'_core.ExactInverseOverlap' object is not callable`).
  - **D2**: every G9 exact-trace run file fails (`g9_*_eigen`, `g9_*_bloch_exact`: `int() argument ... 'NoneType'`).
  Neither run file writes a moment file. Both defects are localised to a call site, and the probes in
  section 5 show that with those call sites changed (in a scratch copy) the outputs are byte-identical.
- **Full gate (`S7_SET=all`, 36 run files, 16 and 1 threads) on 9f8ccc8, run before the change of plan:**
  - G1 to G7 and tutorials 01 to 05 are identical at both thread counts, every file.
  - G9 is identical for the 8 random-phase run files; D2 stops the 4 exact ones.
  - D1 stops all 8 G10 files.
  - Table in section 3.
- **MCP**: g1 and g5 through `run` and `status`. 7/7 byte-identical. `run` returned in 0.14 s and 0.007 s.
  A model dir outside `--root` is refused (`isError=True`), and no output dir is created.
- **Refusals through the CLI**: 6/6 refused, each with exit 1, a message on stderr and no output dir.
  The six are the five of the self-check plus a nested unknown key; both unknown-key files point
  `model.dir` at a missing path, which shows the refusal comes before any load.
- **Smells**: section 8. The serious ones:
  - (S1) `words describe` accepts `Tr[VX delta(E-H) VY G+(E)]`, `... G-(E)` and `... dG+(E)/dE`. It
    approximates all three as Kubo-Bastin instead of refusing them.
  - (S2) tutorials 01, 02 and 05 are generator run files. Their word, expansion and params are
    validated, then ignored: the dispatch calls `run_tutorial0N(out_dir)`.
- **Codex section 2.5** (deleted legacy tree) was **skipped** as instructed. The developer has not
  deleted the tree yet.
- **Incident:** while stopping my own gate at about 11:25 UTC I killed the parent shell of the
  stage-6 tester's gate, pid 3921966 (`check_goldens.sh 6 build6`). It is my error, not the other
  session's; details in section 9.

## 1. Build and self-check

Fresh build (`rm -rf ~/lsquant-unified-build7`), stage-2 cmake line; FFTW found through pkg-config:

```
$ cmake -S . -B ~/lsquant-unified-build7 -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR="$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')"
cmake exit 0      (CMakeCache: pkgcfg_lib_FFTW3_fftw3=/usr/lib/x86_64-linux-gnu/libfftw3.so)
$ cmake --build ~/lsquant-unified-build7 -j 32
build exit 0      [100%] Built target _core
```
Warnings: 101, 0 under `lsquant/`. `mcp` was already installed (2.3.0).

`OMP_NUM_THREADS=16 KPM_SEED=12345 PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_7_check.py`
on 9f8ccc8: `STAGE 7 SELF-CHECK PASS`, exit 0. Every line of (a) to (d) gives the same verdicts as
report.md section 6, and it also passed on 6ec46f1. The self-check covers only g1 and g5, so it
cannot see D1 or D2.

## 2. Harness (`scripts/check_goldens.sh 7`)

New stage-7 section. `S7_SET=all` (default) runs every run file at `S7_THREADS` (default `16 1`).
`S7_SET=quick` runs the closure subset at 16 threads. Stage `7aux` runs only the MCP and refusal
parts. `S7_KEEP=<dir>` keeps the outputs for the figure. Comparisons:

| old harness rows | now driven by | criterion |
|---|---|---|
| stage 1 G1, G2 (`stage1_moments`) | g1.json, g2.json | cmp vs goldens/stage0/chain1d |
| stage 2 G3 (+G1) | g3.json: G3 curve, `DOS_reconstructed.dat`, G1 | cmp |
| stage 3 G5 (2 moments + 4 curves), G4 (moments + 5 curves), G4c (moments + Bastin) | g5.json, g4.json, g4c.json | cmp |
| stage 4 G6 (moments + 4 curves), G7 (3 moments) | g6.json; g7.json and g7_s{x,y,z}.json | cmp |
| stage 5 G9 (no byte golden) | 12 `g9_*.json` | cmp vs `stage5_bases` run with the run file's arguments, plus `figures/stage_5_table.py` on the run-file outputs (full set only) |
| stage 6 G10 (`stage6_nonorth`) | 8 `g10_*.json` | cmp vs `stage6_nonorth` with the run file's arguments (exact applier, run once, at 16 threads, and only if the run file produced moments; it takes about 45 min per word on loaded TCN), plus max rel diff vs G1/G2/G5 <= 1e-10 (stage-6 brief section 6) |
| stage-0 tutorial captures | tutorial01..05.json | every non-PNG file of `goldens/stage0/tutorials/0N` at the same relative path, cmp; `captured/index.tsv` (the stage-0 capture shim's index: name, sha1, legacy temp path) is checked by sha1 of each produced capture |

The following old rows are no longer compared, because no run file produces them:
- the stage-1 batch-16 and random-phase repeat;
- stage-2 form (a);
- the stage-3/4 rows that reconstruct "curves from the golden moment file";
- the stage-5 convention rows (s5_conv.py).

The greps of stages 1 to 5 still run. Two stage-7 greps are new:
- `catalogue|runfile|grammar|Tr\[|mcp` on code lines of `lsquant/core` and `lsquant/io`;
- `eval(|exec(|literal_eval|exec_module` in `python/lsquant/{words,runfile.py,cli.py,mcp_server.py}`.

Decision: G9 real-basis run files carry no `csr_dir`, so they build the real basis in lsquant. They
are therefore compared with `stage5_bases` without `--csr-dir`, which matches the run file. The
stage-5 gate used the wannier2sparse CSR; see S6.

## 3. Per-golden table

### 3a. Quick closure gate (S7_SET=quick, T = 16, head 9f8ccc8, `GATE FAIL`, 92 rows, 13 failures)

| golden | run file(s) | identical | differ | wall (s) |
|---|---|---:|---:|---|
| G1 | g1 | 1 | 0 | 4.5 |
| G2 | g2 | 1 | 0 | 98.6 |
| G3 | g3 | 3 | 0 | 4.4 |
| G4 | g4 | 6 | 0 | 1.3 |
| G4c | g4c | 2 | 0 | 9.8 |
| G5 | g5 | 6 | 0 | 85.8 |
| G7 | g7_sx, g7_sy, g7_sz | 3 | 0 | 51.6, 29.8, 50.2 |
| Tut 01 | tutorial01 | 24 (23 cmp + index.tsv 22 sha1) | 0 | 86.5 |
| Tut 03 | tutorial03 | 9 (8 cmp + index.tsv 4 sha1) | 0 | 42.0 |
| Tut 04 | tutorial04 | 4 (3 cmp + index.tsv 3 sha1) | 0 | 49.0 |
| G9 | g9_graphene_real vs stage5_bases | 8 | 0 | 45.4 |
| G9 | g9_graphene_eigen vs stage5_bases | 0 | 8 (not produced, D2) | 0.8 (exit 1) |
| G10 | g10_haldane_real_generalised vs stage6_nonorth | 0 | 1 (no moment file, D1) | 1.1 (exit 1) |
| G10 vs G5 (1e-10) | g10_haldane_real_generalised | 0 | 2 (not produced) | |
| MCP | g1, g5 through the MCP server | 7 | 0 | |
| greps | stages 1 to 5 and 7 | all empty on code lines | | |

**Deferred to the final validation run** (S7_SET=all): g6, g7 (combined), g9_graphene_{bloch,
bloch_exact, bloch_w1, real_w1}, all g9_haldane_*, g10 other than haldane_real_generalised,
tutorial02 and tutorial05. They were in fact already run once on 9f8ccc8 (3b) before the change of
plan.

### 3b. Full gate (S7_SET=all, T = 16 and T = 1, head 9f8ccc8, 818 rows)

This run used the harness before two corrections:
- `index.tsv` was then cmp'd as a data file, so each tutorial reported one "DIFFERS (not produced)".
  A sha1 check on the kept outputs (`~/s7keep/T16`, `T1`) matches every entry: tutorial 01 22/22,
  02 192/192, 03 4/4, 04 3/3, 05 8/8, at both thread counts. The table counts those rows as identical.
- The stage-7 grep matched a comment line (`state_generator.hpp:232 ... Tr[...]`). The corrected grep
  excludes comment lines and is empty.

| golden | T=16 identical / differ | T=1 identical / differ |
|---|---|---|
| G1 | 1 / 0 | 1 / 0 |
| G2 | 1 / 0 | 1 / 0 |
| G3 | 3 / 0 | 3 / 0 |
| G4 | 6 / 0 | 6 / 0 |
| G4c | 2 / 0 | 2 / 0 |
| G5 | 6 / 0 | 6 / 0 |
| G6 | 5 / 0 | 5 / 0 |
| G7 (g7 + g7_s{x,y,z}) | 6 / 0 | 6 / 0 |
| Tut 01 | 24 / 0 | 24 / 0 |
| Tut 02 | 194 / 0 | 194 / 0 |
| Tut 03 | 9 / 0 | 9 / 0 |
| Tut 04 | 4 / 0 | 4 / 0 |
| Tut 05 | 17 / 0 | 17 / 0 |
| G9 vs stage5_bases | 64 / 32 | 64 / 32 |
| G9 stage-5 table | 0 / 4 | 0 / 4 |
| G10 vs stage6_nonorth | 0 / 8 | 0 / 8 |
| G10 vs G1/G2/G5 (1e-10) | 0 / 16 | 0 / 16 |
| MCP | 7 / 0 | |

- The 64 identical G9 rows are the 8 random-phase run files (real, bloch, real_w1 and bloch_w1 for
  both models), 8 files each. They match `stage5_bases` byte for byte, at both thread counts.
- The 32 differing rows are the 4 exact-trace run files (D2), whose files are not produced.
- The stage-5 table refuses to evaluate when the blochx/eigen outputs are missing, so its 4 FAIL rows
  are D2 as well.
- Wall times at T = 16 / T = 1 for the long ones: g6 84.5 / 23.4 s, tutorial02 85.3 / 94.9 s,
  tutorial05 51.4 / 209.9 s. TCN load was 40 or more throughout.

## 4. MCP transcript (9f8ccc8, `python -m lsquant.mcp_server --root $TMP/s7mcp`, a scratch root holding copies of the chain1d and haldane golden operators; client `mcp` 2.3.0 over stdio)

```
MCP <- list_tools: ["words_describe", "run", "status", "reconstruct", "inspect"]
MCP -> call_tool: {"tool": "run", "arguments": {"runfile_json": "{\"model\": {\"kind\": \"csr\", \"dir\": \"goldens/stage0/chain1d/operators\", \"label\": \"chain1d\"}, \"basis\": {\"name\": \"real\"}, \"metric\": {\"policy\": \"orthogonal\"}, \"word\": \"Tr[delta(E-H)]\", \"expansion\": {\"moments\": 64, \"window\": [-2.0, 2.0], \"trace\": \"exact\", \"states\": 2048}, \"output\": {\"dir\": \"out/g1\", \"legacy_names\": true, \"legacy_layout\": true}}"}}
MCP <- run (isError=False, returned in 0.140 s): {"run_id": "f9f95d45d006", "log": "$TMP/s7mcp/out/g1/run_f9f95d45d006.log", "output_dir": "$TMP/s7mcp/out/g1"}
MCP -> call_tool: {"tool": "status", "arguments": {"run_id": "f9f95d45d006"}}
MCP <- status (isError=False, returned in 0.009 s): {"run_id": "f9f95d45d006", "status": "running", "returncode": null, ...}
MCP <- status poll 4 (final): {"run_id": "f9f95d45d006", "status": "succeeded", "returncode": 0, ...}
MCP -> call_tool: {"tool": "run", "arguments": {"runfile_json": "{\"model\": {\"kind\": \"csr\", \"dir\": \"goldens/stage0/haldane/operators\", \"label\": \"haldane\", \"operators\": [\"VX\", \"VY\"]}, ..., \"word\": \"Tr[VX delta(E-H) VY dG-(E)/dE]\", \"expansion\": {\"moments\": 128, \"window\": [-4.5, 4.5], \"trace\": \"exact\", \"states\": 512}, \"reconstruction\": {\"kernel\": \"jackson\", \"broadening_meV\": 10.0, \"formulas\": [\"bastin\", \"sea\", \"surf\", \"greenwood\"], \"legacy_text_round_trip\": true}, \"output\": {\"dir\": \"out/g5\", ...}}"}}
MCP <- run (isError=False, returned in 0.007 s): {"run_id": "53d8a156c55c", "log": "$TMP/s7mcp/out/g5/run_53d8a156c55c.log", "output_dir": "$TMP/s7mcp/out/g5"}
MCP <- status (isError=False, returned in 0.004 s): {"run_id": "53d8a156c55c", "status": "running", "returncode": null, ...}
MCP <- status poll 127 (final): {"run_id": "53d8a156c55c", "status": "succeeded", "returncode": 0, ...}
MCP -> call_tool: {"tool": "run", "arguments": {"runfile_json": "{\"model\": {\"kind\": \"csr\", \"dir\": \"/etc\", ...}, \"output\": {\"dir\": \"out/outside\", ...}}"}}
MCP <- run (isError=True, returned in 0.792 s): "Error executing tool run"
MCP RESULT confinement (model.dir=/etc): refused (isError=True); out dir created: False
IDENTICAL [MCP] g1.json -> chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [MCP] g5.json -> haldane/NonEqOpVX-{VY,VX}haldaneKPM_M128x128_stateexact.chebmom2D (G5)  (2 files)
IDENTICAL [MCP] g5.json -> haldane/Kubo{Bastin,BastinSea,BastinSurf,GreenWood}_...JACKSON.dat (G5)  (4 files)
```
The full JSON, unabridged, is in the gate logs `~/s7logs/gate7.log` and `gate7_quick.log` on TCN.

On the refused call:
- the agent gets only the generic "Error executing tool run";
- the reason (`ValueError: Access denied: path '/etc' ... outside configured root`) goes to the
  server's stderr as an `UnexpectedToolError` traceback.

On 6ec46f1 the same call came back with `isError=False`. 9f8ccc8 corrected the flag, but not the
message (S7).

## 5. Refusals through the CLI (`python -m lsquant run <bad>.json --out <dir>`)

| case | run file | exit | stderr | output dir |
|---|---|---|---|---|
| unknown operator | g1 + `Tr[FOO delta(E-H)]` | 1 | `Unknown operator 'FOO' in word ... (registered operators: ['1', 'A', 'B', 'O', 'P', 'PX', ... 'VZSZ'])` | not created |
| rank 3 | g2 + `Tr[VX delta VY delta VZ delta]` | 1 | `Unsupported word ... of rank 3 ...; nearest recognised shape is 'Tr[VX delta(E-H) VX delta(E-H)]' (kubo_greenwood)` | not created |
| mixed basis | g1 + `Tr[VX@bloch delta(E-H)]` | 1 | `Mixed basis without explicit transform: operator 'VX@bloch' is in 'bloch' basis while calculation/active basis is 'real' ...` | not created |
| U(t) without time | g6 without `time` + `Tr[VX delta(E-H) U(t)^+ VX U(t)]` | 1 | `Word ... contains time evolution U(t) without 'time' parameters (steps and tmax required)` | not created |
| unknown key (top) | g1 + `unknown_extra_key`, `model.dir` = missing path | 1 | `Unknown key 'unknown_extra_key' in run file (top-level); allowed keys are [...]` | not created |
| unknown key (nested) | g1 + `expansion.windw`, `model.dir` = missing path | 1 | `Unknown key 'windw' in run file (expansion); allowed keys are [...]` | not created |

Both unknown-key run files point `model.dir` at a path that does not exist. The message is still the
key error, which shows that validation precedes any load.

`words describe` without a run file (reported, not gated):
- `FOO`, rank 3 and `@bloch` are refused, with the same messages as above.
- `Tr[VX G+(E) VY delta(E-H)]` is refused with a nearest-shape message.
- `Tr[VX delta(E-H) U(t)^+ VX U(t)]` is accepted (describe uses default time parameters, by design).
- **`Tr[VX delta(E-H) VY G+(E)]`, `Tr[VX delta(E-H) VY G-(E)]` and `Tr[VX delta(E-H) VY dG+(E)/dE]`
  are ACCEPTED as `kubo_bastin`** (formula `bastin`, engine word `VX delta VY delta`): smell S1.

## 6. Defects (localised, not fixed)

**D1: every G10 run file fails.** 8/8 run files, both routes, at T = 16 and 1. The errors are
`'_core.ExactInverseOverlap' object is not callable` (generalised) and
`'_core.ExactInverseSqrtOverlap' object is not callable` (Löwdin).

- **Cause.** `lsquant/api/bindings.cpp:263,265,267` binds `overlap`, `inverse_overlap` and
  `inverse_sqrt_overlap` of `MetricPolicy` as read-only **properties**, but the Python side calls
  them as methods:
  - `python/lsquant/runfile.py:1158` (`metric_obj.inverse_sqrt_overlap()`) and `:1176`
    (`metric_obj.inverse_overlap()`) on 9f8ccc8;
  - `python/lsquant/words/compiler.py:484, 490, 519, 523, 539, 542, 547, 579, 586`
    (`inverse_*_overlap()`) and `:524, 548` (`overlap()`).
- `python/lsquant/_core.pyi:93,95` declares them as methods, so the stub is out of sync with the binding.
- The developer's report lists G10 run files, but the self-check never runs them.
- **Probe.** I made a scratch copy of `python/lsquant` (`/tmp/s7probe*`; the repository is unchanged)
  with only those call sites changed from `x()` to `x`, run at T = 16 against `stage6_nonorth` with
  the run file's arguments:
  - `g10_chain1d_real_lowdin`: 4/4 files byte-identical (.chebmom1D, .chebmom2D, both .chebmom);
  - `g10_chain1d_{real,complex}_generalised` and `g10_haldane_{real,complex}_generalised`: 16/16
    byte-identical.
  - The other Löwdin probes were stopped for time (dense S^-1/2, about 45 min per word on loaded TCN)
    and belong to the final run.
  - With only `inverse_*` changed, the generalised Kubo word still failed at `metric_policy_obj.overlap()`
    (`'_core.SparseMatrixOperator' object is not callable`).

**D2: every G9 exact-trace run file fails.** `g9_{graphene,haldane}_{eigen,bloch_exact}`, T = 16
and 1, with `int() argument must be a string, a bytes-like object or a real number, not 'NoneType'`.

- **Cause.** These run files omit `expansion.states`, so `RunFile` resolves it to `None`. The preview
  compile at `runfile.py` (the `compile_preview` loop, `execute_runfile`) passes `rf.expansion`
  unchanged. `words/compiler.py:401`, `int(exp_cfg.get("states", cfg.get("states", 1)))`, then gets
  `None`, because the key is present.
- **Probe.** With line 401 changed to treat `None` as absent, in a scratch copy, the 4 run files
  against `stage5_bases` with the same arguments give `g9_graphene_{bloch_exact,eigen}` 16/16
  and `g9_haldane_{bloch_exact,eigen}` 32/32 byte-identical in all.

## 7. Figure

`figures/out/stage_7/fig_stage_7.png` (generated with `python3 figures/stage_7.py --runs ~/s7keep/T16`).
The data are the full-run outputs on 9f8ccc8, which include tutorials 02 and 05 because the full set
was already run. Caption (`fig_stage_7.caption.txt`):

> FIG. 7. The six tutorial figures of stage 0 (left, legacy build, figures/out/stage_0) beside the same
> curves drawn from the outputs of `python -m lsquant run runs/goldens/tutorial0N.json` (middle, solid)
> over the stage-0 captures (dashed), with the residual new - golden (right). OMP_NUM_THREADS = 16,
> KPM_SEED = 12345.
> (a) tutorial 01 from runs/goldens/tutorial01.json: 24 of 24 golden files (every non-PNG file;
> captured/index.tsv checked by sha1) byte-identical; max |new - golden| 0.
> (b) tutorial 02 from runs/goldens/tutorial02.json: 194 of 194; 0.
> (c) tutorial 03 from runs/goldens/tutorial03.json: 9 of 9; 0.
> (d) tutorial 04 from runs/goldens/tutorial04.json: 4 of 4; 0.
> (e) tutorial 05 from runs/goldens/tutorial05.json: 17 of 17; 0.
> (f) tutorial 06: stage-0 figure only; there is no runs/goldens/tutorial06.json, so tutorial 06 is not
> reproduced from words.

The middle panels show the curves as written in the files (raw legacy units, e.g. sigma_xy before
the tutorial's normalisation by the plateau). Panels 01, 03 and 04 are the closure set. Panels 02
and 05 come from the same full run and are also identical.

## 8. Smells

- **S1 (a word approximated instead of refused).** `match_catalogue`
  (`words/catalogue.py`, rank-2 branch) maps any second spectral slot in
  `{dG-/dE, dG+/dE, G-, G+}` to `kubo_bastin`, with the same engine word `VX delta VY delta` and
  formula `bastin`.
  - `Tr[VX delta(E-H) VY G+(E)]` and `Tr[VX delta(E-H) VY G-(E)]` are not the Bastin integrand, yet
    they compile to it silently.
  - `dG+/dE` is the conjugate term, also folded in without comment.
  - The brief says a word the catalogue does not recognise is refused, never approximated.
- **S2 (words that do nothing).** Generator run files (tutorial01, tutorial02, tutorial05) validate
  `word`, `expansion`, `time`, `reconstruction` and `model.params`, then
  `_execute_generator_runfile` calls `run_tutorial01(out_dir)` / `run_tutorial02(out_dir)` /
  `run_tutorial05(out_dir)` (runfile.py on 9f8ccc8, lines 767, 771, 792). Every parameter is
  hard-coded in `apps/tutorial0N.py`.
  - Changing the run file's word or moments changes nothing.
  - These goldens are therefore reproduced by the tutorial apps, not "from words".
  - Review-1 item 1 replaced `script` by a generator name; the dispatch itself is unchanged.
- **S3 (an implicit extra word).** When a Bastin `VX-VY` word asks for `greenwood`,
  `execute_runfile` silently appends the companion word `Tr[VX delta(E-H) VX delta(E-H)]`
  (runfile.py, comment "append the companion longitudinal word"). This is how g5 and tutorial03 get
  `KuboGreenWood_NonEqOpVX-VX...`.
- **S4 (a value approximated instead of refused).** For `model.kind = csr`, `expansion.window = "eigen"`
  falls through to Gershgorin (runfile.py, `else: window = window_from_gershgorin(ham)`), with no
  error and no record in `run.resolved.json`.
- **S5 (stub out of sync).** `_core.pyi` declares `overlap`, `inverse_overlap` and
  `inverse_sqrt_overlap` as methods while the binding exposes properties; this is the root of D1.
- **S6 (G9 real basis).** The `g9_*_real*` run files build the real basis in lsquant (no `csr_dir`);
  the stage-5 gate used the wannier2sparse CSR (`--csr-dir`). Bytes match `stage5_bases` without
  `--csr-dir`, but the configuration the stage-5 gate certified is not the one the run files use.
- **S7 (MCP error reporting).** Tools raise a bare `ValueError`. Under mcp 2.3.0 the agent receives
  "Error executing tool run" without the reason, and the server prints a traceback on stderr. A
  `ToolError` would carry the message.
- **S8 (coverage).**
  - There is no run file for tutorial 06.
  - G10 run files exist only for the exact applier; the stage-6 gate also covered ischo and chebyshev.
  - The self-check runs only g1 and g5, which is why D1 and D2 reached the tester.
- **Hygiene.** `tutorial01.json` writes operators, BOUNDS and state files into the output dir
  (expected for a generator). No run wrote into the repository (`git status` clean apart from the
  tester's files).

Checks with nothing found:
- **eval:** no `eval`/`exec`/`literal_eval` in the word language, run file, CLI or MCP server.
  `exec_module` appears only in `apps/tutorial0{1,2,3,5}.py`, which import the example generators by
  path.
- **basis or word under core:** the stage-1 to stage-5 greps and the stage-7 grep are empty on code
  lines.
- **goldens:** nothing re-blessed, nothing under `goldens/` touched (the harness only reads it).
- **legacy tree:** not deleted, as instructed; codex 2.5 skipped.

## 9. Run history and incident

- **6ec46f1, 10:10 UTC.** A full gate started. The G10 `stage6_nonorth` references (dense Löwdin,
  about 45 min per word) made the T = 1 pass infeasible, so I stopped the run (11:25) and changed the
  harness to run those references once, and only when the run file produced moments.
- **Incident.** When I stopped my gate at 11:25 UTC I killed pid 3921966 by mistake. It was
  `bash scripts/check_goldens.sh 6 build6`, the stage-6 tester's gate. Its EXIT trap removed its
  temporary directory, so that run was lost.
  - The pkill that Fable was told about (about 11:26) was this kill of mine, not the other session's.
  - The stage-6 session has since restarted (`check_goldens.sh 6only build6`). I have not touched it again.
  - From then on I killed only processes identified by their full command line as mine.
- **11:25.** A restart on 6ec46f1 was stopped at 11:35 to move to 9f8ccc8, as Fable asked.
- **9f8ccc8, 11:35 to 12:21.** Self-check, then the full gate (`S7_SET=all`, T = 16 and 1), section 3b.
- **12:24 to about 12:45.** Quick closure gate (`S7_SET=quick`, T = 16) with the final harness,
  section 3a.

## 10. Files

- `scripts/check_goldens.sh` (stage 7, 7aux, S7_SET)
- `figures/stage_7.py`
- `figures/out/stage_7/fig_stage_7.png` and `fig_stage_7.caption.txt`
- this `tests.md`

Logs on TCN: `~/s7logs/{gate7.log, gate7_quick.log, selfcheck2.log, probe*.log}`. Kept outputs:
`~/s7keep` (full run, T16 and T1) and `~/s7keep_quick`.
