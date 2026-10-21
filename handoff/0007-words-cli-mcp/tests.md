# Handoff 0007 — tester report, closure run (stage 7: words, run files, CLI, MCP, deleted legacy tree)

**Tester:** Opus subagent driven by Fable. **Machine:** TCN fat node (`fatnode-tacn`).
**Worktree and branch:** `~/lsquant-unified-0007`, branch `handoff/0007-codex`.
**Build:** `~/lsquant-unified-build7c` (fresh).

**Head tested:** developer head **`5e8186d`** (`origin/handoff/0007-agi`):
- `7687487` deleted the legacy tree;
- review-2 items 1 to 8 are in `2055207..12255e9`;
- the report is at `5e8186d`.

It was merged into `handoff/0007-codex` as `9cdbe0f`. The merge touched no harness or figure file.

The first tester run (heads 6ec46f1 and 9f8ccc8, defects D1 and D2, smells S1 to S8) is in the
git history of this file (commit 65f8510).

## 0. Summary

### Byte verdict on the fat node, T = 16, every run file under `runs/goldens/` (40)

All outputs are kept in `~/s7keepc/T16`. Every golden is byte-identical except tutorial 06:

| golden | identical | differ |
|---|---:|---:|
| G1 | 1 | 0 |
| G2 | 1 | 0 |
| G3 | 3 | 0 |
| G4 | 6 | 0 |
| G4c | 2 | 0 |
| G5 | 6 | 0 |
| G6 | 5 | 0 |
| G7 | 6 | 0 |
| tutorial 01 | 24 | 0 |
| tutorial 02 | 194 | 0 |
| tutorial 03 | 9 | 0 |
| tutorial 04 | 4 | 0 |
| tutorial 05 | 17 | 0 |
| **tutorial 06** | **139** | **142** |
| MCP g1 and g5 | 7 | 0 |

- **G9.** All 12 run files run (D2 is fixed).
  - 96/96 files byte-identical to `stage5_bases` run with the same arguments (including `--csr-dir`).
  - The run-file `csr_dir` CSR files are 4/4 identical to a fresh wannier2sparse run.
  - `figures/stage_5_table.py` on the run-file outputs: 16 PASS, 0 FAIL.
- **G10.** All 11 run files run (D1 is fixed), including the new chebyshev and ischo files.
  - 44/44 files byte-identical to `stage6_nonorth` run with the run file's resolved metric.
  - The 11 argument checks pass.
  - 22/22 moment files are within the stage-6 tolerance of G1, G2 and G5 (max 9.5e-12, Löwdin/Chebyshev).
- **Tutorial 06: 142 files differ, all localised.**
  - Every numerator curve (`TimeProjSZ-*`, `TimeProjSX-*`) is identical.
  - All 139 differing curves are the denominator `TimeProj1-P*` (operator `1`). The other 3 are
    `ey_ratelaw.dat` twice, which agrees to 1.8e-16 relative, and `captured/index.tsv`.
  - **Cause:** the legacy pipeline that made the golden loaded `operators/<label>.1.CSR`, which the
    tutorial generators never write. The legacy `BuildOPFromCSRFile` returns silently on a missing
    file, so the golden denominators come from an unbuilt operator.
  - The run file uses the identity, which is correct. This needs Fable's decision (section 6).
- **MCP:** run/status of g1 and g5 gives 7/7 identical. The confinement refusal now carries its reason.
- **Refusals through the CLI:** 11/11. That includes `G+`, `G-` and `dG+/dE`, `window: "eigen"` on a
  csr model, and an unmatched reconstruction formula.
- **Greps:** all empty on code lines.
- **Section 2.5 (deleted tree):** 18/18 PASS, and the `legacy-final` build reproduces G1, G2 and G5
  byte for byte at 16, 1 and 64 threads.
- **Self-check** (`tests/dev/stage_7_check.py`): PASS once a `build/lsquant_core` link exists. It
  hard-codes `REPO/build`; see S4.

### T = 1 versus T = 16
Moved to the Aron node at Fable's direction (section 7). It is still running at hand-off.

## 1. Build and self-check

```
$ cmake -S . -B ~/lsquant-unified-build7c -DLSQUANT_BUILD_NEW_CORE=ON -DLSQUANT_BUILD_PYTHON=ON -Dpybind11_DIR=$(python3 -c 'import pybind11; print(pybind11.get_cmake_dir())')
cmake exit 0
$ cmake --build ~/lsquant-unified-build7c -j 32      (43 s: the legacy targets are gone)
build exit 0     [100%] Built target _core
```

The self-check was run first as
`OMP_NUM_THREADS=16 KPM_SEED=12345 PYTHONPATH=<build>/lsquant_core:python python3 tests/dev/stage_7_check.py`.
It exited 1 with `FileNotFoundError: .../lsquant-unified-0007/build/lsquant_core/stage6_nonorth`,
because it looks for `stage6_nonorth` under `REPO/build` (lines 370 to 372).

With the link `build/lsquant_core -> ~/lsquant-unified-build7c/lsquant_core` (inside the ignored
`build/`), it prints `STAGE 7 SELF-CHECK PASS`, exit 0. Its parts:
- (a) the catalogue round trip;
- (b) 12 refusals and checks, among them the new G+ word, eigen window on csr, unmatched formula,
  unknown generator and unknown generator parameter;
- (c) G1, G5, `g9_graphene_eigen` and `g10_chain1d_real_generalised` byte-identical;
- (d) MCP.

## 2. Harness changes (`scripts/check_goldens.sh 7`)

These are on top of the first run's harness:

- **Tutorials 01 to 06.** Every non-PNG file of `goldens/stage0/tutorials/0N` is compared at the same
  relative path; `captured/index.tsv` is checked by sha1.
- **G9.** Each run file's `model.csr_dir` (resolved against the repository) is passed to
  `stage5_bases` as `--csr-dir`. Its CSR files are cmp'd once against a fresh
  `wannier2sparse <model> 48 48 1 VX`.
- **G10.** The loop runs over every `runs/goldens/g10_*.json`:
  - `stage6_nonorth` gets the run's resolved metric (route, applier, order, coeff_tol, tol, max_iter,
    max_dim from `run.resolved.json`), and a check confirms the arguments match.
  - The reference runs once, at the first thread count, into `S7_REF6`.
  - Each file is compared with G1/G2/G5 at the stage-6 tolerance: exact 1e-10, ischo 1e-9,
    chebyshev 1e-8.
  - **Bug fixed after this run:** the relative difference is now `max_k |new - gold| / max_k |gold|`,
    the stage-6 tester's definition. The pointwise relative difference from `golden_diff.py`, which
    this run used, is about 1 on moments that vanish by symmetry and so reported 22 false FAILs. I
    re-evaluated those 22 rows from the kept outputs with the corrected definition (section 3).
- **Refusals** now also cover `g_plus`, `g_minus`, `dg_plus`, `window_eigen_csr` and
  `formula_unmatched`.
- **`stage7_deleted` (section 2.5)**, below.
- **Options:**
  - `S7_T1_SKIP`: run files that run at the first thread count only;
  - `S7_SET=list` with `S7_LIST`;
  - stages `7deleted` and `7tail` (MCP, refusals, greps, section 2.5).

How this run was executed:
- The 11 `stage6_nonorth` references were pre-computed in parallel into `~/s7ref6c` with
  `S7_REF6_WAIT_LOG`. That took 16:05 to 17:27 UTC; the two dense exact Löwdin chain1d files took
  about 80 min each under the node's load.
- The gate ran at T = 16 over all 40 run files (16:11 to 20:22 UTC). I then stopped it, killing only my
  own process tree from its recorded PID, as Fable asked. The T = 16 outputs were saved to
  `~/s7keepc/T16` first.
- MCP, refusals, greps and section 2.5 then ran as `scripts/check_goldens.sh 7tail`.

## 3. Per-golden table (fat node, T = 16, head 5e8186d)

| golden | run file(s) | identical | differ | notes |
|---|---|---:|---:|---|
| G1 | g1 | 1 | 0 | |
| G2 | g2 | 1 | 0 | |
| G3 | g3 | 3 | 0 | |
| G4 | g4 | 6 | 0 | |
| G4c | g4c | 2 | 0 | |
| G5 | g5 | 6 | 0 | |
| G6 | g6 | 5 | 0 | |
| G7 | g7, g7_sx, g7_sy, g7_sz | 6 | 0 | |
| Tut 01 | tutorial01 | 24 | 0 | index.tsv: 22/22 sha1 |
| Tut 02 | tutorial02 | 194 | 0 | index.tsv: 192/192 sha1 |
| Tut 03 | tutorial03 | 9 | 0 | |
| Tut 04 | tutorial04 | 4 | 0 | |
| Tut 05 | tutorial05 | 17 | 0 | |
| **Tut 06** | tutorial06 | **139** | **142** | 139 `TimeProj1-P*` denominator curves; `ey_ratelaw.dat` ×2 (rel 1.8e-16); index.tsv (128 of 255 entries, the denominators); section 6 |
| G9 vs stage5_bases | 12 g9 files | 96 | 0 | `--csr-dir` = the run file's csr_dir |
| G9 csr_dir vs fresh wannier2sparse | | 4 | 0 | graphene and haldane HAM, VX |
| G9 stage-5 table | | 16 PASS | 0 | real vs bloch moments up to 2.4e-14; bloch(exact) vs eigen up to 3.3e-15; eigen vs G9 up to 4.7e-15 |
| G10 vs stage6_nonorth | 11 g10 files | 44 | 0 | plus 11/11 argument checks |
| G10 vs G1/G2/G5 (stage-6 rel) | 11 g10 files | 22 PASS | 0 | see below |
| MCP | g1, g5 | 7 | 0 | |

Wall times at T = 16, from `~/s7logs/gate7c.log`; the node load was 20 to 95 throughout:

| run file | wall (s) |
|---|---:|
| `g10_chain1d_complex_lowdin` | 2506 |
| `g10_chain1d_real_lowdin` | 1049 |
| `g10_chain1d_real_lowdin_chebyshev` | 1387 |
| `g10_chain1d_real_generalised_chebyshev` | 788 |
| `tutorial06` | about 1 h 55 min |
| every other run file | under 100 |

G10 against the orthogonal goldens, with `rel = max_k |mu - mu_gold| / max_k |mu_gold|`:

| run file | rel vs G1 | rel vs G2 | tolerance |
|---|---:|---:|---|
| `g10_chain1d_real_generalised` | 4.7e-16 | 2.4e-15 | 1e-10 |
| `g10_chain1d_complex_generalised` | 1.0e-15 | 2.6e-15 | 1e-10 |
| `g10_chain1d_real_lowdin` | 8.0e-13 | 1.5e-14 | 1e-10 |
| `g10_chain1d_complex_lowdin` | 6.2e-13 | 1.2e-14 | 1e-10 |
| `g10_chain1d_real_generalised_ischo` | 2.7e-15 | 2.4e-15 | 1e-9 |
| `g10_chain1d_real_generalised_chebyshev` | 1.5e-12 | 1.0e-13 | 1e-8 |
| `g10_chain1d_real_lowdin_chebyshev` | 9.5e-12 | 1.7e-13 | 1e-8 |

| run file | rel vs G5 VX-VY | rel vs G5 VX-VX | tolerance |
|---|---:|---:|---|
| `g10_haldane_real_generalised` | 1.3e-14 | 1.1e-14 | 1e-10 |
| `g10_haldane_complex_generalised` | 1.3e-14 | 1.1e-14 | 1e-10 |
| `g10_haldane_real_lowdin` | 6.0e-14 | 8.3e-14 | 1e-10 |
| `g10_haldane_complex_lowdin` | 4.0e-14 | 5.2e-14 | 1e-10 |

All 22 rows PASS.

## 4. MCP (`python -m lsquant.mcp_server --root $TMP/s7mcp`, mcp 2.3.0)

```
MCP <- list_tools: ["words_describe", "run", "status", "reconstruct", "inspect"]
MCP -> call_tool: {"tool": "run", "arguments": {"runfile_json": "<g1.json, output.dir out/g1>"}}
MCP <- run (isError=False, returned in 0.110 s): {"run_id": "56e9b79b307b", "log": "$TMP/s7mcp/out/g1/run_56e9b79b307b.log", "output_dir": "$TMP/s7mcp/out/g1"}
MCP <- status (isError=False, returned in 0.007 s): {"run_id": "56e9b79b307b", "status": "running", "returncode": null, ...}
MCP RESULT g1: succeeded
MCP -> call_tool: {"tool": "run", "arguments": {"runfile_json": "<g5.json, output.dir out/g5>"}}
MCP <- run (isError=False, returned in 0.009 s): {"run_id": "00f304e902e5", ...}
MCP <- status (isError=False, returned in 0.005 s): {"run_id": "00f304e902e5", "status": "running", ...}
MCP RESULT g5: succeeded
MCP -> call_tool: {"tool": "run", "arguments": {"runfile_json": "<g1.json with model.dir /etc>"}}
MCP <- run (isError=True, returned in 0.028 s): "Error executing tool run: Access denied: path '/etc' resolves to '/etc' outside configured root '$TMP/s7mcp'"
MCP RESULT confinement (model.dir=/etc): refused (isError=True); out dir created: False
IDENTICAL [MCP] g1 -> G1 (1 file); g5 -> G5 (6 files)
```
The full JSON is in `~/s7logs/gate7tail.log`. The refusal reason now reaches the client, so S7 of the
first run is fixed.

## 5. Refusals (CLI, `python -m lsquant run <bad>.json --out <dir>`: exit 1, stderr, output dir not created)

| case | stderr |
|---|---|
| unknown_operator | `Unknown operator 'FOO' in word 'Tr[FOO delta(E-H)]' (registered operators: [...])` |
| rank_3 | `Unsupported word ... of rank 3 ...; nearest recognised shape is 'Tr[VX delta(E-H) VX delta(E-H)]' (kubo_greenwood)` |
| mixed_basis | `Mixed basis without explicit transform: operator 'VX@bloch' is in 'bloch' basis while calculation/active basis is 'real' ...` |
| ut_without_time | `Word 'Tr[VX delta(E-H) U(t)^+ VX U(t)]' contains time evolution U(t) without 'time' parameters ...` |
| unknown_key_top (model.dir missing) | `Unknown key 'unknown_extra_key' in run file (top-level); ...` |
| unknown_key_nested (model.dir missing) | `Unknown key 'windw' in run file (expansion); ...` |
| g_plus | `Unrecognised word shape 'Tr[VX delta(E-H) VY G+(E)]'; nearest recognised shape is 'Tr[VX delta(E-H) VY dG-(E)/dE]' (kubo_bastin)` |
| g_minus | same, for `G-(E)` |
| dg_plus | same, for `dG+(E)/dE` |
| window_eigen_csr | `expansion.window='eigen' is not supported for model.kind='csr' (no k blocks; use gershgorin or an explicit window)` |
| formula_unmatched (g1 with `bastin`) | `Reconstruction formula 'bastin' in reconstruction.formulas does not match any word in ['Tr[delta(E-H)]'] (allowed formulas for given word(s): ['dos', 'dos_curve'])` |

`words describe` refuses the same words. It still accepts `Tr[VX delta(E-H) U(t)^+ VX U(t)]`
without time (default time, by design). S1 of the first run is fixed.

## 6. Tutorial 06: localisation of the 142 differing files

**What differs.** The numerators are fine and only the denominators fail:
- Every numerator curve is byte-identical: `meanTimeProjSZ-PZ*`, `meanTimeProjSX-PX*`, both
  top-level and captured.
- Every denominator curve differs: `meanTimeProj1-P*`, the word `Tr[P U(t)^+ 1 U(t) delta(E-H) P]`.
  First differing line of `captured/0001_meanTimeProj1-PZey_L20_D0.05...dat`: new `0 0.0532345`,
  golden `0 9.98496e-07`.
- The new curve is constant in t and equal to the numerator at t = 0, which is
  `Tr[P delta P]`, the correct value for the identity.
- The golden values are about 1e-6 and noisy.

**Origin.** It is a legacy defect in the golden, not a new defect:
- The stage-0 pipeline (`examples/06_spin_relaxation_square/lsqrelax.py:57`, `rp_spin.py:51,64`) calls
  `inline_compute-kpm-TimeEvProjetedOp <label> 1 <proj> ...`.
- That driver builds its operators from `operators/<label>.<ID>.CSR`
  (`legacy-final:src/inline_compute-kpm-TimeEvProjetedOp.cpp:53-62`).
- `make_ey.py` and `make_rashba.py` never write a `<label>.1.CSR`.
- `legacy-final:include/sparse_matrix.hpp:132-133`: `BuildOPFromCSRFile` returns silently when
  `OPERATOR_FromCSRFile` fails, so `OP[1]` stays an unbuilt matrix.

**Why the stage-0 figure was still right.** The tutorial normalises the curve as
`s = num / mean(den)`, then `s / s[0]`. The denominator's scale cancels:
- `ey_ratelaw.dat` agrees to 1.8e-16 relative (rounding);
- the stage-0 figures matched.

**Decision for Fable.** Two options:
- re-bless the 139 denominator curves, the two `ey_ratelaw.dat` and `captured/index.tsv` from the new
  code, with this reason, as was done for G4c / D2;
- or have the run file emulate the defect (not recommended).

Nothing under `goldens/` was touched.

## 7. T = 1 versus T = 16 on the Aron node (in progress at hand-off)

Following Fable, the thread-count identity is checked on `tcn-aronnode` (96 cores), not on the fat node.

| item | value |
|---|---|
| worktree | `~/lsquant-codex-7` at `5e8186d`, detached |
| build | `~/lsquant-codex-7-build`, gcc 9.4.0, Python 3.10.22 |
| cmake | the `~/lsquant/aron_build.sh` line: `lsq-env` Eigen, FFTW and GSL, Release |
| `build/w2s_csr` | copied from the fat-node worktree; Aron has no wannier2sparse binary (S3) |
| run script | `~/s7logs/aron_runs7.sh`, root PID 433093 |
| outputs | `~/s7aron/T16` and `~/s7aron/T1` |
| log | `~/s7aron/runs.log` |

The script runs all 40 run files at T = 16 one after another, and at T = 1 up to 24 at a time,
longest first.

Status at hand-off (22:29 CEST):
- 37 of 40 T = 1 runs are done, all with exit 0. Still running: `tutorial06` and the two dense
  Löwdin chain1d files.
- The T = 16 queue is on `tutorial06`, its first file.

To finish: `bash ~/s7logs/aron_cmp7.sh` on Aron. It cmp's every file of `T1/<run file>` against
`T16/<run file>` (excluding `run.resolved.json`) and prints IDENTICAL/DIFFERS per run file, with
totals.

Aron is not byte-identical to the fat-node goldens: gcc 9.4, Eigen and glibc 2.31 give last-bit
differences, measured in `~/lsquant/aron_build.log` and `/tmp/aron_g5diff` on Aron. Only T = 1 versus
T = 16 is compared there.

## 8. Section 2.5: the deleted legacy tree (`scripts/check_goldens.sh 7tail`, HEAD 9cdbe0f)

```
PASS  src/, include/, test/ absent from the working tree
PASS  git ls-files src include test: 0 files
PASS  build ~/lsquant-unified-build7c configured from this tree
PASS  CMakeLists.txt: 0 legacy target lines
PASS  scripts/regenerate_goldens.sh removed
PASS  legacy shell scripts under examples/: 0
PASS  tag legacy-final exists (commit 05e2f655); it holds the legacy tree (187 files under src include test); it is an ancestor of HEAD
PASS  git worktree add /tmp/legacy-opus7 legacy-final
PASS  legacy-final: cmake -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON && cmake --build (exit 0)
PASS  legacy binaries inline_compute-kpm-{nonEqOp,spectralOp}, inline_kubo{Bastin,BastinSea,BastinSurf,Greenwood}FromChebmom built
stage 0 with the legacy-final build: IDENTICAL G1, G2 (T = 16, 1, 64), G5 moments and 4 curves (T = 16): 12/12
INFO  worktree /tmp/legacy-opus7 removed
```
The provenance of the goldens stays reproducible from `legacy-final`.

## 9. Behaviour of `ensure_w2s_csr_dir` (review-2 item 6), as Fable asked

The tests are in `~/s7logs/w2stest.log`. Each runs `g9_graphene_real.json` with a different `csr_dir`.

| case | setup | result |
|---|---|---|
| A | committed run file (`csr_dir: build/w2s_csr/graphene`); the directory did not exist; `~/wannier2sparse/build/wannier2sparse` present | Exit 0. The run silently created `<worktree>/build/w2s_csr/graphene/{graphene.HAM.CSR, graphene.VX.CSR, graphene.out}` inside the repository's ignored `build/` (`git status --ignored` shows `!! build/`). Those CSRs are byte-identical to a fresh wannier2sparse run, and the run's outputs are byte-identical to `stage5_bases --csr-dir`. |
| B | csr_dir absent, no binary reachable (`HOME` redirected, `W2S_BIN` unset, not in PATH: `find_w2s_binary() -> None`) | Exit 1 with `RealBasis: '.../csrB/graphene.HAM.CSR' not found in csr_dir '.../csrB'`. The generation step does nothing silently and the load then fails. The message does not say that wannier2sparse was not found or where it was looked for. `run.resolved.json` was already written. |
| C | csr_dir exists, non-empty, holds a 6×4×1 wannier2sparse output | Not regenerated. Exit 1 with `receipt '.../graphene.out' supercell (6, 4, 1) does not match requested supercell (48, 48, 1)`: the receipt check catches it. |
| D | csr_dir exists but is empty | Generated. Exit 0, and moments and curves are identical to A. Only the `source` paths in `*.basis.json` differ. |

## 10. Figure

`figures/out/stage_7/fig_stage_7.png`, from `python3 figures/stage_7.py --runs ~/s7keepc/T16`.
Caption (condensed from `fig_stage_7.caption.txt`):

> FIG. 7. The six tutorial figures of stage 0 (left, legacy build, figures/out/stage_0) beside the
> same curves drawn from the outputs of `python -m lsquant run runs/goldens/tutorial0N.json` (middle,
> solid) over the stage-0 captures (dashed), with the residual new - golden (right). OMP_NUM_THREADS =
> 16, KPM_SEED = 12345, head 5e8186d.
> (a) to (e) tutorials 01 to 05: 24/24, 194/194, 9/9, 4/4 and 17/17 golden files byte-identical
> (every non-PNG file; captured/index.tsv by sha1); max |new - golden| 0.
> (f) tutorial 06 from runs/goldens/tutorial06.json: 139 of 281 byte-identical. The 141 differing data
> files are the 139 denominator curves of the word with operator '1' (TimeProj1-P*),
> captured/0254_ey_ratelaw.dat and ey_ratelaw.dat; captured/index.tsv also differs on those entries.
> Drawn: the EY numerator curves, all byte-identical (residual 0). See tests.md section 6 for the
> legacy defect behind the golden denominators.

## 11. Smells

- **S1 (fixed).** `G+`, `G-` and `dG+/dE` are now refused.
- **S2 (mostly fixed).** Tutorials 01, 02 and 05 now compile the run file's word with its
  parameters. Residual: `apps/tutorial06.py:125-127` builds its own word strings
  (`Tr[P:{proj} U(t)^+ {op} U(t) delta(E-H) P:{proj}]`). The run file's `words` are only validated
  (`tutorial06.py:231-232`).
- **S3 (new, item 6).**
  - `ensure_w2s_csr_dir` writes into the repository tree (`build/w2s_csr/...`) as a side effect of a run.
  - It depends on a machine-local binary path (`~/wannier2sparse/build/wannier2sparse`, `W2S_BIN`, PATH).
    On Aron the binary is absent and the run files only work after the CSRs are copied.
  - Without the binary the error does not name the missing binary (case B).
  - A non-empty `csr_dir` is never regenerated. Only the receipt's supercell check protects against
    stale files (case C).
- **S4 (new).** `tests/dev/stage_7_check.py:370-372` hard-codes `REPO/build[/lsquant_core]/stage6_nonorth`,
  so it fails with any other build directory.
- **S5 (legacy left behind, not in the brief's list).**
  - `legacy/` still tracks 19 legacy driver sources (quarantine of 5b5db0c).
  - Documentation and example scripts still name `inline_*` drivers (AGENTS.md, README.md,
    `agent_template/`, `docs/spectral_scales.md`, `examples/*/lsq*.py`, ...).
  - `examples/06/lsqrelax.py` and `rp_spin.py` still call the deleted drivers.
- **S6 (golden defect).** Tutorial 06 denominators come from a legacy load of a missing operator file
  (section 6). This is a candidate for re-blessing, not a code defect.
- **S7 (cost).**
  - The dense exact Löwdin chain1d run files take 17 to 42 min at 16 threads on a loaded node.
  - `tutorial06` takes about 2 h.
  - Neither is feasible at 1 thread on the shared fat node, hence `S7_T1_SKIP` and the Aron split.
- **Not found:**
  - no `eval`, `exec` or `literal_eval` in the word language, run file, CLI or MCP;
  - no basis or word names under core (greps empty);
  - no golden re-blessed or touched.

## 12. Run history

| time (UTC) | event |
|---|---|
| 16:04 | merge of 5e8186d; fresh build |
| 16:05 to 17:27 | stage6 references pre-computed (`~/s7ref6c`, my PIDs, reniced to 19) |
| 16:11 to 20:22 | gate at T = 16 over all 40 run files; stopped at the start of the T = 1 phase by `~/s7logs/watch.sh`, which killed only the process tree of my recorded PID after saving the T = 16 outputs |
| 20:22 to 20:40 | `7tail` |
| 19:50 onwards | T = 1 versus T = 16 on Aron |
