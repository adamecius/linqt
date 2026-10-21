# `lsquant` — Python Library, CLI, and MCP Server

`lsquant` compiles observable words `Tr[A_0 f_1(H) A_1 ... f_p(H) A_p]` into bytecode `Plan`s for the C++ Kernel Polynomial Method engine (`lsquant._core`) and reconstructs spectral and transport curves.

## Command-Line Interface (`python -m lsquant` / `lsquant`)

```bash
# List recognised word shapes in the catalogue
python -m lsquant words

# Describe a specific observable word (AST, bytecode plan, operators, reconstruction spec)
python -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"

# Validate and execute a run file
python -m lsquant run runs/goldens/g5.json --out /tmp/haldane_out

# Reconstruct a curve from a moment file or output directory
python -m lsquant reconstruct /tmp/haldane_out --formula bastin --broadening-meV 10

# Inspect a moment file, run file, or output directory
python -m lsquant inspect /tmp/haldane_out
```

### Generator Catalogue (`model.kind = "generator"`)

`model.name` is validated against a fixed catalogue in `python/lsquant/runfile.py` and resolved relative to the repository root:

| `model.name` | Module Path |
|---|---|
| `make_chain` | `examples/01_dos_1d_chain/make_chain.py` |
| `make_disordered_chain` | `examples/02_localization_1d_chain/make_disordered_chain.py` |
| `make_haldane` | `examples/03_haldane_hall/make_haldane.py` |
| `make_chain_mag` | `examples/04_spin_precession_magnetic_chain/lsqspin.py` |
| `make_square` | `examples/05_charge_transport_square/make_square.py` |
| `gen_models` | `python/lsquant/models/wannier.py` |

## Run-File Schema

Every calculation is configured by a strictly validated JSON run file (`python/lsquant/runfile.py`). Unknown keys at the top level or inside any sub-object are refused before any file is loaded. Fields marked **[extension]** are schema additions beyond `handoff/0007-words-cli-mcp/brief.md` section 4 needed to reproduce the multi-word, multi-target, and tutorial golden fixtures.

| Section / Field | Allowed Values / Type | Default | Brief / Extension | Description |
|---|---|---|---|---|
| `model.kind` | `"csr"` \| `"wannier"` \| `"generator"` | `"csr"` | `"csr"`, `"wannier"` in brief; `"generator"` **[extension]** | Model source (`"generator"` runs a named model generator from `GENERATOR_CATALOGUE` for Tutorials 01, 02, and 05). |
| `model.dir` | path | `null` | Brief | Directory containing `<label>.<OP>.CSR` (`csr`) or `<seed>_hr.dat`/`.uc`/`.xyz` (`wannier`). |
| `model.label` | string | `null` | Brief | Operator file prefix `<label>.<OP>.CSR` when `model.kind == "csr"`. |
| `model.seed` | string | `null` | Brief | Wannier90 seed name when `model.kind == "wannier"`. |
| `model.operators` | `list[str]` | `[]` | Brief | Static operator names available in `model.dir`. |
| `model.csr_dir` | path | `null` | **[extension]** | Pre-built supercell CSR directory for `model.kind == "wannier"` (G9 real-space runs). |
| `model.name` | `"make_chain"` \| `"make_disordered_chain"` \| `"make_haldane"` \| `"make_square"` \| `"make_chain_mag"` \| `"gen_models"` | `null` | **[extension]** | Generator name from the fixed `GENERATOR_CATALOGUE` in `runfile.py`, resolved from the repository root (Tutorials 01, 02, 05). |
| `model.tutorial` | `1..5` | `null` | **[extension]** | Tutorial driver number when `model.kind == "generator"` (Tutorials 01, 02, 05). |
| `model.params` | `dict` | `{}` | **[extension]** | Generator parameter dictionary (`sizes`, `hopping_t`, `W`, `realisations`, `L`, `disorder_seed`) (Tutorials 01, 02, 05). |
| `basis.name` | `"real"` \| `"bloch"` \| `"eigen"` | `"real"` | Brief | Representation basis. |
| `basis.supercell` | `[N1, N2, N3]` | `[1, 1, 1]` | Brief | Supercell or k-grid dimensions. |
| `basis.disorder` | `{"W": float, "seed": int}` | `{"W": 0.0, "seed": 7}` | Brief | Real-space Anderson box disorder parameters. |
| `metric.policy` | `"orthogonal"` \| `"nonorthogonal"` | `"orthogonal"` | Brief | Basis overlap policy. |
| `metric.route` | `"lowdin"` \| `"generalised"` | `"lowdin"` | Brief | Non-orthogonal evaluation route. |
| `metric.applier` | `"exact"` \| `"chebyshev"` \| `"ischo"` | `"exact"` | Brief (`"ischo"` from Stage 6) | Overlap inverse / inverse-square-root applier. |
| `metric.order` | `int` | `60` | Brief | Chebyshev polynomial order for overlap appliers. |
| `metric.coeff_tol` | `float` | `1e-12` | **[extension]** | Chebyshev coefficient truncation threshold for overlap appliers (G10). |
| `metric.tol` | `float` | `1e-13` | Brief | PCG solver tolerance for `"ischo"`. |
| `metric.max_iter` | `int` | `200` | **[extension]** | Maximum PCG iterations for `"ischo"` (G10). |
| `metric.max_dim` | `int` | `4096` | **[extension]** | Maximum dense matrix dimension for `"exact"` overlap appliers (G10). |
| `word` | string | *(required)* | Brief | Primary observable word `Tr[...]`. |
| `words` | `list[str]` | `[word]` | **[extension]** | Multiple observable words evaluated in a single run (G7, G9, G10 chain1d, Tutorials 04 and 05). |
| `expansion.moments` | `int` | `64` | Brief | Number of Chebyshev moments `M` per expanded slot. |
| `expansion.window` | `"gershgorin"` \| `"eigen"` \| `[lo, hi]` | `"gershgorin"` (`real`), `"eigen"` (`bloch`, `eigen`) | Brief | Spectral window mode or explicit `[lo, hi]` bounds in eV (default is a function of `basis.name`). |
| `expansion.trace` | `"exact"` \| `"random"` | `"exact"` | Brief | Trace evaluation mode. |
| `expansion.states` | `int` \| `null` | `null` | Brief | Number of trace states (`dim` for `exact`, `1` for `random` when omitted). |
| `expansion.seed` | `int` | `12345` | Brief | Random-phase state generator seed (`KPM_SEED`). |
| `expansion.batch` | `int` \| `null` | `moments` | Brief | Engine vector batch size. |
| `expansion.state_label` | string \| `null` | `null` | **[extension]** | Explicit `_state<label>` suffix override when needed by a golden filename. |
| `time.steps` | `int` | — | Brief | Number of time steps (`> 1`). |
| `time.tmax` | `float` \| `null` | `20.0` | Brief | Maximum evolution time in fs. |
| `time.dt` | `float` \| `null` | `null` | **[extension]** | Explicit time step in fs as an alternative to `tmax = dt * (steps - 1)`. |
| `time.bessel` | `"scipy"` \| `"legacy"` | `"scipy"` | Brief (Section 5) | Bessel coefficient backend (`"legacy"` uses GSL values for stage-0 byte identity). |
| `reconstruction.kernel` | `"jackson"` \| `"lorentz"` \| `"none"` | `"jackson"` | `"jackson"` in brief; `"lorentz"`, `"none"` **[extension]** | Damping kernel (`lorentz` and `none` expose the additional kernels in `python/lsquant/kernels.py`). |
| `reconstruction.broadening_meV` | `float` \| `null` | `null` | Brief | Broadening in meV. |
| `reconstruction.fermi_energies` | `list[float]` | `[]` | Brief | Fermi energies in eV for time-correlation reconstruction. |
| `reconstruction.formulas` | `list[str]` \| `null` | shape default | Brief | Reconstruction formulas (`dos`, `greenwood`, `bastin`, `sea`, `surf`, `greenwood_fft`, `bastin_fft`, `time_correlation`, `row_m0`). |
| `reconstruction.num_points` | `int` \| `null` | `null` | **[extension]** | Explicit energy grid point count override. |
| `reconstruction.fft_grid` | `bool` | `false` | **[extension]** | Reconstruct on the Chebyshev-node FFT grid (`Kubo*_JACKSON_FFTgrid.dat` in G4). |
| `reconstruction.legacy_text_round_trip` | `bool` | `false` | Brief (Section 5) | Round-trip moments through the 15-digit legacy ASCII file before reconstruction (Tutorial 03, G4, G5, G6, G10). |
| `reconstruction.per_site_dos` | `bool` | `false` | **[extension]** | Also write standalone per-site `DOS_reconstructed.dat` (G3). |
| `reconstruction.targets` | `list[dict]` | `[]` | **[extension]** | Additional `(formulas, broadening_meV, subdir, also_top_dir)` reconstruction targets from the same moments (G4 `64` meV top-level + `10` meV `sea_surf/` and `fft_grid/`). |
| `output.dir` | path | `"."` | Brief | Output directory (always receives `run.resolved.json`). |
| `output.legacy_names` | `bool` | `true` | Brief | Use legacy stage-0 file naming conventions. |
| `output.legacy_layout` | `bool` | `true` | Brief | Write legacy `.chebmom1D`/`.chebmom2D`/`.chebmomTD` text files. |
| `output.write_chebmom` | `bool` | `false` | **[extension]** | Also write the unified `.chebmom` file alongside legacy moment files (G9 and G10). |
| `output.write_basis_json` | `bool` | `false` | **[extension]** | Write `<stem>.basis.json` provenance metadata (G9). |
| `output.capture_tutorial` | `bool` | `false` | **[extension]** | Copy outputs into `captured/NNNN_*` for tutorial golden verification (Tutorials 01–05). |

## MCP Server (`python -m lsquant.mcp_server`)

`lsquant.mcp_server` exposes five tools over `stdio` using `FastMCP` (`pip install mcp`):

- `words_describe(word)` — parses and compiles `word`, returning the same JSON structure as `python -m lsquant words describe "<word>"`.
- `run(runfile_json)` — validates `runfile_json` and dispatches `python -m lsquant run` in a background subprocess, returning `{"run_id", "log", "output_dir"}` immediately.
- `status(run_id)` — polls a background calculation started by `run`.
- `reconstruct(output_dir, formula, broadening_meV, fermi_energies)` — reconstructs curves from moment files in `output_dir`.
- `inspect(path)` — inspects a moment file, run file, or output directory, returning the same JSON structure as `python -m lsquant inspect <path>`.

All file operations are restricted to `LSQUANT_MCP_ROOT` (defaults to the repository root or `--root <dir>`).

### Claude Code `mcp.json` Configuration

Add the following entry to `.mcp.json` (or `~/.claude/mcp.json`):

```json
{
  "mcpServers": {
    "lsquant": {
      "command": "python3",
      "args": ["-m", "lsquant.mcp_server"],
      "env": {
        "PYTHONPATH": "/path/to/lsquant/build/lsquant_core:/path/to/lsquant/python",
        "LSQUANT_MCP_ROOT": "/path/to/lsquant",
        "OMP_NUM_THREADS": "16",
        "KPM_SEED": "12345"
      }
    }
  }
}
```
