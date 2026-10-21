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
