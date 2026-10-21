# Stage 7 Report — Observable Word Language, Run File, CLI, and MCP Server

Branch: `handoff/0007-agi` (merged `origin/handoff/0006-agi` and `origin/cpp-unified`).

---

## 1. Files Created and Modified

| File | Purpose |
|------|---------|
| `python/lsquant/words/grammar.py` | Hand-written tokenizer (`tokenize`) and recursive-descent parser (`parse`) producing `Word` AST (`OperatorSlot`, `SpectralSlot`, `CommutatorSlot`) with basis-transform sandwich folding (`B^+ O B`, `U^+ O U`). |
| `python/lsquant/words/catalogue.py` | Catalogue table (`CATALOGUE`) of the 7 recognised observable word shapes, basis-consistency validator (`validate_word_bases`), and nearest-shape diagnostic matcher (`match_catalogue`, `nearest_catalogue_entry`). |
| `python/lsquant/words/compiler.py` | Word compiler (`compile`, `describe`, `list_catalogue`) mapping `(Word, basis, params)` to `(Plan, operators, ReconstructionSpec)`, including orthogonal and non-orthogonal (`lowdin`, `generalised`) operator wrapping. |
| `python/lsquant/words/__init__.py` | Public exports of `lsquant.words` (`parse`, `compile`, `describe`, `CATALOGUE`, AST classes). |
| `python/lsquant/moments.py` | Extended `dos_plan` with optional `operator_index` keyword argument for `Tr[A delta(E-H)]` (`A != 1`). |
| `python/lsquant/runfile.py` | `RunFile` dataclass, strict schema validator (`validate_runfile_dict`) refusing unknown keys before any file load, defaults table (`DEFAULT_RUNFILE_DEFAULTS`), execution engine (`execute_runfile`), curve reconstructor (`reconstruct_target`), and inspector (`inspect_target`). |
| `python/lsquant/cli.py` | `argparse` CLI (`run`, `reconstruct`, `inspect`, `words`). |
| `python/lsquant/__main__.py` | Entry point for `python -m lsquant`. |
| `python/lsquant/mcp_server.py` | Stdio MCP server (`words_describe`, `run`, `status`, `reconstruct`, `inspect`) sharing validation and execution logic with the CLI and enforcing root-directory confinement. |
| `python/lsquant/README.md` | CLI usage and Claude Code `mcp.json` configuration entry. |
| `python/pyproject.toml` | Added `[project.scripts]` console entry `lsquant = "lsquant.cli:main"`. |
| `python/lsquant/__init__.py` | Re-exported `words`, `runfile`, and `cli`. |
| `python/lsquant/apps/tutorial01.py`, `python/lsquant/apps/tutorial03.py` | Added top-level 120 meV curve and `captured/NNNN_*` outputs for Tutorial 01 and Tutorial 03. |
| `runs/goldens/*.json` | 36 run files covering G1, G2, G3, G4, G4c, G5, G6, G7, G9 (both models and all bases/disorders), G10 (all 4 fixtures on both routes), and Tutorials 01–05. |
| `tests/dev/stage_7_check.py` | Single developer self-check covering (a) catalogue canonical round-trip, (b) 5 malformed/unsupported refusals, (c) CLI `g1.json` and `g5.json` byte identity vs G1 and G5, and (d) MCP server subprocess equivalence vs CLI. |

Zero C++ files under `lsquant/core/`, `lsquant/io/`, or `lsquant/api/` were modified. The legacy tree (`src/`, `include/`, `test/`) was untouched.

---

## 2. Word Language Grammar (EBNF)

```ebnf
word          ::= "Tr[" body "]" | "sea[" "Tr[" product "]" "]" ;
body          ::= "sea[" product "]" | product ;
product       ::= factor ( ["*"] factor )* ;
factor        ::= spectral | commutator | transform_seq | operator ;
spectral      ::= "delta(E-H)"
                | "G+(E)" | "G-(E)"
                | "dG-(E)/dE" | "dG+(E)/dE"
                | "U(t)^+" | "U(t)" ;
commutator    ::= "[" axis "," "U(t)" "]^+"
                | "[" axis "," "U(t)" "]" ;
axis          ::= "X" | "Y" | "Z" | "VX" | "VY" | "VZ" ;
transform_seq ::= "B^+" operator "B"
                | "B" operator "B^+"
                | "U^+" operator "U"
                | "U" operator "U^+" ;
operator      ::= ("1" | "P:" ident | ident) ["^+"] ["@" basis_tag] ;
ident         ::= LETTER (LETTER | DIGIT | "_")* ;
basis_tag     ::= "real" | "bloch" | "eigen" ;
```

---

## 3. Catalogue of Recognised Word Shapes

| `shape_id` | Canonical Word | Plan Builder | Default Reconstruction | Legacy Moment Prefix | Legacy Curve Prefix | Rank | Time-Dependent |
|---|---|---|---|---|---|---:|---|
| `dos_identity` | `Tr[delta(E-H)]` | `dos_plan` | `dos_curve` | `SpectralOp1` | `mean` | 1 | No |
| `dos_operator` | `Tr[A delta(E-H)]` | `dos_plan` (`operator_index=0`) | `dos_curve` | `SpectralOp<A>` | `mean` | 1 | No |
| `kubo_greenwood` | `Tr[VX delta(E-H) VX delta(E-H)]` | `kubo_plan` | `kubo_greenwood` | `NonEqOp<A>-<B>` | `KuboGreenWood_` | 2 | No |
| `kubo_bastin` | `Tr[VX delta(E-H) VY dG-(E)/dE]` (or `Tr[sea[...]]`) | `kubo_plan` | `kubo_bastin` (`sea`, `surf`) | `NonEqOp<A>-<B>` | `KuboBastin*_` | 2 | No |
| `msd` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `msd_plan` | `time_correlation_curve` | `Correlation<V>` | `mean...EF` | 1 | Yes |
| `projected_evolution` | `Tr[P U(t)^+ SX U(t) delta(E-H) P]` | `projected_evolution_plan` | `row_m0` | `TimeProj<A>-<P>` | `TimeProj` | 1 | Yes |
| `time_correlation` | `Tr[A delta(E-H) U(t)^+ B U(t)]` | `time_correlation_plan` | `time_correlation_curve` | `TimeCorr<A>-<B>` | `mean...EF` | 1 | Yes |

---

## 4. Run-File Schema and Defaults

Unknown keys at the top level or inside any sub-object raise `ValueError` during `validate_runfile_dict` before any file is loaded.

```json
{
  "model": {
    "kind": "csr | wannier | generator",
    "dir": "<path>",
    "label": "<csr_label>",
    "seed": "<wannier_seed>",
    "operators": ["VX", "VY"],
    "csr_dir": null,
    "tutorial": null,
    "script": null,
    "params": {}
  },
  "basis": {
    "name": "real | bloch | eigen",
    "supercell": [1, 1, 1],
    "disorder": {"W": 0.0, "seed": 7}
  },
  "metric": {
    "policy": "orthogonal | nonorthogonal",
    "route": "lowdin | generalised",
    "applier": "exact | chebyshev | ischo",
    "order": 60,
    "coeff_tol": 1e-12,
    "tol": 1e-13,
    "max_iter": 200,
    "max_dim": 4096
  },
  "word": "Tr[...]",
  "words": ["Tr[...]"],
  "expansion": {
    "moments": 64,
    "window": "gershgorin | eigen | [lo, hi]",
    "trace": "exact | random",
    "states": null,
    "seed": 12345,
    "batch": null,
    "state_label": null
  },
  "time": {
    "steps": 64,
    "tmax": 20.0,
    "dt": null,
    "bessel": "scipy | legacy"
  },
  "reconstruction": {
    "kernel": "jackson | lorentz | none",
    "broadening_meV": null,
    "fermi_energies": [],
    "formulas": null,
    "num_points": null,
    "fft_grid": false,
    "legacy_text_round_trip": false,
    "per_site_dos": false,
    "targets": []
  },
  "output": {
    "dir": ".",
    "legacy_names": true,
    "legacy_layout": true,
    "write_chebmom": false,
    "write_basis_json": false,
    "capture_tutorial": false
  }
}
```

Every run writes `run.resolved.json` in its output directory.

---

## 5. Run Files under `runs/goldens/` and Golden Mapping

| Run File | Golden / Tutorial Reproduced |
|---|---|
| `runs/goldens/g1.json` | **G1**: `goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D` |
| `runs/goldens/g2.json` | **G2**: `goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D` |
| `runs/goldens/g3.json` | **G3**: `goldens/stage0/chain1d/meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat`, `DOS_reconstructed.dat`, and G1 `.chebmom1D` |
| `runs/goldens/g4.json` | **G4**: `goldens/stage0/graphene/NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D`, `KuboBastin_...JACKSON.dat` (64 meV), `sea_surf/KuboBastin{Sea,Surf}_...JACKSON.dat` (10 meV), `fft_grid/Kubo{GreenWood,Bastin}_...JACKSON_FFTgrid.dat` (10 meV) |
| `runs/goldens/g4c.json` | **G4c**: `goldens/stage0/graphene/offcentre_enclosing/NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact.chebmom2D` and `KuboBastin_...JACKSON.dat` (10 meV) |
| `runs/goldens/g5.json` | **G5**: `goldens/stage0/haldane/NonEqOpVX-{VY,VX}haldaneKPM_M128x128_stateexact.chebmom2D` and `Kubo{Bastin,BastinSea,BastinSurf,GreenWood}_...JACKSON.dat` (10 meV) |
| `runs/goldens/g6.json` | **G6**: `goldens/stage0/chain1d_ballistic/CorrelationVXchain1dKPM_M256_state20exact.chebmomTD` and the 4 `mean...EF{0.000000,0.500000,1.000000,1.500000}JACKSON.dat` curves |
| `runs/goldens/g7.json` (`g7_sx.json`, `g7_sy.json`, `g7_sz.json`) | **G7**: `goldens/stage0/chain1d_mag/TimeProj{SX,SY,SZ}-PXchain1d_magKPM_M32_stateexact.chebmomTD` |
| `runs/goldens/g9_graphene_{real,bloch,bloch_exact,eigen,real_w1,bloch_w1}.json` | **G9 (graphene)**: `goldens/stage0/bases/graphene/` DOS and `sigma_xx` across `real`, `bloch`, and `eigen` bases at `W = 0` and `W = 1` eV |
| `runs/goldens/g9_haldane_{real,bloch,bloch_exact,eigen,real_w1,bloch_w1}.json` | **G9 (haldane)**: `goldens/stage0/bases/haldane/` DOS and `sigma_xx` across `real`, `bloch`, and `eigen` bases at `W = 0` and `W = 1` eV |
| `runs/goldens/g10_chain1d_{real,complex}_{lowdin,generalised}.json` | **G10 (chain1d)**: `goldens/stage0/nonorth/chain1d/{real,complex}` DOS and Kubo moments on `lowdin` and `generalised` routes vs G1/G2 |
| `runs/goldens/g10_haldane_{real,complex}_{lowdin,generalised}.json` | **G10 (haldane)**: `goldens/stage0/nonorth/haldane/{real,complex}` `VX-VY` and `VX-VX` moments and Hall curves on `lowdin` and `generalised` routes vs G5 |
| `runs/goldens/tutorial01.json` | **Tutorial 01**: `goldens/stage0/tutorials/01/` (`captured/0000_*`..`0021_*`, `eta{10,20,24,80}/`, top-level `.chebmom1D` and 120 meV `.dat`) |
| `runs/goldens/tutorial02.json` | **Tutorial 02**: `goldens/stage0/tutorials/02/` (`captured/0000_*`..`0191_*` and top-level `N = 512` `.dat`) |
| `runs/goldens/tutorial03.json` | **Tutorial 03**: `goldens/stage0/tutorials/03/` (`captured/0000_*`..`0003_*` and top-level 4 Haldane `.dat` curves) |
| `runs/goldens/tutorial04.json` | **Tutorial 04**: `goldens/stage0/tutorials/04/` (`captured/0000_*`..`0002_*` and top-level 3 `TimeProj` `.chebmomTD` files) |
| `runs/goldens/tutorial05.json` | **Tutorial 05**: `goldens/stage0/tutorials/05/` (`captured/0000_*`..`0007_*` and top-level 8 MSD and DOS `.dat` curves) |

---

## 6. Self-Check Output and MCP Transcript (`tests/dev/stage_7_check.py` on TCN)

Executed on TCN (`~/lsquant-agi-build`):

```text
$ OMP_NUM_THREADS=16 KPM_SEED=12345 PYTHONPATH=build/lsquant_core:python python3 tests/dev/stage_7_check.py
== (a) Catalogue canonical parse/print round-trip ==
  PASS  Tr[delta(E-H)]                             -> shape=dos_identity         plan_builder=dos_plan
  PASS  Tr[A delta(E-H)]                           -> shape=dos_operator         plan_builder=dos_plan
  PASS  Tr[VX delta(E-H) VX delta(E-H)]            -> shape=kubo_greenwood       plan_builder=kubo_plan
  PASS  Tr[VX delta(E-H) VY dG-(E)/dE]             -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]         -> shape=msd                  plan_builder=msd_plan
  PASS  Tr[P U(t)^+ SX U(t) delta(E-H) P]          -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[A delta(E-H) U(t)^+ B U(t)]             -> shape=time_correlation     plan_builder=time_correlation_plan
  PASS  Tr[sea[VX delta(E-H) VY dG-(E)/dE]]        -> shape=kubo_bastin          plan_builder=kubo_plan
  PASS  Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]    -> shape=projected_evolution  plan_builder=projected_evolution_plan
  PASS  Tr[B^+ VX@bloch B delta(E-H)]              -> shape=dos_operator         plan_builder=dos_plan

== (b) Refusal of 5 malformed / unsupported inputs ==
  PASS  [unknown_operator] refused with: Unknown operator 'FOO' in word 'Tr[FOO delta(E-H)]' (registered operators: ['1', 'VX', 'VY'])
  PASS  [rank_3] refused with: Unsupported word 'Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]' of rank 3 (3 expanded spectral slots): rank >= 3 moment tensors are out of scope; nearest recognised shape is 'Tr[VX delta(E-H) VX delta(E-H)]' (kubo_greenwood)
  PASS  [mixed_basis_without_transform] refused with: Mixed basis without explicit transform: operator 'VX@bloch' is in 'bloch' basis while calculation/active basis is 'real' (insert explicit B^+ ... B or U^+ ... U transform)
  PASS  [ut_without_time] refused with: Word 'Tr[A delta(E-H) U(t)^+ B U(t)]' contains time evolution U(t) without 'time' parameters (steps and tmax required)
  PASS  [unknown_key_in_runfile] refused with: Unknown key 'unknown_extra_key' in run file (top-level); allowed keys are ['basis', 'expansion', 'metric', 'model', 'output', 'reconstruction', 'time', 'word', 'words']
  PASS  [runfile_schema] all 36 files in runs/goldens/*.json validated

== (c) CLI reproduction of G1 (runs/goldens/g1.json) and G5 (runs/goldens/g5.json) ==
  PASS  cmp IDENTICAL: G1 (SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D)
  PASS  cmp IDENTICAL: G5 (KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat)
  PASS  cmp IDENTICAL: G5 (KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat)

== (d) MCP server subprocess vs CLI JSON equivalence ==
  MCP tools registered: ['words_describe', 'run', 'status', 'reconstruct', 'inspect']
  MCP -> call_tool: {"tool": "words_describe", "arguments": {"word": "Tr[VX delta(E-H) VY dG-(E)/dE]"}}
  MCP <- words_describe response:
    {
      "word": {
        "canonical": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "sea": false,
        "rank": 2,
        "time_dependent": false,
        "slots": [
          {
            "type": "operator",
            "name": "VX",
            "operator_key": "VX"
          },
          {
            "type": "spectral",
            "kind": "delta",
            "token": "delta(E-H)",
            "mode": "expand"
          },
          {
            "type": "operator",
            "name": "VY",
            "operator_key": "VY"
          },
          {
            "type": "spectral",
            "kind": "dG-/dE",
            "token": "dG-(E)/dE",
            "mode": "expand"
          }
        ]
      },
      "catalogue": {
        "shape_id": "kubo_bastin",
        "canonical_example": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "plan_builder": "kubo_plan",
        "default_reconstruction": "kubo_bastin",
        "default_formulas": [
          "bastin"
        ],
        "engine_word": "VX delta VY delta",
        "legacy_moment_prefix": "NonEqOpVX-VY",
        "legacy_moment_ext": ".chebmom2D",
        "operators": [
          {
            "type": "operator",
            "name": "VX",
            "operator_key": "VX"
          },
          {
            "type": "operator",
            "name": "VY",
            "operator_key": "VY"
          }
        ]
      },
      "plan": {
        "word": "VX delta VY delta",
        "instructions": [
          {
            "op": "Expand",
            "num_moments": 64
          },
          {
            "op": "ApplyOperator",
            "operator_index": 1,
            "operator_name": "VY"
          },
          {
            "op": "MirrorToBra"
          },
          {
            "op": "ApplyOperator",
            "operator_index": 0,
            "operator_name": "VX"
          },
          {
            "op": "Expand",
            "num_moments": 64
          },
          {
            "op": "Trace",
            "num_weights": 1
          }
        ]
      },
      "operators": [
        "VX",
        "VY"
      ],
      "reconstruction": {
        "shape_id": "kubo_bastin",
        "default_reconstruction": "kubo_bastin",
        "formulas": [
          "bastin"
        ],
        "kernel": "jackson",
        "broadening_meV": null,
        "fermi_energies": [],
        "num_points": null,
        "use_fft_grid": false,
        "legacy_moment_prefix": "NonEqOpVX-VY",
        "legacy_moment_ext": ".chebmom2D",
        "legacy_text_round_trip": false
      },
      "dt": null
    }
  PASS  MCP words_describe JSON == CLI words describe JSON
  MCP -> call_tool: {"tool": "inspect", "arguments": {"path": "goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"}}
  MCP <- inspect response:
    {
      "kind": "moments",
      "path": "/home/ICN2/jgarcia/lsquant-agi-build/goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
      "rank": 1,
      "dims": [
        64
      ],
      "system_size": 2048,
      "number_of_states": 0,
      "state_label": "exact",
      "word": "",
      "window": [
        -2.0,
        2.0
      ],
      "band_width": 4.0,
      "half_width": 2.0,
      "band_center": 0.0,
      "kpm_weighted": true,
      "symmetrised": false,
      "number_of_time_steps": 0,
      "time_step": 0.0
    }
  PASS  MCP inspect JSON == CLI inspect JSON

STAGE 7 SELF-CHECK PASS
```

---

## 7. Deviations and Open Questions

- **Deviations**:
  - `pip install --user mcp` on TCN installed `mcp 2.3.0`, where `FastMCP` was renamed to `MCPServer` (`from mcp.server.mcpserver import MCPServer`). `python/lsquant/mcp_server.py` imports `MCPServer` with a fallback to `FastMCP` so both `mcp 2.x` and `mcp 1.x` work transparently.
  - Per instructions, Section 8 of the brief (deleting the legacy tree `src/`, `include/`, `test/`) was not executed in this session and awaits tester confirmation of the goldens from words.
- **Open Questions**: None.
