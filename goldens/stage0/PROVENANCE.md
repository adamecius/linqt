# Stage-0 golden masters — provenance

Every file under `goldens/stage0/` was produced by the **legacy** LSQUANT build in one
session on TCN, following `handoff/0000-goldens/codex.md`. No file was edited by hand.
Inputs (operators, `BOUNDS`, model seeds) are copied from the committed tree and listed
as such.

| Item | Value |
|---|---|
| Machine | `fatnode-tacn` (TCN), AMD EPYC 7551, 64 logical CPUs, 2 TB RAM, Ubuntu 22.04.4 LTS, kernel 5.15.0-161 |
| Compiler / tools | g++ 11.4.0 (Ubuntu 11.4.0-1ubuntu1~22.04.3), cmake 3.22.1, Python 3.10.12, numpy 2.2.6, scipy 1.15.3, matplotlib 3.10.9; Eigen 3.4.0, FFTW 3.3.8, GSL 2.7.1 |
| Repository | `adamecius/lsquant`, branch `cpp-unified`, HEAD `8ecbf836677f8d7a13f72dcdd503225a6ffac820` (clone `~/lsquant-unified`) |
| wannier2sparse | `26cab4affff8266f8a76ec73f3201e6c063a1ba6` (scripts/fetch_wannier2sparse.sh) |
| Build | `cmake -S . -B build -DLSQUANT_BUILD_NEW_CORE=OFF -DLSQUANT_BUILD_INPUT_TOOLS=ON` (CMAKE_BUILD_TYPE=Release, default flags); `cmake --build build -j 32` |
| Seed | `KPM_SEED=12345` for every command, except where a tutorial script sets its own `KPM_SEED` (noted below) |
| Threads | `OMP_NUM_THREADS=16` for every generation command (tutorial 06 `rp_validate.py`: 64, as its README says; `rp_spin.py` sets 16 itself) |
| Date | 2026-10-20 (UTC) |

## G1, G2, G3 — `chain1d/` (exact trace, `M = 64`, window `[-2, 2]`)

Produced by `KPM_SEED=12345 OMP_NUM_THREADS=16 bash scripts/regenerate_goldens.sh build`
(wannier2sparse `chain1d 2048 1 1 VX`, then `inline_compute-kpm-nonEqOp chain1d VX VX 64 exact`,
`inline_compute-kpm-spectralOp chain1d 1 64 exact`, and `test/dos_reconstructed.sh`).
All files of this regeneration under `test/golden/chain1d/` came out byte-identical to the
committed ones; they were copied here:

- `operators/chain1d.HAM.CSR`, `operators/chain1d.VX.CSR` (N = 2048), `BOUNDS` (`-2 2`)
- `SpectralOp1chain1dKPM_M64_stateexact.chebmom1D` (G1)
- `NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D` (G2)
- `DOS_reconstructed.dat` (G3: `inline_kpm-DOS-standalone chain1d 1 64 exact`, per-site, nan/inf rows dropped)
- `meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat` (G3 raw curve:
  `inline_spectralFunctionFromChebmom SpectralOp1chain1dKPM_M64_stateexact.chebmom1D 120`, Jackson 120 meV)

## G4 — `graphene/` (random phase, `KPM_SEED=12345`, `M = 64`)

Produced by the same `scripts/regenerate_goldens.sh` run (wannier2sparse
`graphene_golden 8 8 1 all`, `inline_compute-kpm-nonEqOp graphene_golden VX VX 64`,
`inline_kuboBastinFromChebmom <mom> 64`). **These files differ from the committed
`test/golden` ones** (see `handoff/0000-goldens/tests.md`, section 3); the committed
`test/golden` was restored and not changed. The files here are this session's output:

- `operators/graphene_golden.*.CSR` (17 operators, dim 256; byte-identical to the committed CSRs),
  `graphene_golden_hr.dat`, `.uc`, `.xyz` (committed wannier2sparse seed inputs)
- `NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D`
- `KuboBastin_NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefaultJACKSON.dat`
- `fft_grid/`: `inline_kuboBastinFromChebmom_FFTgrid <mom> 10`, `inline_kuboGreenwoodFromChebmom_FFTgrid <mom> 10`
- `sea_surf/`: `inline_kuboBastinSeaFromChebmom <mom> 10`, `inline_kuboBastinSurfFromChebmom <mom> 10`

Note: `scripts/regenerate_goldens.sh` deletes `test/golden/operators/` and wannier2sparse
does not write the two committed sidecars `graphene_golden.HAM.desc` and `.VX.desc`, so the
moments are computed without them (bounds from the Gershgorin + pad safeguard, header
`BandWidth = 15.0366847652085`).

## G5 — `haldane/` (exact trace, `M = 128`, bounds `[-4.5, 4.5]`, 10 meV)

Inputs: `operators/haldane.{HAM,VX,VY}.CSR` copied from `test/golden/haldane/operators`
(dim 512); `BOUNDS` written as `-4.5 4.5` (test/hall_haldane.sh). In a temporary directory
with an exact local-state file over the full basis:

```
inline_compute-kpm-nonEqOp haldane VX VY 128 exact
inline_compute-kpm-nonEqOp haldane VX VX 128 exact
inline_kuboBastinFromChebmom     NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D 10
inline_kuboBastinSeaFromChebmom  NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D 10
inline_kuboBastinSurfFromChebmom NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D 10
inline_kuboGreenwoodFromChebmom  NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D 10
```

## G6 — `chain1d_ballistic/` (exact trace, `M = 256`, 64 times, `tmax = 20`, 10 meV)

Inputs: `operators/chain1d.{HAM,VX}.CSR` from `test/golden/chain1d_ballistic/operators`
(the committed `chain1d.HAM.desc` is deliberately not copied: test/msd_ballistic.sh stages only
the two CSRs and a `BOUNDS` file); `BOUNDS` written as `-2 2`.

```
inline_compute-kpm-MeanSquareDisplacement chain1d VX 256 64 20 exact
inline_timeCorrelationsFromChebmom CorrelationVXchain1dKPM_M256_state20exact.chebmomTD 10 <E>   for E = 0 0.5 1.0 1.5
```

## G7 — `chain1d_mag/` (exact trace, `M = 32`, 96 times, `tmax = 95`)

Inputs: `operators/chain1d_mag.{HAM,SX,SY,SZ,PX}.CSR` from `test/models/chain1d_mag/operators`.

```
inline_compute-kpm-TimeEvProjetedOp chain1d_mag <OP> PX 32 96 95 exact    for OP = SX SY SZ
```

## G8 and the tutorial figure data — `tutorials/01` ... `tutorials/06`

Produced by `figures/stage_0.py` (run with `--out`/`--work` outside the repository), which runs
every README command block of `examples/0N` in a fresh copy of the tutorial (inputs only:
`*.py`, `*.json`, `README.md`, `requirements.txt`) with `LSQUANT_BIN=build`,
`PATH=build:$PATH`, `KPM_SEED=12345`, `OMP_NUM_THREADS=16`, `MPLBACKEND=Agg`. The exact
command list per tutorial is the `TUTORIALS` table in `figures/stage_0.py`.

Per tutorial directory:

- `fig_*.png`: the figures the tutorial scripts write;
- `*.dat`: curve files left in the tutorial directory by the README commands and scripts;
- `captured/NNNN_<name>`: every `.dat` / `.chebmom1D` / `.chebmom2D` / `.chebmomTD` file the
  tutorial Python scripts read (through `numpy.loadtxt` or `open`), in read order, copied
  by the capture shim of `figures/stage_0.py` before the scripts delete them;
  `captured/index.tsv` maps each to its SHA-1 and its source (`<tmpdir>/` = a
  `tempfile.TemporaryDirectory` of the script).

Parameters a README does not fix, and deviations:

- 01: README step 2 "Repeat for chain1d_N1024 and chain1d_N4096" run as a loop over the
  three labels with the same JSON fields.
- 02: `lsqloc.py` sets `KPM_SEED` to the disorder seed (1..16) for each realisation.
- 04: the README says the operators are shipped under `examples/04/operators`; that
  directory does not exist. Operators staged from `test/models/chain1d_mag/operators` and an
  exact local-state file written (as test/spin_precession.sh does).
- 05: README step 3 references a `dos.json` the README never defines; not run standalone
  (`lsqtransport.py` writes and runs it). `lsqtransport.py` sets `KPM_SEED=1`.
- 06: `rp_validate.py` with `OMP_NUM_THREADS=64` (README); `rp_spin.py` sets `KPM_SEED=1..20`
  and `OMP_NUM_THREADS=16` per random vector itself.

## G9 — `bases/graphene/`, `bases/haldane/` (stage 5; independent oracle, not the legacy build)

G9 is the only golden not produced by the legacy LSQUANT build: AGENTS.md section 8 names
`wannier2sparse/tools/hr_exactdiag.py` (exact diagonalisation on a k mesh, independent code)
as its oracle. Generated on TCN on 2026-10-21 by the stage-5 tester (handoff
`0005-bloch-and-eigen-bases`, branch `handoff/0005-codex`) with Python 3.10.12, numpy 2.2.6.

| Item | Value |
|---|---|
| wannier2sparse | `~/wannier2sparse` at `3f553db8dc2f226f44592b6ee2a717d3cd71d8c7` (local edits only in `tools/block_pilot.cpp`, `docs/FINDINGS.md`, `roadmap_plans.md`; `tools/hr_exactdiag.py`, `examples/gen_models.py` and `src/` clean). The pinned `26cab4a` of `scripts/fetch_wannier2sparse.sh` has no `tools/hr_exactdiag.py`. |
| Models | `python3 examples/gen_models.py` (graphene: t = -1 eV, a1 = (1.5, 0.8660254, 0), a2 = (1.5, -0.8660254, 0), tau_B = (1, 0, 0); Haldane: t1 = -1, t2 = 0.15, phi = pi/2). Committed unchanged as `<model>/model/<model>_hr.dat`, `.uc`, `.xyz`. |
| Velocities | `wannier2sparse <model> 1 1 1 VX VY --mode bundle -o bund` (velocity mode `berry_connection`, `-i (R.lat + tau_b - tau_a) H`, eV Angstrom); `bund/<model>.w2sp/operators/VX.hr.dat`, `VY.hr.dat` copied to `<model>/model/<model>_vx_hr.dat`, `_vy_hr.dat` (the file names `hr_exactdiag.py sigma` reads). |
| Window | eigenvalues of `H(k)` on the commensurate 48 x 48 Gamma mesh `k = n/48` (the k mesh of the 48 x 48 x 1 supercell): both models span [-3, 3] eV; padded 10 %: `[-3.3, 3.3]` eV, half-width a = 3.3 eV. |
| Broadening | Gaussian sigma = pi a / M with M = 256 (the Jackson resolution of 256 moments on that window) = 0.0404970928 eV = 40.4970928 meV. |
| Mesh | `--nk 48`: hr_exactdiag's mesh is the shifted (i + 1/2)/48 mesh, not the Gamma mesh of the supercell (part of the finite-size difference of the G9 row). |
| Grid | `--emin -3.3 --emax 3.3 --ngrid 2641` (spacing 2.5 meV). |

Commands (in a scratch copy of the generated models, `<m>` = graphene, haldane):

```
python3 ~/wannier2sparse/tools/hr_exactdiag.py dos   <m> --nk 48 --eta 0.0404970928 --ngrid 2641 --emin -3.3000000000 --emax 3.3000000000 --out goldens/stage0/bases/<m>/G9_dos
python3 ~/wannier2sparse/tools/hr_exactdiag.py sigma <m> --nk 48 --eta 0.0404970928 --ngrid 2641 --emin -3.3000000000 --emax 3.3000000000 --comp xx --out goldens/stage0/bases/<m>/G9_sigma_xx
```

Files per model: `G9_dos.json` (`energy_eV`, `dos` per cell, integral 2.000),
`G9_sigma_xx.json` (`energy_eV`, `sigma_xx_e2h`), `G9_params.txt` (window, a, M, sigma, grid,
nk, pad: read by `scripts/check_goldens.sh 5` for the window and the broadening of the
lsquant runs), `model/` (the inputs above). Units and caveats of the oracle, as used by the
gate (`figures/stage_5_table.py`): the DOS is per cell (lsquant's extensive curve is divided
by N_k = 2304); `sigma_xx_e2h` carries the prefactor 2 pi^2 / A with A = |a1_x a2_y| =
1.299 A^2 (`cmd_sigma`; the true honeycomb cell |a1 x a2| = 2.598 A^2, so for these two
models the printed values are twice the e^2/h value; `cell_area_2d` in the same file
documents this pitfall but `cmd_sigma` does not use it). The gate compares shapes through
the fixed analytic factor 2 pi / (A N_k) and prints the least-squares factor next to it.

The whole generation is `goldens/stage0/bases/gen_g9.sh <repo>` (committed next to the
files; it generates the models in a scratch directory, writes the bundle velocities, derives
the window and the broadening, and runs the two commands above per model). A second run into
a scratch tree reproduced all 16 files byte for byte.

sha256:

```
d7cbacdc4b6bdffb5bd8afcab09695cd8f5c08114b4fb20f6ad02f6706969f9b  bases/graphene/G9_dos.json
b823599cf32f762953f123f1564a6bcb793ceaa0e22d9b8292c48d3395c727fc  bases/graphene/G9_sigma_xx.json
d79065b5bd8708d2520c56f7b235846049bd87fe8dc9eef93da7cc92cb474a43  bases/graphene/G9_params.txt
1dc70fea90f60fc8340ee0ced64639a6d9390f7f636a5eea10c06dc333ab2fde  bases/graphene/model/graphene_hr.dat
0bc2d2218e7805356c364afce524055907c1fccc287ef88cff54e6e4ba022560  bases/graphene/model/graphene.uc
2779e3f2fecf3dcff54fad0817088b4b9ab536293cb7e30f44b525a13563cdf9  bases/graphene/model/graphene.xyz
b3014714d39c0354052b6d0308b6b3a92bbb0da6d0cda0ee6d032a7b443cc2b6  bases/graphene/model/graphene_vx_hr.dat
204b566f8a05dc3405219a49630b204beaed79a2d4cbba06324783951e53d8f5  bases/graphene/model/graphene_vy_hr.dat
8cce83692f3e2820899ad73ff930e2ab1aac566a3e7e5ec7d664a3bd925ce789  bases/haldane/G9_dos.json
203c94cc23096c84d01c15391ebb3fd68d2e290a3c8887e0d1f0db11385d3ddf  bases/haldane/G9_sigma_xx.json
d79065b5bd8708d2520c56f7b235846049bd87fe8dc9eef93da7cc92cb474a43  bases/haldane/G9_params.txt
d4b53b22c4e33e067ad0ca77a6ac66c0aa76d1379d2d265e0278700dcd83d706  bases/haldane/model/haldane_hr.dat
0bc2d2218e7805356c364afce524055907c1fccc287ef88cff54e6e4ba022560  bases/haldane/model/haldane.uc
2779e3f2fecf3dcff54fad0817088b4b9ab536293cb7e30f44b525a13563cdf9  bases/haldane/model/haldane.xyz
b7bdaaad6a393fc5da52d9d0109f0872dc5c86f5771c2e8a688f3523627d05b8  bases/haldane/model/haldane_vx_hr.dat
ea86f1903319c27467003d92f6e332248ae9a7dfcb991a38571c552054f18c10  bases/haldane/model/haldane_vy_hr.dat
```

## Stability

`scripts/check_goldens.sh 0 build` reran G1, G2, G5 (16 threads) and G1, G2 (1 and 64
threads) in fresh temporary directories: every file byte-identical (see
`handoff/0000-goldens/tests.md`, section 4).

## Legacy defects (evidence only)

Files under `chain1d/legacy_defects/` are legacy outputs known to be wrong. They are kept as
evidence and **are not references**: no harness reads them.

- `chain1d/legacy_defects/batch16/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D`:
  G2 recomputed with `BATCH_SIZE=16`
  (`BATCH_SIZE=16 KPM_SEED=12345 OMP_NUM_THREADS=16 inline_compute-kpm-nonEqOp chain1d VX VX 64 exact`,
  legacy build at d9becf6). It is stable at 16 and 1 threads but wrong: legacy
  `Vectors::IterateAll` does not carry the recursion across batches, so block k is shifted
  by k orders. The command, the evidence and the cause are in its `README.md`. The batched
  run of the new core is gated against G2 itself (`scripts/check_goldens.sh 1`).
