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

## G10 (part) — `nonorth/` (rotated chain1d and Haldane, handoff 0006)

Not legacy output: generated by the developer with
`python3 python/lsquant/apps/make_nonorth_fixtures.py` (seed 12345, strength `alpha = 0.55`)
from `chain1d/operators` and `haldane/operators`: `H' = T^+ H T`, `S' = T^+ T`,
`V' = T^+ V T`, variants `real` and `complex`; each directory carries `manifest.json` with
the sha256 of every `.CSR`. The reference for these inputs is the orthogonal G1, G2, G5
(traces are basis-invariant), never a file computed from them. Verified independently by
the tester (`figures/stage_6_fixtures.py`, dense numpy/scipy, log
`figures/out/stage_6/fixture_verification.log`, handoff `0006` `tests.md` section 1): sha256
equal to the manifests; S Hermitian, positive definite, cond(S) 5.22 to 6.70 equal to the
manifests to 2e-14; `eigh(H', S')` equal to `eigvalsh(H)` to 1.5e-14; H', V' exactly
Hermitian. Accepted as committed at `a5108e3`; not regenerated.
