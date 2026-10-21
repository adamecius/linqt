#!/usr/bin/env bash
# Golden-master gate (AGENTS.md section 8). Runs the stage's programs on the committed
# golden inputs in a temporary directory and compares every produced file byte for byte
# with goldens/stage0/. It never writes into goldens/ and never modifies a golden.
#
# usage: scripts/check_goldens.sh <stage> [build_dir]
#
#   stage 0  legacy build (build_dir holds the inline_* drivers): reruns the G1, G2
#            (chain1d, exact trace, M = 64) and G5 (Haldane, exact trace, M = 128)
#            pipelines at OMP_NUM_THREADS=16, then G1 and G2 again at 1 and 64 threads.
#   stage 1  new core (build_dir/lsquant_core/stage1_moments, handoff 0001): G1 and G2
#            at OMP_NUM_THREADS = 1, 8, 64; G2 again with --batch 16 at 8 threads (the
#            batched result must equal the unbatched G2); random-phase determinism
#            (two runs, cmp);
#            grep for "basis|std::cout|getenv" under lsquant/core (must be empty).
#   stage 2  Python binding (build_dir/lsquant_core/_core*.so, handoff 0002): stage 1 in
#            full, then python3 -m lsquant.apps.stage2_dos in its two forms at
#            OMP_NUM_THREADS = 1 and 16 (PYTHONPATH=build_dir/lsquant_core:python):
#            (a) --moments G1 -> JACKSON.dat vs G3; (b) --operators -> .chebmom1D vs G1
#            and JACKSON.dat vs G3; grep for "delta|DOS|Kubo|basis|std::cout|getenv"
#            under lsquant/core and lsquant/io on code lines only (lines whose first
#            non-blank characters are //, * or /* are comments, excluded and listed;
#            the code-line result must be empty).
#   stage 3  rank-2 reconstruction (handoff 0003): stages 1 and 2 in full, then
#            python3 -m lsquant.apps.stage3_kubo at OMP_NUM_THREADS = 16 and 1:
#            G5 moments (Haldane VX-VY and VX-VX, window [-4.5, 4.5], exact trace,
#            512 states, M = 128) vs goldens/stage0/haldane; G5 curves bastin, sea,
#            surf (VX-VY) and greenwood (VX-VX) at 10 meV from the golden moment
#            files; G4 moments (graphene VX-VX, --window gershgorin, random phase
#            seed 12345, 1 state, M = 64); G4 curves bastin at 64 meV, sea and surf
#            at 10 meV (sea_surf/), greenwood_fft and bastin_fft at 10 meV
#            (fft_grid/); grep -rn "compat/" lsquant/core lsquant/io must be empty.
#   stage 4  time evolution (handoff 0004): stages 1 to 3 in full, then
#            python3 -m lsquant.apps.stage4_time at OMP_NUM_THREADS = 16 and 1:
#            G6 moments (chain1d_ballistic, --word msd --velocity VX, window [-2, 2],
#            exact trace, 256 states, M = 256, 64 steps, tmax 20) vs the golden
#            .chebmomTD; the four G6 curves (10 meV, E_F = 0, 0.5, 1.0, 1.5) from the
#            golden moment file and again from the produced one; G7 moments
#            (chain1d_mag, --word projected, SX, SY, SZ with PX, --window gershgorin,
#            exact trace, 128 states, M = 32, 96 steps, tmax 95); grep -rnwE
#            "time|velocity|spin|projector|bessel|gsl" (whole lowercase words) on code
#            lines of lsquant/core and lsquant/io must be empty.
#   stage 5  Bloch and eigen bases (handoff 0005): stages 1 to 4 in full, then
#            python3 -m lsquant.apps.stage5_bases on graphene and Haldane
#            (goldens/stage0/bases/<model>/model, wannier2sparse examples/gen_models.py),
#            supercell 48 48 1, M = 256, window and broadening from G9_params.txt, at
#            OMP_NUM_THREADS = 16 and 1: W = 0 real and Bloch (random phase, 32 states,
#            seed 12345), Bloch and eigen (exact trace), DOS and kuboxx; W = 1 eV Anderson
#            (disorder seed 7) real and Bloch; the FFTW convention of BlochTransform against
#            the dense B and the wannier2sparse CSR (counted); the real-path files of
#            T = 1 and 16 must be byte-identical (Bloch and eigen reported);
#            figures/stage_5_table.py prints the codex.md tolerance rows (PASS/FAIL are
#            counted, INFO is reported); grep -rni fftw lsquant/core/chebyshev lsquant/io
#            must be empty. Needs ~/wannier2sparse/build/wannier2sparse (W2S_BIN) and
#            ~/wannier2sparse/tools (W2S_TOOLS). S5_KEEP=<dir> keeps the stage-5 outputs
#            (for figures/stage_5.py). Stage "5only" runs the stage-5 part alone.
#
# Output: one line per compared file, IDENTICAL or DIFFERS; for DIFFERS the first
# differing bytes (cmp -l | head -5) and the numeric localisation of
# scripts/golden_diff.py (first differing line, element index, both values, max
# relative difference over the file). Exit code 0 only if every file is IDENTICAL
# (and, for stage >= 1, the grep is empty).
#
# KPM_SEED=12345 always. Seeds, tolerances and goldens are never changed here.
set -uo pipefail

usage() { echo "usage: scripts/check_goldens.sh <stage> [build_dir]"; exit 2; }
[ $# -ge 1 ] || usage
STAGE="$1"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${2:-$REPO/build}" 2>/dev/null && pwd)" || { echo "ERROR: build dir '${2:-build}' not found"; exit 2; }
GOLD="$REPO/goldens/stage0"
export KPM_SEED=12345
[ -d "$GOLD" ] || { echo "ERROR: $GOLD missing"; exit 2; }

TMP="$(mktemp -d "${TMPDIR:-/tmp}/check_goldens.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
NFAIL=0
NCMP=0

# compare <produced file> <golden file> <tag>
compare() {
  local new="$1" gold="$2" tag="$3"
  NCMP=$((NCMP + 1))
  if [ ! -f "$new" ]; then
    echo "DIFFERS   $tag  (not produced: $(basename "$new"))"; NFAIL=$((NFAIL + 1)); return
  fi
  if [ ! -f "$gold" ]; then
    echo "DIFFERS   $tag  (no golden: ${gold#$REPO/})"; NFAIL=$((NFAIL + 1)); return
  fi
  if cmp -s "$new" "$gold"; then
    echo "IDENTICAL $tag"
  else
    NFAIL=$((NFAIL + 1))
    echo "DIFFERS   $tag"
    echo "          cmp: $(cmp "$new" "$gold" 2>&1 | head -1)"
    cmp -l "$new" "$gold" 2>/dev/null | head -5 | sed 's/^/          cmp -l: /'
    echo "          $(python3 "$REPO/scripts/golden_diff.py" "$new" "$gold" 2>&1)"
  fi
}

# exact-trace state file over the full basis of a CSR (legacy "local" state set)
write_exact() { local n; n=$(head -1 "$1" | awk '{print $1}'); { echo local; echo "$n"; seq 0 $((n - 1)); } > "$2"; }

# ---------------------------------------------------------------- stage 0 (legacy)
# G1 + G2: the chain1d part of scripts/regenerate_goldens.sh, on the golden operators.
legacy_chain1d() {   # <threads> <workdir>
  local T="$1" W="$2"
  mkdir -p "$W/operators"
  cp "$GOLD/chain1d/operators/"*.CSR "$W/operators/"
  cp "$GOLD/chain1d/BOUNDS" "$W/BOUNDS"
  write_exact "$W/operators/chain1d.HAM.CSR" "$W/exact"
  ( cd "$W"
    OMP_NUM_THREADS=$T "$BUILD/inline_compute-kpm-nonEqOp"   chain1d VX VX 64 exact >nonEq.log 2>&1
    OMP_NUM_THREADS=$T "$BUILD/inline_compute-kpm-spectralOp" chain1d 1 64 exact   >spectral.log 2>&1 )
  compare "$W/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
          "$GOLD/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
          "[T=$T] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)"
  compare "$W/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "[T=$T] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)"
}

# G5: Haldane, exact trace, M = 128, bounds [-4.5, 4.5] (test/hall_haldane.sh), 10 meV.
legacy_haldane() {   # <threads> <workdir>
  local T="$1" W="$2" f
  mkdir -p "$W/operators"
  cp "$GOLD/haldane/operators/"haldane.{HAM,VX,VY}.CSR "$W/operators/"
  cp "$GOLD/haldane/BOUNDS" "$W/BOUNDS"
  write_exact "$W/operators/haldane.HAM.CSR" "$W/exact"
  ( cd "$W"
    export OMP_NUM_THREADS=$T
    "$BUILD/inline_compute-kpm-nonEqOp" haldane VX VY 128 exact >xy.log 2>&1
    "$BUILD/inline_compute-kpm-nonEqOp" haldane VX VX 128 exact >xx.log 2>&1
    XY=NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
    XX=NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
    "$BUILD/inline_kuboBastinFromChebmom"     "$XY" 10 >bastin.log 2>&1
    "$BUILD/inline_kuboBastinSeaFromChebmom"  "$XY" 10 >sea.log 2>&1
    "$BUILD/inline_kuboBastinSurfFromChebmom" "$XY" 10 >surf.log 2>&1
    "$BUILD/inline_kuboGreenwoodFromChebmom"  "$XX" 10 >greenwood.log 2>&1 )
  for f in NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D \
           NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D; do
    compare "$W/$f" "$GOLD/haldane/$f" "[T=$T] haldane/$f (G5)"
  done
  for f in "$GOLD"/haldane/*.dat; do
    compare "$W/$(basename "$f")" "$f" "[T=$T] haldane/$(basename "$f") (G5)"
  done
}

stage0() {
  for b in inline_compute-kpm-nonEqOp inline_compute-kpm-spectralOp inline_kuboBastinFromChebmom \
           inline_kuboBastinSeaFromChebmom inline_kuboBastinSurfFromChebmom inline_kuboGreenwoodFromChebmom; do
    [ -x "$BUILD/$b" ] || { echo "ERROR: legacy driver $BUILD/$b missing (build with the legacy tree)"; exit 2; }
  done
  echo "== stage 0: legacy pipelines vs goldens/stage0 (KPM_SEED=$KPM_SEED, build $BUILD)"
  legacy_chain1d 16 "$TMP/g12_T16"
  legacy_haldane 16 "$TMP/g5_T16"
  legacy_chain1d 1  "$TMP/g12_T1"
  legacy_chain1d 64 "$TMP/g12_T64"
}

# ---------------------------------------------------------------- stage 1 (new core)
stage1() {
  local EXE="$BUILD/lsquant_core/stage1_moments" T W w R
  [ -x "$EXE" ] || EXE="$BUILD/lsquant/stage1_moments"
  [ -x "$EXE" ] || { echo "ERROR: stage1_moments not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  echo "== stage 1: $EXE vs goldens/stage0/chain1d (KPM_SEED=$KPM_SEED)"
  for T in 1 8 64; do
    W="$TMP/s1_T$T"; mkdir -p "$W"
    for w in dos kubo; do
      OMP_NUM_THREADS=$T "$EXE" --operators "$GOLD/chain1d/operators" --label chain1d --window -2 2 \
        --trace exact --states 2048 --moments 64 --word $w --out "$W" >"$W/$w.log" 2>&1 \
        || echo "          [T=$T] stage1_moments --word $w exited non-zero: $(tail -1 "$W/$w.log")"
    done
    compare "$W/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
            "$GOLD/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
            "[T=$T] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)"
    compare "$W/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
            "$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
            "[T=$T] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)"
  done
  # G2 with batching (batchSize 16 < M = 64), T = 8: must equal the unbatched G2
  W="$TMP/s1_batch16"; mkdir -p "$W"
  OMP_NUM_THREADS=8 "$EXE" --operators "$GOLD/chain1d/operators" --label chain1d --window -2 2 \
    --trace exact --seed 12345 --states 2048 --moments 64 --word kubo --batch 16 --out "$W" >"$W/kubo.log" 2>&1 \
    || echo "          [T=8 batch 16] stage1_moments --word kubo exited non-zero: $(tail -1 "$W/kubo.log")"
  compare "$W/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "[T=8, batch 16] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2, batched)"
  # random-phase determinism: no golden; two identical runs are the criterion
  for R in 1 2; do
    W="$TMP/s1_random$R"; mkdir -p "$W"
    for w in dos kubo; do
      OMP_NUM_THREADS=8 "$EXE" --operators "$GOLD/chain1d/operators" --label chain1d --window -2 2 \
        --trace random --seed 12345 --states 4 --moments 64 --word $w --out "$W" >"$W/$w.log" 2>&1 \
        || echo "          [random run $R] stage1_moments --word $w exited non-zero: $(tail -1 "$W/$w.log")"
    done
  done
  local n=0 f
  for f in "$TMP/s1_random1"/*.chebmom*; do
    [ -f "$f" ] || continue; n=$((n + 1))
    compare "$TMP/s1_random2/$(basename "$f")" "$f" "[random, run2 vs run1] $(basename "$f")"
  done
  [ "$n" -gt 0 ] || { echo "DIFFERS   [random] no .chebmom output produced"; NFAIL=$((NFAIL + 1)); }
}

core_grep() {
  local hits
  if [ ! -d "$REPO/lsquant/core" ]; then
    echo "GREP      lsquant/core missing"; NFAIL=$((NFAIL + 1)); return
  fi
  hits="$(cd "$REPO" && grep -rn "basis\|std::cout\|getenv" lsquant/core)"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rn \"basis\\|std::cout\\|getenv\" lsquant/core: empty"
  else
    echo "GREP      grep -rn \"basis\\|std::cout\\|getenv\" lsquant/core: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
}

# ---------------------------------------------------------------- stage 2 (Python binding)
G1F=SpectralOp1chain1dKPM_M64_stateexact.chebmom1D
G3F=meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat

stage2() {
  local PYDIR="$BUILD/lsquant_core" T W rc
  ls "$PYDIR"/_core*.so >/dev/null 2>&1 || PYDIR="$BUILD/lsquant"
  ls "$PYDIR"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  echo "== stage 2: python3 -m lsquant.apps.stage2_dos, PYTHONPATH=$PYDIR:python, vs goldens/stage0/chain1d (KPM_SEED=$KPM_SEED)"
  for T in 1 16; do
    # form (a): reconstruct the G3 curve from the committed G1 moment file
    W="$TMP/s2_T${T}_a"; mkdir -p "$W"
    ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR:$REPO/python" python3 -m lsquant.apps.stage2_dos \
        --moments "$GOLD/chain1d/$G1F" --broadening-meV 120 --out "$W" ) >"$W/run.log" 2>&1; rc=$?
    [ $rc -eq 0 ] || echo "          [T=$T] stage2_dos form (a) exited $rc: $(tail -1 "$W/run.log")"
    sed "s|$W/||; s|^|          [T=$T] (a) printed: |" "$W/run.log"
    compare "$W/$G3F" "$GOLD/chain1d/$G3F" "[T=$T] (a) --moments G1 -> chain1d/$G3F (G3)"
    # form (b): moments through MomentEngine from the golden operators, then reconstruct
    W="$TMP/s2_T${T}_b"; mkdir -p "$W"
    ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR:$REPO/python" python3 -m lsquant.apps.stage2_dos \
        --operators "$GOLD/chain1d/operators" --label chain1d --window -2 2 --trace exact --states 2048 \
        --moments 64 --broadening-meV 120 --out "$W" ) >"$W/run.log" 2>&1; rc=$?
    [ $rc -eq 0 ] || echo "          [T=$T] stage2_dos form (b) exited $rc: $(tail -1 "$W/run.log")"
    sed "s|$W/||; s|^|          [T=$T] (b) printed: |" "$W/run.log"
    compare "$W/$G1F" "$GOLD/chain1d/$G1F" "[T=$T] (b) --operators -> chain1d/$G1F (G1)"
    compare "$W/$G3F" "$GOLD/chain1d/$G3F" "[T=$T] (b) --operators -> chain1d/$G3F (G3)"
  done
}

# physics names and side channels in the engine, on code lines only: a line whose first
# non-blank characters are //, * or /* is a comment (e.g. the notation delta_{m0}) and is
# excluded from the verdict, but listed below it.
stage2_grep() {
  local all hits coms pat='delta\|DOS\|Kubo\|basis\|std::cout\|getenv'
  local com='^[^:]+:[0-9]+:[[:space:]]*(//|\*|/\*)'
  if [ ! -d "$REPO/lsquant/core" ] || [ ! -d "$REPO/lsquant/io" ]; then
    echo "GREP      lsquant/core or lsquant/io missing"; NFAIL=$((NFAIL + 1)); return
  fi
  all="$(cd "$REPO" && grep -rn "$pat" lsquant/core lsquant/io)"
  hits="$(printf '%s\n' "$all" | grep -Ev "$com" | grep -v '^$')"
  coms="$(printf '%s\n' "$all" | grep -E "$com")"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rn \"$pat\" lsquant/core lsquant/io, code lines: empty"
  else
    echo "GREP      grep -rn \"$pat\" lsquant/core lsquant/io, code lines: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
  [ -z "$coms" ] || echo "$coms" | sed 's/^/          excluded comment line: /'
}

# ---------------------------------------------------------------- stage 3 (rank-2 reconstruction)
# One stage3_kubo invocation in its own directory; prints the app's output lines.
run_s3() {   # <threads> <workdir> <tag> <args...>
  local T="$1" W="$2" tag="$3" rc; shift 3
  mkdir -p "$W"
  ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR3:$REPO/python" python3 -m lsquant.apps.stage3_kubo \
      "$@" --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
  [ $rc -eq 0 ] || echo "          [T=$T] $tag exited $rc: $(tail -1 "$W/err.log")"
  sed "s|$W/||; s|^|          [T=$T] $tag printed: |" "$W/run.log"
}

G5XY=NonEqOpVX-VYhaldaneKPM_M128x128_stateexact
G5XX=NonEqOpVX-VXhaldaneKPM_M128x128_stateexact
G4XX=NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault

stage3() {
  local T W G
  PYDIR3="$BUILD/lsquant_core"
  ls "$PYDIR3"/_core*.so >/dev/null 2>&1 || PYDIR3="$BUILD/lsquant"
  ls "$PYDIR3"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  echo "== stage 3: python3 -m lsquant.apps.stage3_kubo, PYTHONPATH=$PYDIR3:python, vs goldens/stage0/{haldane,graphene}"
  for T in 16 1; do
    # 1. G5 moments: Haldane, exact trace over the 512 sites, M = 128, window [-4.5, 4.5]
    W="$TMP/s3_T${T}_g5m"
    run_s3 $T "$W" "G5 VX-VY moments" --operators "$GOLD/haldane/operators" --label haldane --left VX --right VY \
      --window -4.5 4.5 --trace exact --states 512 --moments 128
    run_s3 $T "$W" "G5 VX-VX moments" --operators "$GOLD/haldane/operators" --label haldane --left VX --right VX \
      --window -4.5 4.5 --trace exact --states 512 --moments 128
    compare "$W/$G5XY.chebmom2D" "$GOLD/haldane/$G5XY.chebmom2D" "[T=$T] haldane/$G5XY.chebmom2D (G5 moments)"
    compare "$W/$G5XX.chebmom2D" "$GOLD/haldane/$G5XX.chebmom2D" "[T=$T] haldane/$G5XX.chebmom2D (G5 moments)"
    # 2. G5 curves from the committed golden moment files, 10 meV
    W="$TMP/s3_T${T}_g5c"
    run_s3 $T "$W" "G5 bastin"    --moments "$GOLD/haldane/$G5XY.chebmom2D" --formula bastin    --broadening-meV 10
    run_s3 $T "$W" "G5 sea"       --moments "$GOLD/haldane/$G5XY.chebmom2D" --formula sea       --broadening-meV 10
    run_s3 $T "$W" "G5 surf"      --moments "$GOLD/haldane/$G5XY.chebmom2D" --formula surf      --broadening-meV 10
    run_s3 $T "$W" "G5 greenwood" --moments "$GOLD/haldane/$G5XX.chebmom2D" --formula greenwood --broadening-meV 10
    for G in KuboBastin_$G5XY KuboBastinSea_$G5XY KuboBastinSurf_$G5XY KuboGreenWood_$G5XX; do
      compare "$W/${G}JACKSON.dat" "$GOLD/haldane/${G}JACKSON.dat" "[T=$T] haldane/${G}JACKSON.dat (G5 curve)"
    done
    # 3. G4 moments: graphene, random phase seed 12345, 1 state, M = 64, Gershgorin window (pad 0.10)
    W="$TMP/s3_T${T}_g4m"
    run_s3 $T "$W" "G4 VX-VX moments" --operators "$GOLD/graphene/operators" --label graphene_golden --left VX --right VX \
      --window gershgorin --trace random --seed 12345 --states 1 --moments 64
    compare "$W/$G4XX.chebmom2D" "$GOLD/graphene/$G4XX.chebmom2D" "[T=$T] graphene/$G4XX.chebmom2D (G4 moments)"
    # 4. G4 curves from the committed golden moment file
    W="$TMP/s3_T${T}_g4c"
    run_s3 $T "$W" "G4 bastin 64 meV"     --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula bastin        --broadening-meV 64
    run_s3 $T "$W" "G4 sea 10 meV"        --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula sea           --broadening-meV 10
    run_s3 $T "$W" "G4 surf 10 meV"       --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula surf          --broadening-meV 10
    run_s3 $T "$W" "G4 greenwood_fft 10 meV" --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula greenwood_fft --broadening-meV 10
    run_s3 $T "$W" "G4 bastin_fft 10 meV" --moments "$GOLD/graphene/$G4XX.chebmom2D" --formula bastin_fft    --broadening-meV 10
    compare "$W/KuboBastin_${G4XX}JACKSON.dat" "$GOLD/graphene/KuboBastin_${G4XX}JACKSON.dat" \
            "[T=$T] graphene/KuboBastin_${G4XX}JACKSON.dat (G4 curve, 64 meV)"
    compare "$W/KuboBastinSea_${G4XX}JACKSON.dat" "$GOLD/graphene/sea_surf/KuboBastinSea_${G4XX}JACKSON.dat" \
            "[T=$T] graphene/sea_surf/KuboBastinSea_${G4XX}JACKSON.dat (G4 curve, 10 meV)"
    compare "$W/KuboBastinSurf_${G4XX}JACKSON.dat" "$GOLD/graphene/sea_surf/KuboBastinSurf_${G4XX}JACKSON.dat" \
            "[T=$T] graphene/sea_surf/KuboBastinSurf_${G4XX}JACKSON.dat (G4 curve, 10 meV)"
    compare "$W/KuboGreenWood_${G4XX}JACKSON_FFTgrid.dat" "$GOLD/graphene/fft_grid/KuboGreenWood_${G4XX}JACKSON_FFTgrid.dat" \
            "[T=$T] graphene/fft_grid/KuboGreenWood_${G4XX}JACKSON_FFTgrid.dat (G4 curve, 10 meV)"
    compare "$W/KuboBastin_${G4XX}JACKSON_FFTgrid.dat" "$GOLD/graphene/fft_grid/KuboBastin_${G4XX}JACKSON_FFTgrid.dat" \
            "[T=$T] graphene/fft_grid/KuboBastin_${G4XX}JACKSON_FFTgrid.dat (G4 curve, 10 meV)"
  done
}

# lsquant/compat holds legacy formulas for the gate only; nothing in core/ or io/ may include it.
stage3_grep() {
  local hits
  hits="$(cd "$REPO" && grep -rn "compat/" lsquant/core lsquant/io)"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rn \"compat/\" lsquant/core lsquant/io: empty"
  else
    echo "GREP      grep -rn \"compat/\" lsquant/core lsquant/io: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
}

# ---------------------------------------------------------------- stage 4 (time evolution)
# One stage4_time invocation in its own directory; prints the app's output lines.
run_s4() {   # <threads> <workdir> <tag> <args...>
  local T="$1" W="$2" tag="$3" rc; shift 3
  mkdir -p "$W"
  ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR4:$REPO/python" python3 -m lsquant.apps.stage4_time \
      "$@" --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
  [ $rc -eq 0 ] || echo "          [T=$T] $tag exited $rc: $(tail -1 "$W/err.log")"
  sed "s|$W/||; s|^|          [T=$T] $tag printed: |" "$W/run.log"
  : >"$W/run.log"
}

G6M=CorrelationVXchain1dKPM_M256_state20exact
G6EF="0.000000 0.500000 1.000000 1.500000"

stage4() {
  local T W ef OP s e
  PYDIR4="$BUILD/lsquant_core"
  ls "$PYDIR4"/_core*.so >/dev/null 2>&1 || PYDIR4="$BUILD/lsquant"
  ls "$PYDIR4"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  echo "== stage 4: python3 -m lsquant.apps.stage4_time, PYTHONPATH=$PYDIR4:python, vs goldens/stage0/{chain1d_ballistic,chain1d_mag}"
  for T in 16 1; do
    # 1. G6 moments: chain1d ballistic MSD, VX, window [-2, 2], exact trace 256 states, M = 256, 64 steps, tmax 20
    W="$TMP/s4_T${T}_g6m"; s=$(date +%s)
    run_s4 $T "$W" "G6 msd moments" --operators "$GOLD/chain1d_ballistic/operators" --label chain1d --velocity VX \
      --window -2 2 --trace exact --states 256 --moments 256 --steps 64 --tmax 20 --word msd
    e=$(date +%s); echo "          [T=$T] G6 msd moments wall $((e - s)) s"
    compare "$W/$G6M.chebmomTD" "$GOLD/chain1d_ballistic/$G6M.chebmomTD" "[T=$T] chain1d_ballistic/$G6M.chebmomTD (G6 moments)"
    # 1b. G6 curves (10 meV, E_F = 0, 0.5, 1.0, 1.5) from the committed golden moment file
    W="$TMP/s4_T${T}_g6c_gold"
    run_s4 $T "$W" "G6 curves (golden moments)" --moments "$GOLD/chain1d_ballistic/$G6M.chebmomTD" \
      --broadening-meV 10 --fermi-energy 0 0.5 1.0 1.5
    for ef in $G6EF; do
      compare "$W/mean${G6M}EF${ef}JACKSON.dat" "$GOLD/chain1d_ballistic/mean${G6M}EF${ef}JACKSON.dat" \
              "[T=$T] chain1d_ballistic/mean${G6M}EF${ef}JACKSON.dat (G6 curve from golden moments)"
    done
    # 1c. the same four curves from the moment file just produced in step 1
    W="$TMP/s4_T${T}_g6c_new"
    run_s4 $T "$W" "G6 curves (new moments)" --moments "$TMP/s4_T${T}_g6m/$G6M.chebmomTD" \
      --broadening-meV 10 --fermi-energy 0 0.5 1.0 1.5
    for ef in $G6EF; do
      compare "$W/mean${G6M}EF${ef}JACKSON.dat" "$GOLD/chain1d_ballistic/mean${G6M}EF${ef}JACKSON.dat" \
              "[T=$T] chain1d_ballistic/mean${G6M}EF${ef}JACKSON.dat (G6 curve from new moments)"
    done
    # 2. G7 moments: chain1d_mag, S_a with projector PX, Gershgorin window (legacy SpectralBounds
    #    without BOUNDS file: [-2.31, 2.31]), exact trace 128 states, M = 32, 96 steps, tmax 95
    W="$TMP/s4_T${T}_g7m"
    for OP in SX SY SZ; do
      s=$(date +%s)
      run_s4 $T "$W" "G7 $OP-PX moments" --operators "$GOLD/chain1d_mag/operators" --label chain1d_mag --operator $OP \
        --projector PX --window gershgorin --trace exact --states 128 --moments 32 --steps 96 --tmax 95 --word projected
      e=$(date +%s); echo "          [T=$T] G7 $OP-PX moments wall $((e - s)) s"
    done
    for OP in SX SY SZ; do
      compare "$W/TimeProj$OP-PXchain1d_magKPM_M32_stateexact.chebmomTD" \
              "$GOLD/chain1d_mag/TimeProj$OP-PXchain1d_magKPM_M32_stateexact.chebmomTD" \
              "[T=$T] chain1d_mag/TimeProj$OP-PXchain1d_magKPM_M32_stateexact.chebmomTD (G7 moments)"
    done
  done
}

# physics names in the engine (stage 4), on code lines only, whole lowercase words; the
# identifiers TimeAxis, EvolveRegister, EvolveDerived, ApplyDerivedSeries are allowed (they
# are not lowercase whole words, so -w with a case-sensitive pattern never matches them).
stage4_grep() {
  local all hits coms pat='time|velocity|spin|projector|bessel|gsl'
  local com='^[^:]+:[0-9]+:[[:space:]]*(//|\*|/\*)'
  all="$(cd "$REPO" && grep -rnwE "$pat" lsquant/core lsquant/io)"
  hits="$(printf '%s\n' "$all" | grep -Ev "$com" | grep -v '^$')"
  coms="$(printf '%s\n' "$all" | grep -E "$com")"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rnwE \"$pat\" lsquant/core lsquant/io, code lines: empty"
  else
    echo "GREP      grep -rnwE \"$pat\" lsquant/core lsquant/io, code lines: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
  [ -z "$coms" ] || echo "$coms" | sed 's/^/          excluded comment line: /'
}

# ---------------------------------------------------------------- stage 5 (Bloch and eigen bases)
# G9: graphene and Haldane (wannier2sparse examples/gen_models.py, committed under
# goldens/stage0/bases/<model>/model/), supercell 48 x 48 x 1, M = 256, window and
# broadening from goldens/stage0/bases/<model>/G9_params.txt (eigenvalue window of the
# 48 x 48 Gamma mesh padded 10 %, broadening = Jackson resolution pi * a / M).
S5N=48
S5_W2S="${W2S_BIN:-$HOME/wannier2sparse/build/wannier2sparse}"

# Python helpers of the stage-5 gate, written into $TMP (they never touch the repository).
s5_helpers() {
  mkdir -p "$TMP/s5py"
  # 1. format adapter for the stage-5 reader (python/lsquant/models/wannier.py), which
  #    requires (a) the dense Wannier90 _hr.dat layout (every R carries all W*W entries;
  #    gen_models.py writes only the non-zero hoppings) and (b) a comment line after the
  #    orbital count in .xyz (gen_models.py and wannier2sparse have none). Same model:
  #    s5_conv.py checks the lsquant CSR built from the copy against the wannier2sparse CSR
  #    built from the committed originals.
  #    usage: s5_densify.py <src_hr> <dst_hr> <src_xyz> <dst_xyz>
  cat >"$TMP/s5py/s5_densify.py" <<'PY'
import sys, collections
src, dst, sx, dx = sys.argv[1:5]
X = open(sx).read().split('\n')
open(dx, 'w').write('\n'.join([X[0], 'orbital positions (Angstrom), comment line added for the stage-5 gate'] + X[1:]))
L = [l for l in open(src).read().split('\n')]
nw = int(L[1]); nr = int(L[2]); i = 3; deg = []
while len(deg) < nr:
    deg += L[i].split(); i += 1
acc = collections.OrderedDict()
for l in L[i:]:
    s = l.split()
    if len(s) < 7: continue
    R = (int(s[0]), int(s[1]), int(s[2]))
    acc.setdefault(R, {})[(int(s[3]), int(s[4]))] = (s[5], s[6])
assert all(d == '1' for d in deg), "densify: ndegen != 1 not supported"
with open(dst, 'w') as f:
    f.write(L[0] + ' (dense copy for the stage-5 gate)\n%d\n%d\n' % (nw, len(acc)))
    f.write(' '.join('1' for _ in acc) + '\n')
    for R, d in acc.items():
        for b in range(1, nw + 1):
            for a in range(1, nw + 1):
                re, im = d.get((a, b), ('0.0000000000', '0.0000000000'))
                f.write('%d %d %d  %d %d  %s %s\n' % (R[0], R[1], R[2], a, b, re, im))
PY
  # 2. stage5_bases with the Anderson disorder seed forced to 7 (the CLI draws the disorder
  #    with the same --seed as the random-phase states; the brief's --disorder-seed is missing)
  #    (the helper names start with s5_ so that none shadows a standard-library module)
  cat >"$TMP/s5py/s5_dis.py" <<'PY'
import sys
import lsquant.bases as bases
from lsquant.apps import stage5_bases
_orig = bases.make_anderson_disorder
def _seed7(dimension, disorder_w, seed=12345):
    return _orig(dimension, disorder_w, 7)
bases.make_anderson_disorder = _seed7
sys.exit(stage5_bases.main(sys.argv[1:]))
PY
  # 3. FFTW convention of BlochTransform against the dense B (6 x 4 x 1, non-square on
  #    purpose) and against the wannier2sparse supercell CSR at 6 x 4 x 1 and 48 x 48 x 1;
  #    prints PASS/FAIL rows (counted by the gate) and INFO rows (localisation)
  cat >"$TMP/s5py/s5_conv.py" <<'PY'
import sys, numpy as np
from lsquant import _core
from lsquant.models.wannier import load_wannier_model, read_uc, read_xyz
mdir, seed, w2s_dir, w2s_small = sys.argv[1:5]
m = load_wannier_model(mdir, seed=seed)
def gate(ok, text): print(("PASS  " if ok else "FAIL  ") + text)
uc = read_uc(mdir + "/%s.uc" % seed); tau = read_xyz(mdir + "/%s.xyz" % seed)
rng = np.random.default_rng(1)
def rv(n): return rng.standard_normal(n) + 1j * rng.standard_normal(n)
def rel(a, b): return np.max(np.abs(a - b)) / np.max(np.abs(b))
# dense B on a 6 x 4 x 1 supercell, index (R, alpha) / (k, alpha) flat iR * W + alpha
N = (6, 4, 1); W = m.num_wann; Nk = 24
k = m.k_mesh(N); ph = m.bloch_phases(k)
Bf = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.RealToBloch)
Bb = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.BlochToReal)
n = np.array([(a, b, c) for a in range(N[0]) for b in range(N[1]) for c in range(N[2])], float)
Rc = n @ uc
for sgn in (-1, +1):
    Bd = np.zeros((Nk * W, Nk * W), complex)
    for al in range(W):
        Bd[al::W, al::W] = np.exp(sgn * 1j * (k @ (Rc + tau[al]).T)) / np.sqrt(Nk)
    xs = [rv(Nk * W) for _ in range(4)]
    e_f = max(rel(Bf.multiply(x), Bd @ x) for x in xs)
    e_b = max(rel(Bb.multiply(x), Bd.conj().T @ x) for x in xs)
    tag = "implemented sign, gated" if sgn < 0 else "brief.md sign, reported"
    line = ("CONV %s 6x4x1 dense B = exp(%si k.(R+tau))/sqrt(Nk) (%s): RealToBloch vs B x %.3e, BlochToReal vs B^+ x %.3e"
            % (seed, '+' if sgn > 0 else '-', tag, e_f, e_b))
    if sgn < 0:
        gate(max(e_f, e_b) <= 1e-12, line + " (tol 1e-12)")
    else:
        print("INFO  " + line)
def csr_checks(N, d, Hk=None, Vk=None, Bf=None, Bb=None):
    """lsquant supercell CSR vs the wannier2sparse CSR of the same model (brief: the real-space
    basis IS the wannier2sparse CSR) and, at 48 x 48, B H_w2s B^+ = H_k, B VX_w2s B^+ = v_x."""
    dim = N[0] * N[1] * N[2] * W
    Hr = _core.SparseMatrixOperator(_core.ReadCsrText(d + "/%s.HAM.CSR" % seed))
    Vr = _core.SparseMatrixOperator(_core.ReadCsrText(d + "/%s.VX.CSR" % seed))
    Hp = _core.SparseMatrixOperator(m.hamiltonian_csr(N)); Vp = _core.SparseMatrixOperator(m.velocity_csr(N, "VX"))
    xs = [rv(dim) for _ in range(4)]
    s = "%dx%dx%d" % N
    eh = max(rel(Hp.multiply(x), Hr.multiply(x)) for x in xs); ev = max(rel(Vp.multiply(x), Vr.multiply(x)) for x in xs)
    gate(max(eh, ev) <= 1e-12, "CONV %s %s lsquant supercell CSR vs wannier2sparse CSR: H %.3e, VX %.3e (tol 1e-12)" % (seed, s, eh, ev))
    # localisation: the same comparison with the lsquant cell index permuted to n1 fastest (iR = n1 + N1 n2)
    perm = np.array([((n1 + N[0] * n2) * W + a) for n1 in range(N[0]) for n2 in range(N[1]) for a in range(W)])
    def pm(op, x):            # P^T op P x with P mapping row-major (lsquant) to n1-fastest (wannier2sparse)
        y = np.empty_like(x); xp = np.empty_like(x); xp[np.arange(dim)] = x[perm]
        z = op.multiply(xp); y[perm] = z; return y
    eh2 = max(rel(pm(Hp, x), Hr.multiply(x)) for x in xs); ev2 = max(rel(pm(Vp, x), Vr.multiply(x)) for x in xs)
    print("INFO  CONV %s %s localisation: lsquant CSR with the cell index renumbered n1-fastest (n1 + N1 n2) vs wannier2sparse CSR: H %.3e, VX %.3e"
          % (seed, s, eh2, ev2))
    if Hk is not None:
        e1 = max(rel(Bf.multiply(Hr.multiply(Bb.multiply(x))), Hk.multiply(x)) for x in xs)
        e2 = max(rel(Bf.multiply(Vr.multiply(Bb.multiply(x))), Vk.multiply(x)) for x in xs)
        gate(e1 <= 1e-10, "CONV %s %s B H_w2s B^+ x vs H_k x: %.3e (tol 1e-10)" % (seed, s, e1))
        gate(e2 <= 1e-10, "CONV %s %s B VX_w2s B^+ x vs v_x(k) x: %.3e (tol 1e-10, brief.md section 3)" % (seed, s, e2))
csr_checks((6, 4, 1), w2s_small)
# 48 x 48: B^+ B = 1, B H_w2s B^+ = H_k, B VX_w2s B^+ = v_x with the wannier2sparse CSR
N = (48, 48, 1); Nk = 48 * 48
k = m.k_mesh(N); ph = m.bloch_phases(k)
Bf = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.RealToBloch)
Bb = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.BlochToReal)
Hk = _core.BlockDiagonalOperator(m.Hk(k)); Vk = _core.BlockDiagonalOperator(m.dHk(k, "VX"))
xs = [rv(Nk * W) for _ in range(4)]
e = max(rel(Bb.multiply(Bf.multiply(x)), x) for x in xs)
gate(e <= 1e-12, "CONV %s 48x48x1 max|B^+ B x - x|/max|x| = %.3e (tol 1e-12)" % (seed, e))
csr_checks(N, w2s_dir, Hk, Vk, Bf, Bb)
PY
}

# One stage5_bases run in its own directory (wall time printed)
run_s5() {   # <threads> <workdir> <tag> <entry: cli|dis> <args...>
  local T="$1" W="$2" tag="$3" entry="$4" rc s e; shift 4
  mkdir -p "$W"; s=$(date +%s.%N)
  if [ "$entry" = cli ]; then
    ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR5:$REPO/python" python3 -m lsquant.apps.stage5_bases \
        "$@" --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
  else
    ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR5:$REPO/python" python3 "$TMP/s5py/s5_dis.py" \
        "$@" --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
  fi
  e=$(date +%s.%N)
  [ $rc -eq 0 ] || { echo "          [T=$T] $tag exited $rc: $(tail -1 "$W/err.log")"; NFAIL=$((NFAIL + 1)); }
  printf "          [T=%s] %-34s wall %7.2f s  -> %s\n" "$T" "$tag" "$(echo "$e - $s" | bc)" "$(ls "$W" | grep -v '\.log$' | tr '\n' ' ')"
}

stage5() {
  local T m w W P LO HI A SIG BMEV LO1 HI1 B1 C rel f g
  PYDIR5="$BUILD/lsquant_core"
  ls "$PYDIR5"/_core*.so >/dev/null 2>&1 || PYDIR5="$BUILD/lsquant"
  ls "$PYDIR5"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  [ -x "$S5_W2S" ] || { echo "ERROR: wannier2sparse binary $S5_W2S missing (set W2S_BIN)"; exit 2; }
  echo "== stage 5: python3 -m lsquant.apps.stage5_bases, PYTHONPATH=$PYDIR5:python, vs goldens/stage0/bases (G9)"
  s5_helpers
  for m in graphene haldane; do
    P="$GOLD/bases/$m/G9_params.txt"
    LO=$(awk '$1=="window_lower"{print $2}' "$P"); HI=$(awk '$1=="window_upper"{print $2}' "$P")
    A=$(awk '$1=="half_width"{print $2}' "$P");    SIG=$(awk '$1=="sigma_eV"{print $2}' "$P")
    BMEV=$(python3 -c "print('%.7f' % ($SIG * 1000))")
    # W = 1 eV: window half-width widened by W/2 before the 10 % pad, broadening = pi * a_W / M
    read LO1 HI1 B1 <<<"$(python3 -c "import math; a=($A/1.10+0.5)*1.10; print('%.10f %.10f %.7f' % (-a, a, math.pi*a/256*1000))")"
    mkdir -p "$TMP/s5models/$m" "$TMP/s5w2s/$m"
    python3 "$TMP/s5py/s5_densify.py" "$GOLD/bases/$m/model/${m}_hr.dat" "$TMP/s5models/$m/${m}_hr.dat" \
      "$GOLD/bases/$m/model/$m.xyz" "$TMP/s5models/$m/$m.xyz"
    cp "$GOLD/bases/$m/model/$m.uc" "$TMP/s5models/$m/"
    mkdir -p "$TMP/s5w2s/${m}_6x4"
    ( cd "$GOLD/bases/$m/model" && "$S5_W2S" $m $S5N $S5N 1 VX -o "$TMP/s5w2s/$m" --quiet --no-log-file \
      && "$S5_W2S" $m 6 4 1 VX -o "$TMP/s5w2s/${m}_6x4" --quiet --no-log-file ) >"$TMP/s5w2s/$m.log" 2>&1 \
      || echo "          wannier2sparse $m VX failed: $(tail -1 "$TMP/s5w2s/$m.log")"
    echo "   $m: window [$LO, $HI] (W=0), [$LO1, $HI1] (W=1); broadening $BMEV meV (W=0), $B1 meV (W=1); M=256, 32 random-phase states seed 12345, disorder seed 7"
    # FFTW convention and the wannier2sparse CSR (brief.md section 3): PASS/FAIL rows are counted
    ( cd "$REPO" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR5:$REPO/python" python3 "$TMP/s5py/s5_conv.py" "$TMP/s5models/$m" $m \
        "$TMP/s5w2s/$m" "$TMP/s5w2s/${m}_6x4" ) >"$TMP/s5_conv_$m.txt" 2>&1
    sed 's/^/          /' "$TMP/s5_conv_$m.txt"
    NCMP=$((NCMP + $(grep -c '^\(PASS\|FAIL\)' "$TMP/s5_conv_$m.txt")))
    NFAIL=$((NFAIL + $(grep -c '^FAIL' "$TMP/s5_conv_$m.txt")))
    grep -q '^\(PASS\|FAIL\)' "$TMP/s5_conv_$m.txt" || { echo "DIFFERS   $m: convention checks not produced"; NFAIL=$((NFAIL + 1)); }
    for T in 16 1; do
      for w in dos kuboxx; do
        C="--model-dir $TMP/s5models/$m --seed $m --seed 12345 --supercell $S5N $S5N 1 --word $w --moments 256"
        W="$TMP/s5/T$T/$m/W0"
        run_s5 $T "$W/real/$w"   "$m W=0 real  $w (random 32)"  cli $C --basis real  --trace random --states 32 --window $LO $HI --broadening-meV $BMEV
        run_s5 $T "$W/bloch/$w"  "$m W=0 bloch $w (random 32)"  cli $C --basis bloch --trace random --states 32 --window $LO $HI --broadening-meV $BMEV
        run_s5 $T "$W/blochx/$w" "$m W=0 bloch $w (exact)"      cli $C --basis bloch --trace exact               --window $LO $HI --broadening-meV $BMEV
        run_s5 $T "$W/eigen/$w"  "$m W=0 eigen $w (exact)"      cli $C --basis eigen --trace exact               --window $LO $HI --broadening-meV $BMEV
        W="$TMP/s5/T$T/$m/W1"
        run_s5 $T "$W/real/$w"   "$m W=1 real  $w (random 32)"  dis $C --basis real  --trace random --states 32 --disorder-w 1 --window $LO1 $HI1 --broadening-meV $B1
        run_s5 $T "$W/bloch/$w"  "$m W=1 bloch $w (random 32)"  dis $C --basis bloch --trace random --states 32 --disorder-w 1 --window $LO1 $HI1 --broadening-meV $B1
      done
    done
    # thread counts: the real path must be byte-identical between T = 1 and 16 (counted);
    # the Bloch and eigen paths are reported (FFTW / block reduction order may differ)
    for W in "$TMP/s5/T16/$m"/W*/*/*; do
      rel="${W#$TMP/s5/T16/}"
      for f in "$W"/*.chebmom "$W"/*.dat; do
        [ -f "$f" ] || continue
        g="$TMP/s5/T1/$rel/$(basename "$f")"
        case "$rel" in
          */real/*) compare "$g" "$f" "[T=1 vs T=16] $rel/$(basename "$f")" ;;
          *) if cmp -s "$g" "$f"; then echo "INFO      [T=1 vs T=16] $rel/$(basename "$f") identical"
             else echo "INFO      [T=1 vs T=16] $rel/$(basename "$f") differs: $(python3 "$REPO/scripts/golden_diff.py" "$g" "$f" 2>&1)"; fi ;;
        esac
      done
    done
  done
  # the tolerance table of codex.md section 2, at T = 16 and T = 1
  for T in 16 1; do
    echo "   tolerance table, T=$T"
    ( cd "$REPO" && PYTHONPATH="$PYDIR5:$REPO/python" python3 "$REPO/figures/stage_5_table.py" "$TMP/s5/T$T" "$GOLD/bases" ) >"$TMP/s5_table_T$T.txt" 2>&1
    sed 's/^/          /' "$TMP/s5_table_T$T.txt"
    NCMP=$((NCMP + $(grep -c '^\(PASS\|FAIL\)' "$TMP/s5_table_T$T.txt")))
    NFAIL=$((NFAIL + $(grep -c '^FAIL' "$TMP/s5_table_T$T.txt")))
    grep -q '^\(PASS\|FAIL\)' "$TMP/s5_table_T$T.txt" || { echo "DIFFERS   [T=$T] stage-5 table not produced"; NFAIL=$((NFAIL + 1)); }
  done
  [ -z "${S5_KEEP:-}" ] || { rm -rf "$S5_KEEP"; cp -r "$TMP/s5" "$S5_KEEP"; echo "   stage-5 outputs kept in $S5_KEEP"; }
}

# FFTW belongs to the algebra layer (BlochTransform) only; the Chebyshev engine and io never see it.
stage5_grep() {
  local hits
  hits="$(cd "$REPO" && grep -rni "fftw" lsquant/core/chebyshev lsquant/io)"
  if [ -z "$hits" ]; then
    echo "GREP      grep -rni \"fftw\" lsquant/core/chebyshev lsquant/io: empty"
  else
    echo "GREP      grep -rni \"fftw\" lsquant/core/chebyshev lsquant/io: FOUND"
    echo "$hits" | sed 's/^/          /'
    NFAIL=$((NFAIL + 1))
  fi
}

case "$STAGE" in
  0) stage0 ;;
  1) stage1; core_grep ;;
  2) stage1; core_grep; stage2; stage2_grep ;;
  3) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep ;;
  4) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage4; stage4_grep ;;
  5) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage4; stage4_grep; stage5; stage5_grep ;;
  5only) stage5; stage5_grep ;;
  *) echo "ERROR: stage $STAGE has no gate defined yet"; exit 2 ;;
esac

echo "== $NCMP files compared; failures: $NFAIL"
[ "$NFAIL" -eq 0 ] && { echo "GATE PASS (stage $STAGE)"; exit 0; }
echo "GATE FAIL (stage $STAGE)"; exit 1
