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
#   stage 6  non-orthogonal metric (handoff 0006): stages 1, 2 and 3 in full, then
#            build_dir/lsquant_core/stage6_nonorth on the rotated fixtures
#            goldens/stage0/nonorth/{chain1d,haldane}/{real,complex}: routes lowdin and
#            generalised (naive as evidence only), appliers exact, ischo, chebyshev,
#            moments vs G1, G2 (chain1d dos, kubo) and G5 (Haldane kuboxy = VX-VY,
#            kubo = VX-VX) at the brief's tolerances (not bytes; see the stage-6 block),
#            the two routes against each other, and the Haldane Hall plateau C from the
#            stage-3 Bastin reconstruction of every VX-VY file (|C - 1| <= 0.05).
#            S6_KEEP=<dir> keeps the runs; S6_THREADS (default 16) sets OMP_NUM_THREADS.
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

# ---------------------------------------------------------------- stage 6 (non-orthogonal metric)
# Rotated models (goldens/stage0/nonorth/<model>/<variant>, fixtures of handoff 0006) must
# reproduce the orthogonal goldens G1, G2 (chain1d) and the G5 moments (Haldane VX-VY,
# VX-VX). These runs cannot be byte-identical (the applier is a factorisation, an
# iteration or a truncated series), so the criterion is the brief's tolerance (section 6):
#   exact 1e-10; ischo (--tol 1e-12) 1e-9; chebyshev (ctor truncation coeffTol 1e-12) 1e-8.
# Measure (tester's decision, see tests.md):
#   rel = max_k |mu_new[k] - mu_gold[k]| / max_k |mu_gold[k]|  over all complex moments;
# the per-element relative difference of scripts/golden_diff.py is printed as well, but it
# is meaningless on moments that vanish by symmetry (e.g. the odd DOS moments).
# The two header lines (N BandWidth BandCenter / M) must be numerically identical.
# Rows: PASS / FAIL (counted); route naive = the legacy symmetric estimator, EVIDENCE (not
# counted); lowdin + ischo = N/A (ISCHO provides S^-1 only; stage6_nonorth rejects it).
# Then lowdin vs generalised with the same applier (generalised/ischo vs lowdin/exact) at
# the same levels, and the Hall plateau C = (mean sigma_xy over |E| < 0.5) / A * 2 pi,
# A = (dim/2) 3 sqrt(3)/2 (test/hall_haldane.sh, tutorial03.chern_from_plateau), from the
# stage-3 Bastin reconstruction at 10 meV of every produced VX-VY file: PASS if |C-1| <= 0.05.
# With S6_KEEP=<dir> the run directories, gate_table.tsv and chern.tsv are kept there
# (figures/stage_6.py reads them). S6_THREADS sets OMP_NUM_THREADS (default 16).
S6CMP=""
S6_LAST=""
s6_tools() {
  S6CMP="$TMP/s6cmp.py"
  cat > "$S6CMP" <<'PYEOF'
import sys
def load(p):
    L = open(p).read().split("\n")
    h = [L[0].split(), L[1].split()]
    v = []
    for l in L[2:]:
        t = l.split()
        if t:
            v.append(complex(float(t[0]), float(t[1]) if len(t) > 1 else 0.0))
    return h, v
a, b = load(sys.argv[1]), load(sys.argv[2])
hdr = [[float(x) for x in r] for r in a[0]] == [[float(x) for x in r] for r in b[0]]
if len(a[1]) != len(b[1]):
    print("nan header_ok=%s length %d vs %d" % (hdr, len(a[1]), len(b[1])))
    sys.exit()
g = max(abs(x) for x in b[1])
d = [abs(x - y) for x, y in zip(a[1], b[1])]
k = max(range(len(d)), key=lambda i: d[i])
print("%.3e header_ok=%s max_abs=%.3e at_element=%d max|gold|=%.3e" % (d[k] / g, hdr, d[k], k, g))
PYEOF
  S6CHERN="$TMP/s6chern.py"
  cat > "$S6CHERN" <<'PYEOF'
import sys, math
import numpy as np
d = np.loadtxt(sys.argv[1])
dim = int(sys.argv[2])
if not np.all(np.isfinite(d)):
    print("nan nonfinite")
    sys.exit()
E, s = d[:, 0], d[:, 1]
m = np.abs(E) < 0.5
A = (dim / 2.0) * 3.0 * math.sqrt(3.0) / 2.0
p = float(np.mean(s[m]))
print("%.6f plateau=%.6f points=%d" % (p / A * 2 * math.pi, p, int(m.sum())))
PYEOF
}

# s6_judge <new> <gold> <tol> <tag> <gate|evidence>; sets S6_LAST="<rel> <verdict>"
s6_judge() {
  local new="$1" gold="$2" tol="$3" tag="$4" kind="$5" out rel hdr verdict
  if [ ! -f "$new" ] || [ ! -f "$gold" ]; then
    verdict=FAIL; out="nan (missing: $([ -f "$new" ] || basename "$new") $([ -f "$gold" ] || echo "${gold#$TMP/}"))"
  else
    out="$(python3 "$S6CMP" "$new" "$gold" 2>&1)"
    rel="${out%% *}"; hdr="$(echo "$out" | grep -o 'header_ok=[A-Za-z]*')"
    if python3 -c "import sys; r=float('$rel'); sys.exit(0 if (r == r and r <= $tol and '$hdr' == 'header_ok=True') else 1)" 2>/dev/null
    then verdict=PASS; else verdict=FAIL; fi
  fi
  rel="${out%% *}"
  if [ "$kind" = evidence ]; then
    echo "EVIDENCE  $tag  rel=$out  (would be $verdict at tol $tol; not counted)"
  else
    NCMP=$((NCMP + 1))
    [ "$verdict" = PASS ] || NFAIL=$((NFAIL + 1))
    printf '%-9s %s  rel=%s  tol=%s\n' "$verdict" "$tag" "$out" "$tol"
  fi
  if [ -f "$new" ] && [ -f "$gold" ]; then
    echo "          golden_diff.py: $(python3 "$REPO/scripts/golden_diff.py" "$new" "$gold" 2>&1)"
  fi
  S6_LAST="$rel $verdict"
}

s6_field() { echo "$1" | grep -o " $2=[^ ]*" | head -1 | cut -d= -f2; }

stage6() {
  local EXE="$BUILD/lsquant_core/stage6_nonorth" PYDIR6 model variant word route applier W tol
  local args lo hi M gold tag line rc res it eo f b P C kind
  local XY=NonEqOpVX-VYhaldaneKPM_M128x128_stateexact
  [ -x "$EXE" ] || EXE="$BUILD/lsquant/stage6_nonorth"
  [ -x "$EXE" ] || { echo "ERROR: stage6_nonorth not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  PYDIR6="$(dirname "$EXE")"
  # Invoked through a symlink in $TMP (same binary): on TCN another session runs a loop
  # `pkill -f "stage6_nonorth.*--route lowdin"` that would kill the Lowdin runs (tests.md).
  local RUNEXE="$TMP/s6_nonorth_exe"; ln -sf "$EXE" "$RUNEXE"
  s6_tools
  local KEEP="${S6_KEEP:-}" TSV="$TMP/s6_gate_table.tsv" CTSV="$TMP/s6_chern.tsv" T6="${S6_THREADS:-16}"
  printf 'model\tvariant\tword\troute\tapplier\tmax_rel\ttol\tverdict\tresidual\titerations\teffective_order\n' > "$TSV"
  echo "== stage 6: $EXE on goldens/stage0/nonorth vs G1, G2, G5 moments (OMP_NUM_THREADS=$T6)"
  for model in chain1d haldane; do
    if [ $model = chain1d ]; then lo=-2; hi=2; M=64; else lo=-4.5; hi=4.5; M=128; fi
    for variant in real complex; do
      for word in $([ $model = chain1d ] && echo "dos kubo" || echo "kuboxy kubo"); do
        case "$model/$word" in
          chain1d/dos)    gold="$GOLD/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" ;;
          chain1d/kubo)   gold="$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" ;;
          haldane/kuboxy) gold="$GOLD/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D" ;;
          haldane/kubo)   gold="$GOLD/haldane/NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D" ;;
        esac
        b="$(basename "$gold")"; P="$TMP/s6/${model}_${variant}_${word}"
        for route in lowdin generalised naive; do
          for applier in exact ischo chebyshev; do
            [ $route = naive ] && [ $applier != exact ] && continue
            tag="$model/$variant $word $route/$applier"
            if [ $route = lowdin ] && [ $applier = ischo ]; then
              echo "N/A       $tag  (ISCHO provides S^-1 only; the Lowdin route needs S^-1/2)"
              printf '%s\t%s\t%s\t%s\t%s\tN/A\t-\tN/A\t-\t-\t-\n' $model $variant $word $route $applier >> "$TSV"
              continue
            fi
            case $applier in
              exact) tol=1e-10; args="" ;;
              ischo) tol=1e-9; args="--tol 1e-12" ;;
              chebyshev) tol=1e-8; args="" ;;
            esac
            W="${P}_${route}_${applier}"; mkdir -p "$W"
            ( cd "$REPO" && OMP_NUM_THREADS=$T6 "$RUNEXE" --fixture "$GOLD/nonorth/$model/$variant" \
                --label $model --window $lo $hi --route $route --applier $applier $args --trace exact \
                --moments $M --word $word --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
            [ $rc -eq 0 ] || echo "          [$tag] stage6_nonorth exited $rc: $(tail -1 "$W/err.log")"
            line="$(grep '^stage6_nonorth:' "$W/run.log" | tail -1)"
            echo "          stdout: ${line//$W\//}"
            kind=gate; [ $route = naive ] && kind=evidence
            s6_judge "$W/$b" "$gold" $tol "$tag" $kind
            res="$(s6_field "$line" residual)"; it="$(s6_field "$line" iterations)"; eo="$(s6_field "$line" effective_order)"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $model $variant $word $route $applier \
              "${S6_LAST% *}" $tol "$([ $kind = evidence ] && echo "EVIDENCE(${S6_LAST#* })" || echo "${S6_LAST#* }")" \
              "${res:--}" "${it:--}" "${eo:--}" >> "$TSV"
          done
        done
        # the two routes against each other
        s6_judge "${P}_generalised_exact/$b" "${P}_lowdin_exact/$b" 1e-10 \
          "$model/$variant $word generalised/exact vs lowdin/exact" gate
        printf '%s\t%s\t%s\tgeneralised_vs_lowdin\texact\t%s\t1e-10\t%s\t-\t-\t-\n' $model $variant $word ${S6_LAST} >> "$TSV"
        s6_judge "${P}_generalised_chebyshev/$b" "${P}_lowdin_chebyshev/$b" 1e-8 \
          "$model/$variant $word generalised/chebyshev vs lowdin/chebyshev" gate
        printf '%s\t%s\t%s\tgeneralised_vs_lowdin\tchebyshev\t%s\t1e-8\t%s\t-\t-\t-\n' $model $variant $word ${S6_LAST} >> "$TSV"
        s6_judge "${P}_generalised_ischo/$b" "${P}_lowdin_exact/$b" 1e-9 \
          "$model/$variant $word generalised/ischo vs lowdin/exact" gate
        printf '%s\t%s\t%s\tgeneralised_ischo_vs_lowdin_exact\tischo\t%s\t1e-9\t%s\t-\t-\t-\n' $model $variant $word ${S6_LAST} >> "$TSV"
      done
    done
  done
  # Hall plateau from every produced Haldane VX-VY file (stage-3 Bastin, 10 meV)
  echo "== stage 6: Hall plateau, python3 -m lsquant.apps.stage3_kubo --formula bastin --broadening-meV 10"
  printf 'source\tC\tverdict\n' > "$CTSV"
  C="$(python3 "$S6CHERN" "$GOLD/haldane/KuboBastin_${XY}JACKSON.dat" 512 2>&1)"
  echo "REFERENCE Hall golden G5 KuboBastin (legacy)  C=$C"
  printf 'golden_G5\t%s\tREF\n' "${C%% *}" >> "$CTSV"
  for W in "$TMP"/s6/haldane_*_kuboxy_*; do
    tag="$(basename "$W")"
    if [ -f "$W/$XY.chebmom2D" ]; then
      ( cd "$REPO" && OMP_NUM_THREADS=$T6 PYTHONPATH="$PYDIR6:$REPO/python" python3 -m lsquant.apps.stage3_kubo \
          --moments "$W/$XY.chebmom2D" --formula bastin --broadening-meV 10 --out "$W" ) >"$W/bastin.log" 2>&1 \
        || echo "          [$tag] stage3_kubo exited non-zero: $(tail -1 "$W/bastin.log")"
      C="$(python3 "$S6CHERN" "$W/KuboBastin_${XY}JACKSON.dat" 512 2>&1 | tail -1)"
    else
      C="nan (no VX-VY moments produced)"
    fi
    if python3 -c "import sys; c=float('${C%% *}'); sys.exit(0 if abs(c - 1) <= 0.05 else 1)" 2>/dev/null
    then rc=PASS; else rc=FAIL; fi
    case "$tag" in
      *_naive_*) echo "EVIDENCE  Hall $tag  C=$C  (would be $rc; not counted)"; rc="EVIDENCE($rc)" ;;
      *) NCMP=$((NCMP + 1)); [ $rc = PASS ] || NFAIL=$((NFAIL + 1))
         printf '%-9s Hall %s  C=%s  (|C - 1| <= 0.05)\n' $rc "$tag" "$C" ;;
    esac
    printf '%s\t%s\t%s\n' "$tag" "${C%% *}" "$rc" >> "$CTSV"
  done
  if [ -n "$KEEP" ]; then
    mkdir -p "$KEEP"; rm -rf "$KEEP/runs"; cp -r "$TMP/s6" "$KEEP/runs"
    cp "$TSV" "$KEEP/gate_table.tsv"; cp "$CTSV" "$KEEP/chern.tsv"
    echo "          kept: $KEEP/{runs,gate_table.tsv,chern.tsv}"
  fi
}

case "$STAGE" in
  0) stage0 ;;
  1) stage1; core_grep ;;
  2) stage1; core_grep; stage2; stage2_grep ;;
  3) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep ;;
  6) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage6 ;;
  *) echo "ERROR: stage $STAGE has no gate defined yet"; exit 2 ;;
esac

echo "== $NCMP files compared; failures: $NFAIL"
[ "$NFAIL" -eq 0 ] && { echo "GATE PASS (stage $STAGE)"; exit 0; }
echo "GATE FAIL (stage $STAGE)"; exit 1
