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
#            at OMP_NUM_THREADS = 1, 8, 64; random-phase determinism (two runs, cmp);
#            grep for "basis|std::cout|getenv" under lsquant/core (must be empty).
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

case "$STAGE" in
  0) stage0 ;;
  1) stage1; core_grep ;;
  *) echo "ERROR: stage $STAGE has no gate defined yet"; exit 2 ;;
esac

echo "== $NCMP files compared; failures: $NFAIL"
[ "$NFAIL" -eq 0 ] && { echo "GATE PASS (stage $STAGE)"; exit 0; }
echo "GATE FAIL (stage $STAGE)"; exit 1
