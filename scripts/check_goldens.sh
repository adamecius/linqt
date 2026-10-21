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
#            (fft_grid/); G4c (graphene VX-VY, window from goldens/stage0/graphene/
#            offcentre_enclosing/BOUNDS = -8.5 7.5, exact trace, 256 states, M = 128)
#            moments (legacy), and its Bastin curve at 10 meV from the golden and from
#            the new moments vs the curve blessed from the new code (legacy defect D2) (a DIFFERS curve also gets an informational value/energy column
#            note); grep -rn "compat/" lsquant/core lsquant/io must be empty.
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
#   stage 5  Bloch and eigen bases (handoff 0005, codex.md revised after run 1): stages 1 to
#            4 in full, then python3 -m lsquant.apps.stage5_bases on graphene and Haldane
#            (goldens/stage0/bases/<model>/model, unchanged wannier2sparse gen_models.py files),
#            supercell 48 48 1, M = 256, --window eigen, 40 meV, real basis = the wannier2sparse
#            CSR (--csr-dir), at OMP_NUM_THREADS = 16 and 1: W = 0 real and Bloch (random phase,
#            32 states, seed 12345), Bloch and eigen (exact trace), DOS and kuboxx; W = 1 eV
#            (--disorder 1 --disorder-seed 7) real and Bloch; FFTW convention and the
#            wannier2sparse CSR rows (counted); real-path files of T = 1 and 16 byte-identical
#            (Bloch and eigen reported); figures/stage_5_table.py: PASS/FAIL rows counted
#            (real vs Bloch moments 1e-9, Bloch vs eigen curves 1e-8, EigenBasis exact vs G9
#            1e-8), INFO rows reported (KPM vs G9, finite size); grep -rni fftw
#            lsquant/core/chebyshev lsquant/io must be empty. Needs W2S_BIN and W2S_TOOLS
#            (default ~/wannier2sparse). S5_KEEP=<dir> keeps the outputs for figures/stage_5.py.
#            Stage "5only" runs the stage-5 part alone.
#   stage 7  words, run files, CLI and MCP (handoff 0007): for every run file under
#            runs/goldens/, `python3 -m lsquant run <file> --out <tmp>` at OMP_NUM_THREADS =
#            16 and 1 (S7_THREADS), then cmp of every produced file against its golden: G1, G2,
#            G3, G4, G4c, G5, G6, G7 (goldens/stage0), tutorials 01 to 05 (every non-PNG file
#            of goldens/stage0/tutorials/0N, same relative path); G9 run files against
#            stage5_bases run with the run file's arguments (byte) plus figures/stage_5_table.py
#            on the run-file outputs; G10 run files against stage6_nonorth with the run file's
#            arguments (byte) plus max relative difference vs G1/G2/G5 <= 1e-10 (stage-6 brief
#            section 6, exact applier). Then the MCP server (g1.json, g5.json through the `run`
#            and `status` tools, outputs cmp'd; a model dir outside --root must be refused), the
#            CLI refusals (non-zero exit, stderr, no output dir), the greps of stages 1 to 5 and
#            the stage-7 greps. Totals per golden and thread count at the end. S7_KEEP=<dir>
#            keeps the run-file outputs (figures/stage_7.py). Stage "7aux" runs only the MCP
#            and refusal parts. S7_SET=quick runs only the closure subset of handoff 0007
#            (g1 g2 g3 g4 g4c g5 g7_s{x,y,z} g9_graphene_{real,eigen} g10_haldane_real_generalised
#            tutorial01 tutorial03 tutorial04) at 16 threads; S7_SET=all (default) runs every run file.
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
G4CXY=NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact

# Informational, not a verdict: for a two-column curve that DIFFERS, say whether the value
# column is byte-identical as text and what the energy-column offset (golden - new) is.
columns_note() {   # <new> <golden>
  [ -f "$1" ] && [ -f "$2" ] || return 0
  cmp -s "$1" "$2" && return 0
  local v; cmp -s <(awk '{print $2}' "$1") <(awk '{print $2}' "$2") && v="identical" || v="DIFFERS"
  echo "          note: value column $v as text; energy column golden - new: $(paste "$1" "$2" | awk '{d=$3-$1; s+=d; if(NR==1||d<mn)mn=d; if(NR==1||d>mx)mx=d} END{printf "min %.6g max %.6g mean %.6g over %d lines", mn, mx, s/NR, NR}')"
}

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
    # 5. G4c: graphene VX-VY, off-centre window from the BOUNDS file (-8.5 7.5), exact trace
    #    256 states, M = 128: moments vs the legacy moments; Bastin at 10 meV, from the golden
    #    and from the produced moments, vs the curve blessed from the new code (legacy defect
    #    D2, goldens/stage0/legacy_defects/offcentre_axis/).
    for G in offcentre_enclosing; do
      W="$TMP/s3_T${T}_g4c_$G"
      run_s3 $T "$W" "G4c $G VX-VY moments" --operators "$GOLD/graphene/operators" --label graphene_golden --left VX --right VY \
        --window $(cat "$GOLD/graphene/$G/BOUNDS") --trace exact --states 256 --moments 128
      compare "$W/$G4CXY.chebmom2D" "$GOLD/graphene/$G/$G4CXY.chebmom2D" "[T=$T] graphene/$G/$G4CXY.chebmom2D (G4c moments)"
      run_s3 $T "$W/cg" "G4c $G bastin (golden moments)" --moments "$GOLD/graphene/$G/$G4CXY.chebmom2D" --formula bastin --broadening-meV 10
      run_s3 $T "$W/cn" "G4c $G bastin (new moments)" --moments "$W/$G4CXY.chebmom2D" --formula bastin --broadening-meV 10
      for c in cg cn; do
        compare "$W/$c/KuboBastin_${G4CXY}JACKSON.dat" "$GOLD/graphene/$G/KuboBastin_${G4CXY}JACKSON.dat" \
                "[T=$T] graphene/$G/KuboBastin_${G4CXY}JACKSON.dat (G4c curve, 10 meV, from $([ $c = cg ] && echo golden || echo new) moments)"
        columns_note "$W/$c/KuboBastin_${G4CXY}JACKSON.dat" "$GOLD/graphene/$G/KuboBastin_${G4CXY}JACKSON.dat"
      done
    done
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
# goldens/stage0/bases/<model>/model/), supercell 48 x 48 x 1, M = 256, window
# `--window eigen` (the app's default), broadening 40 meV (= G9's Gaussian eta).
S5N=48
S5_W2S="${W2S_BIN:-$HOME/wannier2sparse/build/wannier2sparse}"
S5_BMEV=40

# Python helper of the stage-5 gate, written into $TMP (it never touches the repository).
s5_helpers() {
  mkdir -p "$TMP/s5py"
  # FFTW convention of BlochTransform against the dense B (6 x 4 x 1, non-square on purpose),
  # and the real basis against the wannier2sparse supercell CSR at 6 x 4 x 1 and 48 x 48 x 1
  # (brief section 3: the real basis IS the wannier2sparse CSR, B VX_w2s B^+ = v_x to 1e-10).
  # Prints PASS/FAIL rows (counted by the gate) and INFO rows. The helper name starts with s5_
  # so that it never shadows a standard-library module.
  cat >"$TMP/s5py/s5_conv.py" <<'PY'
import sys, numpy as np
from lsquant import _core
from lsquant.models.wannier import load_wannier_model
from lsquant.bases import make_basis
mdir, seed, w2s_dir, w2s_small = sys.argv[1:5]
m = load_wannier_model(mdir, seed=seed)
uc, tau, W = m.uc, m.tau, m.num_wann
def gate(ok, text): print(("PASS  " if ok else "FAIL  ") + text)
rng = np.random.default_rng(1)
def rv(n): return rng.standard_normal(n) + 1j * rng.standard_normal(n)
def rel(a, b): return np.max(np.abs(a - b)) / np.max(np.abs(b))
# dense B on 6 x 4 x 1, real index (R, alpha) -> iR * W + alpha with iR = n1 + N1 n2 + N1 N2 n3
N = (6, 4, 1); Nk = 24
k = m.k_mesh(N); ph = m.bloch_phases(k)
Bf = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.RealToBloch)
Bb = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.BlochToReal)
n = np.array([(n1, n2, n3) for n3 in range(N[2]) for n2 in range(N[1]) for n1 in range(N[0])], float)
Rc = n @ uc
for sgn in (-1, +1):
    Bd = np.zeros((Nk * W, Nk * W), complex)
    for al in range(W):
        Bd[al::W, al::W] = np.exp(sgn * 1j * (k @ (Rc + tau[al]).T)) / np.sqrt(Nk)
    xs = [rv(Nk * W) for _ in range(4)]
    e_f = max(rel(Bf.multiply(x), Bd @ x) for x in xs)
    e_b = max(rel(Bb.multiply(x), Bd.conj().T @ x) for x in xs)
    line = ("CONV %s 6x4x1 dense B = exp(%si k.(R+tau))/sqrt(Nk), iR = n1 + N1 n2: RealToBloch vs B x %.3e, BlochToReal vs B^+ x %.3e"
            % (seed, '+' if sgn > 0 else '-', e_f, e_b))
    if sgn < 0:
        gate(max(e_f, e_b) <= 1e-12, line + " (brief.md section 3, tol 1e-12)")
    else:
        print("INFO  " + line + " (opposite sign, reported)")
def csr_checks(N, d, Bf=None, Bb=None):
    dim = N[0] * N[1] * N[2] * W
    Hr = _core.SparseMatrixOperator(_core.ReadCsrText(d + "/%s.HAM.CSR" % seed))
    Vr = _core.SparseMatrixOperator(_core.ReadCsrText(d + "/%s.VX.CSR" % seed))
    s = "%dx%dx%d" % N
    xs = [rv(dim) for _ in range(4)]
    # the real basis as stage5_bases builds it: with the wannier2sparse CSRs (--csr-dir) and without
    for lab, cd in (("loaded from --csr-dir", d), ("built by lsquant", None)):
        rb = make_basis("real", m, N, csr_dir=cd) if cd else make_basis("real", m, N)   # model dir holds no CSR
        Hp = rb.hamiltonian(); Vp = rb.operator("VX")
        eh = max(rel(Hp.multiply(x), Hr.multiply(x)) for x in xs); ev = max(rel(Vp.multiply(x), Vr.multiply(x)) for x in xs)
        gate(max(eh, ev) <= 1e-12, "CONV %s %s RealBasis (%s) vs wannier2sparse CSR: H %.3e, VX %.3e (tol 1e-12)" % (seed, s, lab, eh, ev))
    Hp = _core.SparseMatrixOperator(m.hamiltonian_csr(N)); Vp = _core.SparseMatrixOperator(m.velocity_csr(N, "VX"))
    eh = max(rel(Hp.multiply(x), Hr.multiply(x)) for x in xs); ev = max(rel(Vp.multiply(x), Vr.multiply(x)) for x in xs)
    gate(max(eh, ev) <= 1e-12, "CONV %s %s lsquant supercell builder (WannierModel.hamiltonian_csr, velocity_csr) vs wannier2sparse CSR: H %.3e, VX %.3e (tol 1e-12)" % (seed, s, eh, ev))
    if Bf is not None:
        Hk = _core.BlockDiagonalOperator(m.Hk(m.k_mesh(N))); Vk = _core.BlockDiagonalOperator(m.dHk(m.k_mesh(N), "VX"))
        e1 = max(rel(Bf.multiply(Hr.multiply(Bb.multiply(x))), Hk.multiply(x)) for x in xs)
        e2 = max(rel(Bf.multiply(Vr.multiply(Bb.multiply(x))), Vk.multiply(x)) for x in xs)
        gate(e1 <= 1e-10, "CONV %s %s B H_w2s B^+ x vs H_k x: %.3e (tol 1e-10)" % (seed, s, e1))
        gate(e2 <= 1e-10, "CONV %s %s B VX_w2s B^+ x vs v_x(k) x: %.3e (tol 1e-10, brief.md section 3)" % (seed, s, e2))
csr_checks((6, 4, 1), w2s_small)
N = (48, 48, 1); Nk = 48 * 48
k = m.k_mesh(N); ph = m.bloch_phases(k)
Bf = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.RealToBloch)
Bb = _core.BlochTransform(*N, W, ph.ravel(), _core.BlochDirection.BlochToReal)
xs = [rv(Nk * W) for _ in range(4)]
e = max(rel(Bb.multiply(Bf.multiply(x)), x) for x in xs)
gate(e <= 1e-12, "CONV %s 48x48x1 max|B^+ B x - x|/max|x| = %.3e (tol 1e-12)" % (seed, e))
csr_checks(N, w2s_dir, Bf, Bb)
PY
}

# One stage5_bases run in its own directory (wall time printed)
run_s5() {   # <threads> <workdir> <tag> <args...>
  local T="$1" W="$2" tag="$3" rc s e; shift 3
  mkdir -p "$W"; s=$(date +%s.%N)
  ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR5:$REPO/python" python3 -m lsquant.apps.stage5_bases \
      "$@" --out "$W" ) >"$W/run.log" 2>"$W/err.log"; rc=$?
  e=$(date +%s.%N)
  [ $rc -eq 0 ] || { echo "          [T=$T] $tag exited $rc: $(tail -1 "$W/err.log")"; NFAIL=$((NFAIL + 1)); }
  printf "          [T=%s] %-34s wall %7.2f s  -> %s\n" "$T" "$tag" "$(echo "$e - $s" | bc)" "$(ls "$W" | grep -v '\.log$' | tr '\n' ' ')"
}

stage5() {
  local T m w W C MD CSR rel f g
  PYDIR5="$BUILD/lsquant_core"
  ls "$PYDIR5"/_core*.so >/dev/null 2>&1 || PYDIR5="$BUILD/lsquant"
  ls "$PYDIR5"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  [ -x "$S5_W2S" ] || { echo "ERROR: wannier2sparse binary $S5_W2S missing (set W2S_BIN)"; exit 2; }
  echo "== stage 5: python3 -m lsquant.apps.stage5_bases, PYTHONPATH=$PYDIR5:python, vs goldens/stage0/bases (G9)"
  s5_helpers
  for m in graphene haldane; do
    MD="$GOLD/bases/$m/model"; CSR="$TMP/s5w2s/$m"
    mkdir -p "$CSR" "$TMP/s5w2s/${m}_6x4"
    ( cd "$MD" && "$S5_W2S" $m $S5N $S5N 1 VX -o "$CSR" --quiet --no-log-file \
      && "$S5_W2S" $m 6 4 1 VX -o "$TMP/s5w2s/${m}_6x4" --quiet --no-log-file ) >"$TMP/s5w2s/$m.log" 2>&1 \
      || echo "          wannier2sparse $m VX failed: $(tail -1 "$TMP/s5w2s/$m.log")"
    echo "   $m: model $MD (unchanged gen_models.py files), real basis = wannier2sparse CSR ($m $S5N $S5N 1 VX) via --csr-dir;"
    echo "      --window eigen (default), --broadening-meV $S5_BMEV, M=256, 32 random-phase states seed 12345, W=1: --disorder 1 --disorder-seed 7"
    ( cd "$REPO" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR5:$REPO/python" python3 "$TMP/s5py/s5_conv.py" "$MD" $m \
        "$CSR" "$TMP/s5w2s/${m}_6x4" ) >"$TMP/s5_conv_$m.txt" 2>&1
    sed 's/^/          /' "$TMP/s5_conv_$m.txt"
    NCMP=$((NCMP + $(grep -c '^\(PASS\|FAIL\)' "$TMP/s5_conv_$m.txt")))
    NFAIL=$((NFAIL + $(grep -c '^FAIL' "$TMP/s5_conv_$m.txt")))
    grep -q '^\(PASS\|FAIL\)' "$TMP/s5_conv_$m.txt" || { echo "DIFFERS   $m: convention checks not produced"; NFAIL=$((NFAIL + 1)); }
    for T in 16 1; do
      for w in dos kuboxx; do
        C="--model $MD --seed $m --seed 12345 --supercell $S5N $S5N 1 --word $w --moments 256 --broadening-meV $S5_BMEV --csr-dir $CSR"
        W="$TMP/s5/T$T/$m/W0"
        run_s5 $T "$W/real/$w"   "$m W=0 real  $w (random 32)" $C --basis real  --trace random --states 32
        run_s5 $T "$W/bloch/$w"  "$m W=0 bloch $w (random 32)" $C --basis bloch --trace random --states 32
        run_s5 $T "$W/blochx/$w" "$m W=0 bloch $w (exact)"     $C --basis bloch --trace exact
        run_s5 $T "$W/eigen/$w"  "$m W=0 eigen $w (exact)"     $C --basis eigen --trace exact
        W="$TMP/s5/T$T/$m/W1"
        run_s5 $T "$W/real/$w"   "$m W=1 real  $w (random 32)" $C --basis real  --trace random --states 32 --disorder 1 --disorder-seed 7
        run_s5 $T "$W/bloch/$w"  "$m W=1 bloch $w (random 32)" $C --basis bloch --trace random --states 32 --disorder 1 --disorder-seed 7
      done
    done
    # thread counts: the real path must be byte-identical between T = 1 and 16 (counted);
    # the Bloch and eigen paths are reported
    for W in "$TMP/s5/T16/$m"/W*/*/*; do
      rel="${W#$TMP/s5/T16/}"
      for f in "$W"/*.chebmom "$W"/*.chebmom1D "$W"/*.chebmom2D "$W"/*.dat "$W"/*.basis.json; do
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
    ( cd "$REPO" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR5:$REPO/python" python3 "$REPO/figures/stage_5_table.py" "$TMP/s5/T$T" "$GOLD/bases" ) >"$TMP/s5_table_T$T.txt" 2>&1
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

# ---------------------------------------------------------------- stage 7 (words, run files, CLI, MCP)
# Every golden from words: `python -m lsquant run runs/goldens/<name>.json --out <tmp>` and a
# cmp of each produced file with its golden. The old per-stage sections (stage1..stage5) are
# not run at stage 7; the mapping from their compared files to the run files is:
#   stage 1 G1, G2 (stage1_moments)          -> g1.json, g2.json
#   stage 2 G3 (+ G1 via --operators)        -> g3.json (G3 curve, DOS_reconstructed.dat, G1)
#   stage 3 G5 moments + 4 curves            -> g5.json;  G4 moments + 5 curves -> g4.json;
#           G4c moments + Bastin curve       -> g4c.json
#   stage 4 G6 moments + 4 curves            -> g6.json;  G7 3 moments -> g7.json (and g7_s{x,y,z}.json)
#   stage 5 G9 (stage5_bases, no byte golden)-> g9_<model>_<cfg>.json, cmp'd against stage5_bases run
#           with the run file's own arguments, and figures/stage_5_table.py on the run-file outputs
#   stage 6 G10 (stage6_nonorth, exact)      -> g10_<model>_<variant>_<route>.json, cmp'd against
#           stage6_nonorth with the run file's arguments, and vs G1/G2/G5 at 1e-10 (stage-6 brief s. 6)
#   tutorials 01..05 (stage-0 captures)      -> tutorial0N.json, every non-PNG file under
#           goldens/stage0/tutorials/0N at the same relative path
# Not compared from the old sections: stage-1 batch-16 and random-phase determinism, stage-2
# form (a), the stage-3/4 "curves from the golden moment file" rows (no run file reconstructs
# from a golden file; `lsquant reconstruct` is not part of the gate), the stage-5 convention
# rows (s5_conv.py). Totals are printed per golden and per thread count.
S7_THREADS="${S7_THREADS:-16 1}"
S7_RUNS="$REPO/runs/goldens"
G1F2=NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
declare -A S7I S7D S7ORDER
# S7_SET=all (default): every run file. S7_SET=quick: the closure subset of handoff 0007 (16 threads
# by default; g6, g7, the other g9 and g10 files, tutorial02 and tutorial05 deferred to the final run).
S7_SET="${S7_SET:-all}"
S7_QUICK="g1 g2 g3 g4 g4c g5 g7_sx g7_sy g7_sz g9_graphene_real g9_graphene_eigen g10_haldane_real_generalised tutorial01 tutorial03 tutorial04"
[ "$S7_SET" = quick ] && S7_THREADS="${S7_THREADS_QUICK:-16}"
s7want() { [ "$S7_SET" = all ] && return 0; case " $S7_QUICK " in *" $1 "*) return 0 ;; esac; return 1; }

# captured/index.tsv: every line "<NNNN_name> <sha1> <source>" must name a produced capture with that sha1
s7idx() {   # <group> <produced captured dir> <golden index.tsv> <tag>
  local g="$1" d="$2" idx="$3" tag="$4" bad=0 n=0 name sha
  [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
  NCMP=$((NCMP + 1))
  while IFS=$'\t' read -r name sha _; do
    n=$((n + 1))
    [ -f "$d/$name" ] && [ "$(sha1sum "$d/$name" | cut -d' ' -f1)" = "$sha" ] || { bad=$((bad + 1)); echo "          index.tsv line $n: $name sha1 $sha not matched by $d/$name"; }
  done <"$idx"
  if [ $bad -eq 0 ] && [ $n -gt 0 ]; then echo "IDENTICAL $tag ($n entries, every sha1 matches)"; S7I[$g]=$(( ${S7I[$g]:-0} + 1 ))
  else echo "DIFFERS   $tag ($bad of $n entries do not match)"; NFAIL=$((NFAIL + 1)); S7D[$g]=$(( ${S7D[$g]:-0} + 1 )); fi
}
S7GROUPS=()

s7cmp() {   # <group> <new> <golden> <tag>   (compare(), plus a per-group tally)
  local g="$1" before=$NFAIL
  [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
  compare "$2" "$3" "$4"
  if [ "$NFAIL" -eq "$before" ]; then S7I[$g]=$(( ${S7I[$g]:-0} + 1 )); else S7D[$g]=$(( ${S7D[$g]:-0} + 1 )); fi
}

# max relative difference (scripts/golden_diff.py definition), 0 when cmp-identical
s7relval() { cmp -s "$1" "$2" && { echo 0; return; }; python3 "$REPO/scripts/golden_diff.py" "$1" "$2" 2>&1 | grep -o 'max rel diff [^ ]*' | awk '{print $4}'; }

s7rel() {   # <group> <new> <golden> <tol|INFO> <tag>
  local g="$1" new="$2" gold="$3" tol="$4" tag="$5" r
  [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
  if [ ! -f "$new" ]; then
    [ "$tol" = INFO ] && { echo "INFO      $tag  (not produced)"; return; }
    echo "FAIL      $tag  (not produced)"; NFAIL=$((NFAIL + 1)); NCMP=$((NCMP + 1)); S7D[$g]=$(( ${S7D[$g]:-0} + 1 )); return
  fi
  r="$(s7relval "$new" "$gold")"; [ -n "$r" ] || r=inf
  if [ "$tol" = INFO ]; then echo "INFO      $tag  max rel diff $r"; return; fi
  NCMP=$((NCMP + 1))
  if awk -v r="$r" -v t="$tol" 'BEGIN{exit !(r <= t)}'; then
    echo "PASS      $tag  max rel diff $r <= $tol"; S7I[$g]=$(( ${S7I[$g]:-0} + 1 ))
  else
    echo "FAIL      $tag  max rel diff $r > $tol"; NFAIL=$((NFAIL + 1)); S7D[$g]=$(( ${S7D[$g]:-0} + 1 ))
  fi
}

run_s7() {   # <threads> <run file name> <outdir>
  local T="$1" n="$2" W="$3" rc s e
  rm -rf "$W"; mkdir -p "$(dirname "$W")"; s=$(date +%s.%N)
  ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR7:$REPO/python" python3 -m lsquant run "runs/goldens/$n.json" --out "$W" ) \
    >"$W.out" 2>"$W.err"; rc=$?
  e=$(date +%s.%N)
  [ $rc -eq 0 ] || { echo "          [T=$T] lsquant run runs/goldens/$n.json exited $rc: $(tail -1 "$W.err")"; NFAIL=$((NFAIL + 1)); }
  printf "          [T=%s] python -m lsquant run runs/goldens/%-36s wall %8.2f s, %4d files written\n" "$T" "$n.json" "$(echo "$e - $s" | bc)" "$(find "$W" -type f 2>/dev/null | wc -l)"
}

s7_ref5() {   # <threads> <model> <cfg> <outdir>: stage5_bases with the run file's arguments
  local T="$1" m="$2" c="$3" O="$4" w a
  case "$c" in
    real) a="--basis real --trace random --states 32" ;;
    bloch) a="--basis bloch --trace random --states 32" ;;
    bloch_exact) a="--basis bloch --trace exact" ;;
    eigen) a="--basis eigen --trace exact" ;;
    real_w1) a="--basis real --trace random --states 32 --disorder 1 --disorder-seed 7" ;;
    bloch_w1) a="--basis bloch --trace random --states 32 --disorder 1 --disorder-seed 7" ;;
  esac
  for w in dos kuboxx; do
    mkdir -p "$O/$w"
    ( cd "$REPO" && OMP_NUM_THREADS=$T PYTHONPATH="$PYDIR7:$REPO/python" python3 -m lsquant.apps.stage5_bases \
        --model "$GOLD/bases/$m/model" --seed $m --seed 12345 --supercell 48 48 1 --word $w --moments 256 \
        --broadening-meV 40 $a --out "$O/$w" ) >"$O/$w.log" 2>&1 || echo "          [T=$T] stage5_bases $m $c $w failed: $(tail -1 "$O/$w.log")"
  done
}

s7_ref6() {   # <threads> <model> <variant> <route> <outdir>: stage6_nonorth with the run file's arguments
  local T="$1" m="$2" v="$3" r="$4" O="$5" w win st M
  if [ "$m" = chain1d ]; then win="-2 2"; st=2048; M=64; set -- dos kubo; else win="-4.5 4.5"; st=512; M=128; set -- kuboxy kubo; fi
  mkdir -p "$O"
  for w in "$@"; do
    OMP_NUM_THREADS=$T "$EXE6" --fixture "$GOLD/nonorth/$m/$v" --label $m --window $win --route $r --applier exact \
      --order 60 --tol 1e-13 --trace exact --states $st --moments $M --word $w --out "$O" >>"$O.log" 2>&1 \
      || echo "          [T=$T] stage6_nonorth $m $v $r $w failed: $(tail -1 "$O.log")"
  done
}

s7_goldens() {   # <threads>: every run file, every golden
  local T="$1" B="$TMP/s7/T$1" n f g m c v r W R
  for f in "$S7_RUNS"/*.json; do n="$(basename "$f" .json)"; s7want $n && run_s7 $T $n "$B/$n"; done
  echo "   [T=$T] per-golden comparison"
  s7want g1 && s7cmp "T$T G1" "$B/g1/$G1F" "$GOLD/chain1d/$G1F" "[T=$T] g1.json -> chain1d/$G1F (G1)"
  s7want g2 && s7cmp "T$T G2" "$B/g2/$G1F2" "$GOLD/chain1d/$G1F2" "[T=$T] g2.json -> chain1d/$G1F2 (G2)"
  for g in $G3F DOS_reconstructed.dat $G1F; do
    s7want g3 && s7cmp "T$T G3" "$B/g3/$g" "$GOLD/chain1d/$g" "[T=$T] g3.json -> chain1d/$g (G3)"
  done
  if s7want g4; then
    s7cmp "T$T G4" "$B/g4/$G4XX.chebmom2D" "$GOLD/graphene/$G4XX.chebmom2D" "[T=$T] g4.json -> graphene/$G4XX.chebmom2D (G4)"
    s7cmp "T$T G4" "$B/g4/KuboBastin_${G4XX}JACKSON.dat" "$GOLD/graphene/KuboBastin_${G4XX}JACKSON.dat" "[T=$T] g4.json -> graphene/KuboBastin_${G4XX}JACKSON.dat (G4, 64 meV)"
    for g in sea_surf/KuboBastinSea_${G4XX}JACKSON.dat sea_surf/KuboBastinSurf_${G4XX}JACKSON.dat \
             fft_grid/KuboGreenWood_${G4XX}JACKSON_FFTgrid.dat fft_grid/KuboBastin_${G4XX}JACKSON_FFTgrid.dat; do
      s7cmp "T$T G4" "$B/g4/$g" "$GOLD/graphene/$g" "[T=$T] g4.json -> graphene/$g (G4, 10 meV)"
    done
  fi
  if s7want g4c; then
    s7cmp "T$T G4c" "$B/g4c/$G4CXY.chebmom2D" "$GOLD/graphene/offcentre_enclosing/$G4CXY.chebmom2D" "[T=$T] g4c.json -> graphene/offcentre_enclosing/$G4CXY.chebmom2D (G4c)"
    s7cmp "T$T G4c" "$B/g4c/KuboBastin_${G4CXY}JACKSON.dat" "$GOLD/graphene/offcentre_enclosing/KuboBastin_${G4CXY}JACKSON.dat" "[T=$T] g4c.json -> graphene/offcentre_enclosing/KuboBastin_${G4CXY}JACKSON.dat (G4c)"
  fi
  for g in $G5XY.chebmom2D $G5XX.chebmom2D KuboBastin_${G5XY}JACKSON.dat KuboBastinSea_${G5XY}JACKSON.dat \
           KuboBastinSurf_${G5XY}JACKSON.dat KuboGreenWood_${G5XX}JACKSON.dat; do
    s7want g5 && s7cmp "T$T G5" "$B/g5/$g" "$GOLD/haldane/$g" "[T=$T] g5.json -> haldane/$g (G5)"
  done
  if s7want g6; then
    s7cmp "T$T G6" "$B/g6/$G6M.chebmomTD" "$GOLD/chain1d_ballistic/$G6M.chebmomTD" "[T=$T] g6.json -> chain1d_ballistic/$G6M.chebmomTD (G6)"
    for g in $G6EF; do
      s7cmp "T$T G6" "$B/g6/mean${G6M}EF${g}JACKSON.dat" "$GOLD/chain1d_ballistic/mean${G6M}EF${g}JACKSON.dat" "[T=$T] g6.json -> chain1d_ballistic/mean${G6M}EF${g}JACKSON.dat (G6)"
    done
  fi
  for g in SX SY SZ; do
    f=TimeProj$g-PXchain1d_magKPM_M32_stateexact.chebmomTD
    s7want g7 && s7cmp "T$T G7" "$B/g7/$f" "$GOLD/chain1d_mag/$f" "[T=$T] g7.json -> chain1d_mag/$f (G7)"
    n=g7_$(echo $g | tr 'A-Z' 'a-z')
    s7want $n && s7cmp "T$T G7" "$B/$n/$f" "$GOLD/chain1d_mag/$f" "[T=$T] $n.json -> chain1d_mag/$f (G7)"
  done
  # tutorials: every non-PNG file of the stage-0 tutorial golden, at the same relative path; the
  # capture index captured/index.tsv of the stage-0 shim (name, sha1, source path in the legacy
  # pipeline's temporary directory) is metadata, checked by sha1 against the produced captures
  for n in 01 02 03 04 05; do
    s7want tutorial$n || continue
    W="$B/tutorial$n"
    while IFS= read -r f; do
      case "$f" in
        captured/index.tsv) s7idx "T$T Tut$n" "$W/captured" "$GOLD/tutorials/$n/$f" "[T=$T] tutorial$n.json -> tutorials/$n/$f (Tut $n, sha1 column vs the produced captures)" ;;
        *) s7cmp "T$T Tut$n" "$W/$f" "$GOLD/tutorials/$n/$f" "[T=$T] tutorial$n.json -> tutorials/$n/$f (Tut $n)" ;;
      esac
    done < <(cd "$GOLD/tutorials/$n" && find . -type f ! -name '*.png' | sed 's|^\./||' | sort)
    while IFS= read -r f; do
      [ -f "$GOLD/tutorials/$n/$f" ] || echo "INFO      [T=$T] tutorial$n.json wrote $f (no golden at that path)"
    done < <(cd "$W" 2>/dev/null && find . -type f | sed 's|^\./||' | sort)
  done
  # G9: run files vs stage5_bases with the same arguments (byte), then the stage-5 table on the run files
  for m in graphene haldane; do
    for c in real bloch bloch_exact eigen real_w1 bloch_w1; do
      s7want g9_${m}_$c || continue
      R="$TMP/s7ref5/T$T/$m/$c"; s7_ref5 $T $m $c "$R"
      for f in "$R"/dos/* "$R"/kuboxx/*; do
        [ -f "$f" ] || continue
        s7cmp "T$T G9" "$B/g9_${m}_$c/$(basename "$f")" "$f" "[T=$T] g9_${m}_$c.json -> stage5_bases $c $(basename "$(dirname "$f")")/$(basename "$f") (G9)"
      done
    done
    for c in real:W0/real bloch:W0/bloch bloch_exact:W0/blochx eigen:W0/eigen real_w1:W1/real bloch_w1:W1/bloch; do
      R="$TMP/s7tab/T$T/$m/${c#*:}"; mkdir -p "$R/dos" "$R/kuboxx"
      for f in "$B/g9_${m}_${c%%:*}"/*; do
        case "$(basename "$f")" in
          *SpectralOp1*) ln -sf "$f" "$R/dos/" ;;
          *NonEqOp*) ln -sf "$f" "$R/kuboxx/" ;;
        esac
      done
    done
  done
  if [ "$S7_SET" = all ]; then
  echo "   [T=$T] figures/stage_5_table.py on the run-file outputs"
  ( cd "$REPO" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR7:$REPO/python" python3 "$REPO/figures/stage_5_table.py" "$TMP/s7tab/T$T" "$GOLD/bases" ) >"$TMP/s7_table_T$T.txt" 2>&1
  sed 's/^/          /' "$TMP/s7_table_T$T.txt"
  g="T$T G9 table"; [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
  S7I[$g]=$(grep -c '^PASS' "$TMP/s7_table_T$T.txt"); S7D[$g]=$(grep -c '^FAIL' "$TMP/s7_table_T$T.txt")
  NCMP=$((NCMP + S7I[$g] + S7D[$g])); NFAIL=$((NFAIL + S7D[$g]))
  [ $((S7I[$g] + S7D[$g])) -gt 0 ] || { echo "DIFFERS   [T=$T] stage-5 table not produced"; NFAIL=$((NFAIL + 1)); }
  else
    echo "INFO      [T=$T] figures/stage_5_table.py skipped (S7_SET=$S7_SET: the table needs all six G9 configurations of both models)"
  fi
  # G10: run files vs stage6_nonorth with the same arguments (byte), then vs the orthogonal goldens
  for m in chain1d haldane; do
    for v in real complex; do
      for r in lowdin generalised; do
        s7want g10_${m}_${v}_$r || continue
        # stage6_nonorth (dense exact S^-1/2 on 2048 sites: ~45 min per word on a loaded TCN) is run
        # once, at the first thread count, and only when the run file produced its moment files
        n=g10_${m}_${v}_$r; R="$TMP/s7ref6/$n"
        if ! ls "$B/$n"/*.chebmom* >/dev/null 2>&1; then
          g="T$T G10"; [ -n "${S7ORDER[$g]:-}" ] || { S7ORDER[$g]=1; S7GROUPS+=("$g"); }
          echo "DIFFERS   [T=$T] $n.json produced no moment file ($(tail -1 "$B/$n.err")); stage6_nonorth reference not run"
          NFAIL=$((NFAIL + 1)); NCMP=$((NCMP + 1)); S7D[$g]=$(( ${S7D[$g]:-0} + 1 ))
        else
          [ -d "$R" ] || s7_ref6 $T $m $v $r "$R"
        fi
        for f in "$R"/*; do
          [ -f "$f" ] && ls "$B/$n"/*.chebmom* >/dev/null 2>&1 || continue
          s7cmp "T$T G10" "$B/$n/$(basename "$f")" "$f" "[T=$T] $n.json -> stage6_nonorth $(basename "$f") (G10)"
        done
        if [ $m = chain1d ]; then
          s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G1F" "$GOLD/chain1d/$G1F" 1e-10 "[T=$T] $n.json $G1F vs G1"
          s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G1F2" "$GOLD/chain1d/$G1F2" 1e-10 "[T=$T] $n.json $G1F2 vs G2"
          s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G3F" "$GOLD/chain1d/$G3F" INFO "[T=$T] $n.json $G3F vs G3 (curve, 6 digits)"
        else
          s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G5XY.chebmom2D" "$GOLD/haldane/$G5XY.chebmom2D" 1e-10 "[T=$T] $n.json $G5XY.chebmom2D vs G5"
          s7rel "T$T G10 vs G1/G2/G5" "$B/$n/$G5XX.chebmom2D" "$GOLD/haldane/$G5XX.chebmom2D" 1e-10 "[T=$T] $n.json $G5XX.chebmom2D vs G5"
          for g in KuboBastin_${G5XY} KuboBastinSea_${G5XY} KuboBastinSurf_${G5XY} KuboGreenWood_${G5XX}; do
            s7rel "T$T G10 vs G1/G2/G5" "$B/$n/${g}JACKSON.dat" "$GOLD/haldane/${g}JACKSON.dat" INFO "[T=$T] $n.json ${g}JACKSON.dat vs G5 (curve, 6 digits)"
          done
        fi
      done
    done
  done
}

# the MCP server: run g1.json and g5.json through `run`, poll `status`, cmp the outputs
s7_mcp() {
  local R="$TMP/s7mcp" g
  mkdir -p "$R/goldens/stage0/chain1d" "$R/goldens/stage0/haldane"
  cp -r "$GOLD/chain1d/operators" "$R/goldens/stage0/chain1d/"
  cp -r "$GOLD/haldane/operators" "$R/goldens/stage0/haldane/"
  cat >"$TMP/s7_mcp_client.py" <<'PY'
import asyncio, json, os, sys, time
def show(prefix, obj):
    print(prefix + " " + json.dumps(obj, separators=(", ", ": "))); sys.stdout.flush()
async def main(root, runfiles):
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client
    params = StdioServerParameters(command=sys.executable, args=["-m", "lsquant.mcp_server", "--root", root], env=dict(os.environ))
    async with stdio_client(params) as (rd, wr):
        async with ClientSession(rd, wr) as s:
            await s.initialize()
            show("MCP <- list_tools:", [t.name for t in (await s.list_tools()).tools])
            async def call(name, args):
                show("MCP -> call_tool:", {"tool": name, "arguments": args})
                t0 = time.time(); res = await s.call_tool(name, args); dt = time.time() - t0
                text = res.content[0].text if res.content else ""
                try: body = json.loads(text)
                except ValueError: body = text
                err = bool(getattr(res, "isError", False) or getattr(res, "is_error", False))
                show("MCP <- %s (isError=%s, returned in %.3f s):" % (name, err, dt), body)
                return body, err
            for rfp in runfiles:
                name = os.path.splitext(os.path.basename(rfp))[0]
                raw = json.load(open(rfp, encoding="utf-8"))
                raw.setdefault("output", {})["dir"] = "out/" + name
                body, err = await call("run", {"runfile_json": json.dumps(raw)})
                if err or not isinstance(body, dict):
                    print("MCP RESULT %s: run refused" % name); continue
                st, _ = await call("status", {"run_id": body["run_id"]}); n = 1
                while st["status"] == "running":
                    await asyncio.sleep(2.0)
                    st = json.loads((await s.call_tool("status", {"run_id": body["run_id"]})).content[0].text); n += 1
                show("MCP <- status poll %d (final):" % n, st)
                print("MCP RESULT %s: %s" % (name, st["status"]))
            raw = json.load(open(runfiles[0], encoding="utf-8"))
            raw["model"]["dir"] = "/etc"; raw["output"]["dir"] = "out/outside"
            body, err = await call("run", {"runfile_json": json.dumps(raw)})
            refused = err or (isinstance(body, str) and "error" in body.lower())
            print("MCP RESULT confinement (model.dir=/etc): %s (isError=%s)" % ("refused" if refused else "NOT refused", err))
            print("MCP RESULT out dir created: %s" % os.path.exists(os.path.join(root, "out", "outside")))
asyncio.run(main(os.path.realpath(sys.argv[1]), sys.argv[2:]))
PY
  echo "== stage 7 MCP: python -m lsquant.mcp_server --root $R (copies of the chain1d and haldane golden operators), OMP_NUM_THREADS=16"
  ( cd "$R" && OMP_NUM_THREADS=16 PYTHONPATH="$PYDIR7:$REPO/python" python3 "$TMP/s7_mcp_client.py" "$R" "$S7_RUNS/g1.json" "$S7_RUNS/g5.json" ) >"$TMP/s7_mcp.txt" 2>"$TMP/s7_mcp.err"
  echo "          (server stderr: $(wc -l <"$TMP/s7_mcp.err") lines; first error line: $(grep -m1 -E "Error|ValueError" "$TMP/s7_mcp.err" | sed "s/  */ /g"))"
  sed "s|$TMP|\$TMP|g; s/^/          /" "$TMP/s7_mcp.txt"
  { grep -q "MCP RESULT confinement (model.dir=/etc): refused" "$TMP/s7_mcp.txt" && grep -q "MCP RESULT out dir created: False" "$TMP/s7_mcp.txt"; } || { echo "FAIL      MCP confinement"; NFAIL=$((NFAIL + 1)); }
  NCMP=$((NCMP + 1))
  s7cmp "MCP" "$R/out/g1/$G1F" "$GOLD/chain1d/$G1F" "[MCP] g1.json -> chain1d/$G1F (G1)"
  for g in $G5XY.chebmom2D $G5XX.chebmom2D KuboBastin_${G5XY}JACKSON.dat KuboBastinSea_${G5XY}JACKSON.dat \
           KuboBastinSurf_${G5XY}JACKSON.dat KuboGreenWood_${G5XX}JACKSON.dat; do
    s7cmp "MCP" "$R/out/g5/$g" "$GOLD/haldane/$g" "[MCP] g5.json -> haldane/$g (G5)"
  done
}

# refusals through the CLI: non-zero exit, message on stderr, nothing written
s7_refuse() {   # <tag> <args...>
  local tag="$1" rc O="$TMP/s7ref/out_$1"; shift
  ( cd "$REPO" && PYTHONPATH="$PYDIR7:$REPO/python" python3 -m lsquant "$@" ) >"$TMP/s7ref/$tag.out" 2>"$TMP/s7ref/$tag.err"; rc=$?
  NCMP=$((NCMP + 1))
  if [ $rc -ne 0 ] && [ -s "$TMP/s7ref/$tag.err" ] && [ ! -e "$O" ]; then
    echo "REFUSED   [$tag] exit $rc, stderr: $(cat "$TMP/s7ref/$tag.err")"
  else
    echo "FAIL      [$tag] exit $rc, stdout: $(head -c 300 "$TMP/s7ref/$tag.out" | tr '\n' ' ') stderr: $(cat "$TMP/s7ref/$tag.err"); output dir exists: $([ -e "$O" ] && echo yes || echo no)"
    NFAIL=$((NFAIL + 1))
  fi
}

s7_refusals() {
  local R="$TMP/s7ref" f
  mkdir -p "$R"
  echo "== stage 7 refusals through the CLI (python -m lsquant run <bad>.json --out <dir>; the dir must not appear)"
  # a run file built from g1.json / g2.json / g6.json with one change; jq-free, by python
  mk() {   # <name> <base> <python expression on d>
    python3 - "$S7_RUNS/$2.json" "$R/$1.json" "$3" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
for stmt in sys.argv[3].split(";"):
    k, v = stmt.split("=", 1); keys = k.strip().split("."); t = d
    for kk in keys[:-1]: t = t[kk]
    if v.strip() == "DEL": t.pop(keys[-1], None)
    else: t[keys[-1]] = json.loads(v)
json.dump(d, open(sys.argv[2], "w"), indent=1)
PY
  }
  mk unknown_operator g1 'word="Tr[FOO delta(E-H)]"'
  mk rank_3 g2 'model.operators=["VX","VY","VZ"]; word="Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]"'
  mk mixed_basis g1 'model.operators=["VX"]; word="Tr[VX@bloch delta(E-H)]"'
  mk ut_without_time g6 'time=DEL; reconstruction=DEL; word="Tr[VX delta(E-H) U(t)^+ VX U(t)]"'
  mk unknown_key_top g1 'unknown_extra_key=1; model.dir="does/not/exist"'
  mk unknown_key_nested g1 'expansion.windw=[-2,2]; model.dir="does/not/exist"'
  for f in unknown_operator rank_3 mixed_basis ut_without_time unknown_key_top unknown_key_nested; do
    s7_refuse $f run "$R/$f.json" --out "$R/out_$f"
  done
  # the same words through `words describe` (no run file): reported
  for f in "Tr[FOO delta(E-H)]" "Tr[VX delta(E-H) VY delta(E-H) VZ delta(E-H)]" "Tr[VX@bloch delta(E-H)]" \
           "Tr[VX delta(E-H) U(t)^+ VX U(t)]" "Tr[VX delta(E-H) VY G+(E)]" "Tr[VX delta(E-H) VY dG+(E)/dE]" \
           "Tr[VX delta(E-H) VY G-(E)]" "Tr[VX G+(E) VY delta(E-H)]"; do
    ( cd "$REPO" && PYTHONPATH="$PYDIR7:$REPO/python" python3 -m lsquant words describe "$f" ) >"$R/d.out" 2>"$R/d.err"
    echo "INFO      words describe \"$f\": exit $?, $( [ -s "$R/d.err" ] && cat "$R/d.err" || python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); c=d.get("catalogue",{}); print("ACCEPTED as", c.get("shape_id"), "formulas", c.get("default_formulas"), "engine word", repr(c.get("engine_word")))' "$R/d.out")"
  done
}

stage7_setup() {
  PYDIR7="$BUILD/lsquant_core"
  ls "$PYDIR7"/_core*.so >/dev/null 2>&1 || PYDIR7="$BUILD/lsquant"
  ls "$PYDIR7"/_core*.so >/dev/null 2>&1 || { echo "ERROR: _core*.so not found under $BUILD/lsquant_core or $BUILD/lsquant"; exit 2; }
  EXE6="$BUILD/lsquant_core/stage6_nonorth"
  [ -x "$EXE6" ] || { echo "ERROR: $EXE6 missing"; exit 2; }
}

s7_totals() {
  local g
  echo "== stage 7 totals per golden (IDENTICAL or PASS / DIFFERS or FAIL)"
  for g in "${S7GROUPS[@]}"; do printf "   %-26s %4d identical  %4d differ\n" "$g" "${S7I[$g]:-0}" "${S7D[$g]:-0}"; done
}

stage7() {
  local T
  stage7_setup
  echo "== stage 7: python3 -m lsquant run runs/goldens/*.json ($(ls "$S7_RUNS"/*.json | wc -l) run files, S7_SET=$S7_SET), PYTHONPATH=$PYDIR7:python, threads: $S7_THREADS"
  [ "$S7_SET" = all ] || echo "   deferred (not run in this set): $(for f in "$S7_RUNS"/*.json; do n=$(basename "$f" .json); s7want $n || printf "%s " $n; done)"
  for T in $S7_THREADS; do s7_goldens $T; done
  s7_mcp
  s7_refusals
  s7_totals
  [ -z "${S7_KEEP:-}" ] || { rm -rf "$S7_KEEP"; cp -r "$TMP/s7" "$S7_KEEP"; echo "   stage-7 outputs kept in $S7_KEEP"; }
}

# words and run files live in python/lsquant; the engine knows none of it, and nothing evaluates text
stage7_grep() {
  local all hits coms pat='catalogue|runfile|grammar|Tr\[|mcp'
  local com='^[^:]+:[0-9]+:[[:space:]]*(//|\*|/\*)'
  all="$(cd "$REPO" && grep -rnE "$pat" lsquant/core lsquant/io)"
  hits="$(printf '%s\n' "$all" | grep -Ev "$com" | grep -v '^$')"
  coms="$(printf '%s\n' "$all" | grep -E "$com")"
  if [ -z "$hits" ]; then echo "GREP      grep -rnE \"$pat\" lsquant/core lsquant/io, code lines: empty"
  else echo "GREP      grep -rnE \"$pat\" lsquant/core lsquant/io, code lines: FOUND"; echo "$hits" | sed 's/^/          /'; NFAIL=$((NFAIL + 1)); fi
  [ -z "$coms" ] || echo "$coms" | sed 's/^/          excluded comment line: /'
  hits="$(cd "$REPO" && grep -rnE --include='*.py' "\b(eval|exec)\(|literal_eval|exec_module" python/lsquant/words python/lsquant/runfile.py python/lsquant/cli.py python/lsquant/mcp_server.py)"
  if [ -z "$hits" ]; then echo "GREP      grep -rnE \"eval(|exec(|literal_eval|exec_module\" python/lsquant/{words,runfile.py,cli.py,mcp_server.py}: empty"
  else echo "GREP      eval/exec in the word language or run file: FOUND"; echo "$hits" | sed 's/^/          /'; NFAIL=$((NFAIL + 1)); fi
  hits="$(cd "$REPO" && grep -rnE --include='*.py' "\bexec_module\b" python/lsquant)"
  [ -z "$hits" ] || { echo "INFO      exec_module elsewhere in python/lsquant (tutorial generators imported by path):"; echo "$hits" | sed 's/^/          /'; }
}
case "$STAGE" in
  0) stage0 ;;
  1) stage1; core_grep ;;
  2) stage1; core_grep; stage2; stage2_grep ;;
  3) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep ;;
  4) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage4; stage4_grep ;;
  5) stage1; core_grep; stage2; stage2_grep; stage3; stage3_grep; stage4; stage4_grep; stage5; stage5_grep ;;
  5only) stage5; stage5_grep ;;
  7) stage7; core_grep; stage2_grep; stage3_grep; stage4_grep; stage5_grep; stage7_grep ;;
  7aux) stage7_setup; s7_mcp; s7_refusals; s7_totals ;;
  *) echo "ERROR: stage $STAGE has no gate defined yet"; exit 2 ;;
esac

echo "== $NCMP files compared; failures: $NFAIL"
[ "$NFAIL" -eq 0 ] && { echo "GATE PASS (stage $STAGE)"; exit 0; }
echo "GATE FAIL (stage $STAGE)"; exit 1
