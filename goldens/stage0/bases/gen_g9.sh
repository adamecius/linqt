#!/bin/bash
# Generate the G9 goldens (stage 5) from the independent oracle wannier2sparse/tools/hr_exactdiag.py.
# usage: gen_g9.sh <repo>   (writes <repo>/goldens/stage0/bases/<model>/)
set -euo pipefail
REPO="$1"; W2S=~/wannier2sparse; NK=48; M=256; PAD=0.10
OUT="$REPO/goldens/stage0/bases"
WORK=$(mktemp -d /tmp/g9.XXXXXX); cd "$WORK"
python3 "$W2S/examples/gen_models.py" >/dev/null
for m in graphene haldane; do
  cd "$WORK/$m"
  "$W2S/build/wannier2sparse" "$m" 1 1 1 VX VY --mode bundle -o bund --quiet >/dev/null
  cp bund/$m.w2sp/operators/VX.hr.dat ${m}_vx_hr.dat
  cp bund/$m.w2sp/operators/VY.hr.dat ${m}_vy_hr.dat
  # window from the eigenvalues of H(k) on the commensurate N x N Gamma mesh (k = n/N), pad 10 %
  read LO HI A SIG NG <<< "$(python3 - "$m" $NK $PAD $M <<'PY'
import sys, math, numpy as np
sys.path.insert(0, __import__('os').path.expanduser('~/wannier2sparse/tools'))
from hr_exactdiag import read_hr, opk
m, N, pad, M = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
iR, H = read_hr(f"{m}_hr.dat")
n = np.arange(N) / N
K = np.stack(np.meshgrid(n, n, indexing='ij'), -1).reshape(-1, 2)
kf = np.c_[K, np.zeros(len(K))]
w = np.linalg.eigvalsh(0.5 * (opk(iR, H, kf) + np.conj(np.transpose(opk(iR, H, kf), (0, 2, 1)))))
emin, emax = w.min(), w.max(); c = 0.5 * (emin + emax); a = 0.5 * (emax - emin) * (1 + pad)
sig = math.pi * a / M            # Jackson resolution of M moments on the window [c-a, c+a]
ng = int(round(2 * a / 0.0025)) + 1
print("%.10f %.10f %.10f %.10f %d" % (c - a, c + a, a, sig, ng))
PY
)"
  D="$OUT/$m"; mkdir -p "$D/model"
  cp ${m}_hr.dat $m.uc $m.xyz ${m}_vx_hr.dat ${m}_vy_hr.dat "$D/model/"
  python3 "$W2S/tools/hr_exactdiag.py" dos   $m --nk $NK --eta $SIG --ngrid $NG --emin $LO --emax $HI --out "$D/G9_dos"
  python3 "$W2S/tools/hr_exactdiag.py" sigma $m --nk $NK --eta $SIG --ngrid $NG --emin $LO --emax $HI --comp xx --out "$D/G9_sigma_xx"
  printf "window_lower %s\nwindow_upper %s\nhalf_width %s\nmoments %s\nsigma_eV %s\nngrid %s\nnk %s\npad %s\n" \
    $LO $HI $A $M $SIG $NG $NK $PAD > "$D/G9_params.txt"
  echo "$m: window [$LO, $HI] a=$A sigma=$SIG eV ngrid=$NG"
done
rm -rf "$WORK"
