# Stage 6b implementation report — 2026-10-22

Developer: Codex, standing in for Agi. Sources were edited in
`/root/lsquant-codex`, whose branch remains `handoff/0006b-dev`, initially
`b57b810` (the `ff1ad43` code plus handoff documents). No push, build, golden
gate, material calculation, or figure was run. The supplied environment has
no network, no compiled extension, and no SciPy, SymPy, or sisl. The user
explicitly assigned builds and substantial numerical runs to the tester on
Aron; that instruction overrides the usual pre-commit golden requirement.

The brief acquired sections 11 and 12 during implementation. Both were read
and incorporated: the SIESTA path has no dependency on Wannier inputs, and
`metric_forms.py` supplies the three symbolic metric representations. Those
concurrent edits to the brief, and the concurrent `.codex_prompt.txt`, were
not authored, staged, or committed by this developer.

## Commit transport limitation

The workspace mounts `.git` read-only. The first ordinary commit attempt
failed before staging with:

```text
fatal: Unable to create '/root/lsquant-codex/.git/index.lock': Read-only file system
```

Consequently the actual branch reference cannot be advanced here. Coherent
commits are recorded using detached temporary Git metadata at
`/tmp/lsquant-stage6b.git`, based on the original branch HEAD, with the same
working tree and read-only access to the original object store. No second
branch was created. A portable bundle is supplied at
`/tmp/lsquant-stage6b.bundle`; it contains the new commit chain and depends on
`b57b810`. On a writable checkout already on `handoff/0006b-dev` at that base,
fetch the bundle's `HEAD` and fast-forward `FETCH_HEAD`. Concurrent handoff
document edits should remain with their author. Nothing was pushed.

Code commits before this report:

- `7f91ad1`: strict CSR validation, generic Bloch sums, velocity/spin formulas,
  and the corrected Wannier position consumer.
- `f5b4845`: SIESTA input/quadrature and cached Bloch metric adapters.
- `3c044c2`: symbolic metric forms, three numerical routes, and spin Hall oracle.
- `d9370f7`: app/run-file integration, G11 descriptors, and the single self-check.
- `01fb908`: reuse the cached overlap square root for FHS and validate app options early.

The report itself is a final documentation commit; the final response lists
the complete sequence.

## Files and responsibilities

| Files | Change |
|---|---|
| `python/lsquant/models/bloch.py` | Format-independent reciprocal lattice, n1-fastest mesh, and atom/cell Bloch sums |
| `python/lsquant/models/siesta.py` | sisl reader, real-space spinor blocks, orbital metadata, SIESTA-only radial dipole quadrature, midpoint conversion |
| `python/lsquant/models/nonorthogonal.py` | Single velocity ladder, exact block metric factorisation, physical/applied spin and spin-current matrices, strict Python CSR reader |
| `python/lsquant/models/wannier.py` | `_r.dat` reader with paired-image Hermitisation and matching `ndegen`; explicit `WannierModel.velocity_blocks` consumer |
| `python/lsquant/bases.py` | `BlochMetricBasis`, `EigenMetricBasis`, cached exact metric and eigenpairs, eigenphase/dual-bra plan adapter, model protocol without a Wannier import |
| `python/lsquant/metric_forms.py` | SymPy expressions, Hermitian conjugation, metric normal forms, symbolic trace proofs, anticommutator counterexample |
| `python/lsquant/oracle/spin_hall.py` | NumPy-only generalised eigensolve, curvature, Gaussian occupation, eigenphase trace, projected-spin FHS invariant |
| `python/lsquant/apps/stage6b_siesta.py` | Shared word executor; CLI constructs the same run file used by `lsquant run` |
| `python/lsquant/reconstruct.py` | `shc` normalisation in one place |
| `python/lsquant/runfile.py`, `words/catalogue.py`, `words/compiler.py`, `cli.py`, `__init__.py` | SIESTA run-file dispatch, three route choices, SHC reconstruction and operator names |
| `lsquant/io/csr_text.{hpp,cpp}`, `lsquant/api/bindings.cpp`, `python/lsquant/_core.pyi` | C++ CSR diagnostics/refusal and explicit `legacy_fold=False` binding |
| `runs/goldens/g11_*.json` | 24 descriptors: two materials, three velocity modes, three metric forms, plus DOS for each route in berry_connection mode |
| `tests/dev/stage_6b_check.py` | The one requested rotated Kane–Mele developer self-check |

No engine, recursion, legacy, compatibility, golden data, golden harness, or
figure code was changed. The existing C++ `BlockDiagonalOperator` already
parallelises application over k blocks, so no new C++ metric applier was
needed. Dense inverse and inverse square root construction is in NumPy,
once per basis, before the engine runs.

## SIESTA conventions and checks still requiring Aron

The implementation targets sisl 0.16.4; its conventions have **not** been
verified against a real sisl installation in this environment. The reader
includes checks against public sisl `Hk`/`Sk` in both gauges at a nonzero k,
and refuses a layout or phase mismatch. These checks, their residuals, the
sisl version, and orbital metadata are persisted in `basis.json` when run.

- The FDF reader selects the available TSHS/HSX and attaches species/basis
  metadata. The binary reader is used separately for printed Fermi metadata.
  Stored H is already relative to E_F=0; printed E_F is never subtracted.
- Spinor order is `(orbital up, orbital down)`. Real SOC storage components
  are read as `(Re Huu, Re Hdd, Re Hud, Im Hud, Im Huu, Im Hdd, Re Hdu, Im Hdu)`.
  The public-API reconstruction checks explicitly test this assumption.
- The scalar spatial overlap is lifted as `S_orb tensor I_2`. The public
  overlap may already be spin-duplicated; both public shapes are handled.
- Lattice rows and centres are in Angstrom; H is in eV; k is Cartesian in
  inverse Angstrom. Reciprocal rows are `2*pi*inv(lattice).T`.
- Atom-gauge phase is `exp(+i k.(R+tau_j-tau_i))`; Cartesian derivatives
  multiply by `+i(R+tau_j-tau_i)_a`. Bare/cell derivatives multiply by `+iR_a`
  but are dressed back into the atom gauge of H and S before multiplication.
- TSHS round-off is a separate issue from the defective CSR velocities.
  SIESTA Bloch sums measure/accept relative round-off up to `1e-7`, then
  remove it by Hermitian projection, following notes section 11. That is
  an explicit input policy; it does **not** relax the CSR `1e-10` threshold.
  Derivatives and velocities are built from the same accepted matrices.

## Dipole definition and quadrature

A necessary correction to the brief's mixed definitions is explicit in
`midpoint_connection` and in the velocity docstring. Let

```text
Q_a,ij(R) = <0i|r_a|Rj>
D_a,ij(R) = Q_a,ij(R) - tau_i,a S_ij(R)
d_ij(R)   = R + tau_j - tau_i.
```

Translation and conjugation give

```text
Q_a(-R)^+ = Q_a(R) - R_a S(R)
D_a(R) - D_a(-R)^+ = d_a(R) S(R).
```

Thus D, the row-centred integral written in brief section 3, is not a
Hermitian operator. Its Bloch anti-Hermitian part satisfies
`D(k)-D(k)^+ = -i dS_a(k)`. Equation (4), which already contains the symmetric
`dS` correction, needs its Hermitian part:

```text
A_a,ij(R) = Q_a,ij(R) - (tau_i+tau_j+R)_a S_ij(R)/2
A_a(R)   = A_a(-R)^+.
```

Using the unsymmetrised D in addition to the symmetric-dS correction counts
its anti-Hermitian contribution again. The code never guesses A from S or
centres. SIESTA covariant mode requires radial functions from `.ion.xml`;
missing functions or uncertified integrals cause an error.

For orbitals relative to their own centres, the implemented integral is

```text
S_ij(R) = integral phi_i(u)^* phi_j(u-d) d^3u
A_a,ij(R) = integral phi_i(u)^* (u-d/2)_a phi_j(u-d) d^3u.
```

Tensor Gauss–Legendre quadrature covers the intersection of the two radial
cutoff bounding boxes. `psi` vanishes outside its cutoff sphere. A plane at
a time bounds scratch storage by O(order^2). Defaults are order 24, doubling
to at most 96, with absolute refinement tolerance `1e-5` for the overlap and
`1e-5 Angstrom` for A. The same tolerance checks the quadrature overlap
against the stored SIESTA overlap. All stored H/S pairs are included,
including parity-zero S entries. Independently integrated opposite-image
pairs are compared before averaging. Maximum refinement, overlap, and
pair-Hermiticity errors are recorded in the certificate. This is an error
estimate, not a rigorous numerical integration bound. No material
quadrature certificate was generated here; convergence and cost remain
for Aron.

The fixture starts with the requested full `Q'=T^+ diag(tau) T`. It converts
Q' to A' through the same midpoint formula. In cell gauge T is constant;
in atom gauge the corresponding transformation is
`T(k)=diag(exp(-ik.tau)) T diag(exp(+ik.tau))`. The self-check verifies both
representations. It does not discard the centre phases to force a match.

## Numerical formulas and metric forms

`metric_forms.velocity_forms`, `spin_forms`, `spin_current_forms`, and
`hamiltonian_forms` are the symbolic definitions cited by the numerical
implementations. With X=S^-1, Y=S^-1/2:

```text
bare: V = dH_cell - (dS_cell X H + H X dS_cell)/2
berry_connection: V = dH - (dS X H + H X dS)/2
covariant: V = berry_connection + i(H X A - A^+ X H)
Sigma_b = S q_b, q_b = I_orb tensor sigma_b
J^b_a = (V_a q_b + q_b V_a)/2
O_c = X O; O_L = Y O Y; O_sym = (X O + O X)/2.
```

All physical matrices are Hermitian. `O_c` is generally only S-Hermitian:
`O_c^+ S = S O_c`. Spin is in Pauli units (+/-1), without hbar/2. The
spin-scalar overlap assumption is checked against all three Pauli
operators; a non-scalar overlap requires supplied physical spin matrices.

The exact metric uses `eigh(S)=U diag(s) U^+` and forms X and Y once, refusing
nonpositive s without clipping. Its residuals and the cached generalised
eigenvector residual are reported. The generalised Hamiltonian is rescaled
as `a X H + b I`, exactly `X(aH+bS)`. Lowdin is the default. The third,
`anticommutator`, route is explicitly marked approximate and uses
`H_sym=(XH+HX)/2` and the analogous applied observable matrices. Its window
comes from **its own** spectrum, never the spectrum of a different matrix.

The symbolic proof uses Z=S^1/2, YZ=ZY=I and
`H_L=Z (XH) Y`, hence `f(H_L)=Z f(XH) Y`. After substitution, expansion and
cyclic trace reduction, the generalised and Lowdin words have the same
normal form. This applies to arbitrary spectral functions, not only a
finite list of Chebyshev orders. DOS, operator DOS, Greenwood, Bastin,
TD/projected words and expanded MSD terms are covered; sums, projectors and
commutators follow by linearity.

The anticommutator has the same **single-operator** trace by cyclicity.
Products already fail in general for two static operators; the symbolic
pair residual is returned instead of asserting a nonexistent identity.
For the spectral catalogue it is generally an approximation. A sufficient
condition for exact equality of whole words is that S commute with H and
with every static operator in the word. The self-check prints the three
forms and their Hermiticity verdicts, proves the exact routes agree, and
reports the approximate plateau without requiring equality.

The existing engine's `MirrorToBra` copies the **prepared ket seed**, rather
than creating a fresh r. Therefore the adapter prepends C on the ket and S
on the copied bra, producing C r and S C r, with `C^+ S C=I`. It then executes
the unchanged stage-6 asymmetric metric plan. Lowdin uses U on the ket and
identity on the copied bra. Each sampled word is the same as its
`C^+ O C` eigenbasis word, not merely equal in ensemble mean. Random phase
vectors use the engine's `1/sqrt(dimension)` normalisation; the existing
reconstruction dimension factor restores the physical trace. The new
`EigenMetricBasis` uses direct eigenbasis phases. The historical stage-5
`EigenBasis` real-space-correlated sampling remains available unchanged.

## Oracle and units

The NumPy-only oracle evaluates

```text
H C = S C diag(E), C^+ S C=I
J_nm = (C^+ J C)_nm, V_nm = (C^+ VY C)_nm
Omega_n = -2 Im sum_(m != n) J_nm V_mn/(E_n-E_m)^2
sigma^z_xy(E_F) = 2*pi/(N_k*A_cell) sum_kn f(E_kn;E_F) Omega_kn.
```

Occupation broadening is a Gaussian standard deviation in eV,
`f=erfc((E-E_F)/(sqrt(2)*eta))/2`. Exactly degenerate multiplet-internal
terms are excluded; the occupied multiplet sum is gauge invariant when the
Fermi energy is in a gap. The curvature/occupation trace is diagonal in the
eigenbasis, so a unit-modulus eigenphase sample gives that trace exactly.
KPM rank-2 moments retain the interband sample fluctuations; the app writes
both production and eigenphase engine outputs using the same seed/count.

FHS acts on the positive and negative sectors of `C_occ^+ Sigma_z C_occ`.
Their frames are embedded in a common orthonormal representation using
`S^1/2 C_occ`; atom phases are removed **before** periodic boundary links.
Singular links, a missing insulating gap, or a closed projected-spin gap
are reported as undefined topology, never replaced by an integer.

The requested spin-Chern convention is `C_spin=C_+-C_-`, which is 2 for
Kane–Mele. The usual Prodan half-difference is also written; it is 1 there,
and `Z2=half_difference mod 2`. The brief's simultaneous requests for
`C_spin=2` and `Z2=C_spin mod 2` are inconsistent; the latter would incorrectly
assign Kane–Mele Z2=0. No absolute value is used to choose a Hall sign.
With spin mixing the physical spin Hall conductivity need not be an
integer even while the projected-spin invariant is well defined. The
oracle reports both quantities independently.

`reconstruct.shc` implements exactly
`2*pi*kubo_bastin(mu)/(N_k*A_cell)`. The stage-3 dimension factor is already
inside `kubo_bastin`; there is no extra division by the orbital count. Every
SHC curve header says `e^2/h, sigma_z Pauli`. DOS curves are per cell.

## Wannier position and CSR changes

The optional Wannier consumer implements notes section 11 before any
commutator: `A(R)=[A_raw(R)+A_raw(-R)^+]/(2*ndegen(R))`. Degeneracies are
read from the matching `_hr.dat`, since `_r.dat` has no degeneracy table.
Missing opposite images, inconsistent weights, malformed/truncated input,
and non-Hermitian image splitting are refused. Optional wsvec splitting
uses the same image weights as H. The corrected arrays are cached on the
model before evaluating `-dH_cell-i[H,A_full]`. The established stage-5
Wannier velocity sign is preserved; its simultaneous reversal in both
slots leaves a Hall word unchanged. The SIESTA reader never imports this
consumer or uses any Wannier file. Existing stage-5 default paths are
preserved; the new explicit consumer is `WannierModel.velocity_blocks`.

Both CSR readers report absolute `||(M-M^+)/2||_F` and its ratio to `||M||_F`
for **every** read, then refuse a relative value above `1e-10`, unless
`legacy_fold=True` was explicitly passed. Accepted matrices still receive
the established fold, preserving its arithmetic order. Python's older
`fold_hermitian=False` now controls the final fold only; it does not bypass
validation. Malformed CSR indices and non-finite norms are also refused.
No existing caller is silently opted into legacy folding.

## Lightweight checks actually performed

No complete developer self-check was executed. These are editing checks
and tiny formula probes only, with no golden verdict:

```text
Python 3.10 syntax: OK; G11 JSON syntax: OK (24 run files)
git diff --check: no output
Tiny 4x4 orientation probe (not convergence): sigma sign 1.0 sector difference 2.0
Single-k covariance formula probe: 1.303699802838727e-16
Single-k inverse residual: 9.995861419348e-16
```

The orientation probe loaded the NumPy-only oracle and the fixture arrays,
without the extension. The single-k probe evaluated the actual NumPy
function/class bodies extracted from the source AST because the package's
SciPy/extension dependencies are unavailable. It is not an import/integration
check, nor evidence that the full self-check or any material gate passes.
SymPy formulas were syntax checked but not executed (SymPy is absent).

## Commands and pending material runs for Aron

Build with the existing offline Aron recipe and its installed dependency
paths; no new CMake dependency is introduced. Do not fetch pybind11 from
this environment. After building, set the Python search path to the source
and the directory containing the newly built `_core` extension, for example:

```bash
export PATH="$HOME/lsq-env/bin:$PATH"
export PYTHONPATH="$PWD/python:$PWD/build/lsquant_core:$PWD/build"
KPM_SEED=12345 python tests/dev/stage_6b_check.py
python -m lsquant run runs/goldens/g11_bi_covariant_lowdin_shc.json
python -m lsquant run runs/goldens/g11_pdse2_covariant_lowdin_shc.json
```

The self-check contains all section-8 assertions: rotated cell/atom Bloch
sums (1e-12), covariance (1e-10), ladder difference, all-mode and current
Hermiticity (1e-13), Hellmann–Feynman (1e-6), dense-mesh oracle plateau
2.000000 and spin Chern 2 (1e-6), and short M=256 exact-trace KPM on 12x12
against its own finite-mesh oracle (1e-2). It also prints symbolic forms,
checks the Lowdin oracle, checks the Lowdin short KPM, reports the
anticommutator plateau, checks same-sample Bloch/eigen moment agreement,
and reports eigenphase versus exact differences. A single eigenphase DOS
sample is checked against exact moments. The script's assertions have not
been run here.

The app writes legacy-layout and rank-header moments, `shc.dat`/`dos.dat`,
`run.resolved.json`, `basis.json`, eigen reference moments/curves, and an
`eigen_comparison.json`. With oracle enabled it also writes
`exact_shc.json` and `spin_chern.json`. G11 descriptors enable the oracle,
use 48x48x1, M=1024, 16 eigenbasis phase vectors, seed 12345 and 10 meV.
The anticommutator descriptors are experiments, not blessed physical
references; their eigen reference is for the approximate operators, while
the curvature/FHS oracle still uses the physical SIESTA operators.

| Material | Velocity | Quadrature / KPM / oracle / spin Chern |
|---|---|---|
| Bi(111) | bare | Not run; numbers unavailable |
| Bi(111) | berry_connection | Not run; numbers unavailable |
| Bi(111) | covariant | Not run; numbers unavailable |
| PdSe2 | bare | Not run; numbers unavailable |
| PdSe2 | berry_connection | Not run; numbers unavailable |
| PdSe2 | covariant | Not run; numbers unavailable |

No validated material velocity mode has been selected. The descriptors
cover every mode and exact route so the tester can make that selection on
evidence. Material oracles, original-input/pseudopotential provenance,
Santiago's unit comparison, all previous golden regressions, and figures
remain the tester/reviewer's work. The coarse PdSe2 Wannier model has no
plateau claim; it remains only an optional comparison.

## Open questions for review

1. Confirm the targeted sisl APIs, FDF basis attachment, SOC component order,
   and `.ion.xml` species filenames with real Bi/PdSe2 inputs. The built-in
   public-API checks intentionally fail rather than infer another layout.
2. Run and inspect the complete self-check, including SymPy normal forms,
   engine seed handling, Hall sign and the M=256 reconstruction error.
3. Measure radial quadrature convergence, overlap agreement and total cost;
   increase the explicit quadrature order only if the certificate requires
   it. No uncertified dipole or zero-dipole fallback is allowed.
4. Review the midpoint connection derivation against the brief's raw
   row-centred definition, and the spin-Chern/Z2 convention correction.
5. Assess material spin mixing before expecting a physical SHC plateau to
   equal a projected-spin integer; report both rather than forcing equality.
6. Regression-check the strict CSR policy on previous goldens. Any rejected
   operator needs diagnosis or an explicit legacy reproduction request;
   this implementation never changes a golden or enables the flag for it.
7. Import the detached commits into the writable requested branch from the
   bundle; the read-only Git metadata prevented that final reference update
   here.

## Review 1 item 5: static cost audit (no executable profile here)

The user explicitly substitutes source inspection for profiling in this editing
sandbox: no network or build. Before: the tester reports 71 minutes at M=1024,
16 vectors, 2304 blocks of W=72, OMP_NUM_THREADS=48, 264% CPU and 95 threads.
A short M=128/R=2 cProfile before/after was **not run**; both measurements and
the after wall time/CPU remain unavailable. These are not inferred speedups.

The inspected path is `stage6b_siesta.prepare_moments` -> word compiler ->
`MomentEngine.run`, including `BlochMetricBasis.bind_eigen_trace`,
`EigenMetricBasis.operator/states`, `moments.py`, metric composites, state
bindings, and `reconstruct.shc/kubo_bastin`. Contrary to the suspected cause,
there is **no Python per-k loop in the per-moment path at b51a8c6** to remove.
The Python adapter loops only over plan instructions once when compiling;
C and S are whole `_core.BlockDiagonalOperator` applications. The cached
inverse/half-inverse and eigen transforms already operate on full stacked
arrays outside the recursion. The stage-6 `LowdinOperator` and
`MetricInverseProduct` call C++ block appliers. `BlockLocalStates` loops over
k in C++, and reconstruction contracts the moment tensor, with no k axis.
Python damping/grid loops in reconstruction are over moments/energies, not k.
No Python loop location or before/after profile numbers have been invented.

The cost removed is in `BlochMetricBasis.rescaled_hamiltonian`: previously
each Hamiltonian application traversed C++ LowdinOperator (three block
matvecs and two N-vector temporaries) or MetricInverseProduct (two block
matvecs and one N-vector temporary). It now reuses the already computed
Lowdin Hamiltonian, or caches the stacked generalised product once, and
makes **one existing C++ BlockDiagonalOperator call per Hamiltonian
application**. Its OpenMP loop owns all k blocks. Observable metric factors
continue to use the existing stage-6 C++ appliers. `MomentEngine.run` releases
the GIL for the whole native calculation. No recursion, weighting, seed
normalisation, or contraction formula was replaced. Caching changes floating
point association, so numerical agreement must be measured by the tester.

For comparable profiles, use separate writable worktrees at b51a8c6 and the
review chain, build offline on Aron, and run the same Bi Bloch/berry_connection
short descriptor with M=128, trace=random, states=2, seed=12345, oracle=false,
OMP_NUM_THREADS=48 and identical BLAS thread settings. Invoke
`python -m cProfile -o <unique-output>.prof -m lsquant run <short-run.json>`.
The phase receipt separates metric/operators, production moments,
reconstruction and eigen reference; cProfile cannot resolve native C++ costs.
Report the native run cumulative time, phase seconds, wall time and CPU before
and after. No performance or golden verdict is asserted here.

## Review 1 applied

This entry supersedes the earlier eigenphase gate/production description above.
The complete review was read, then AGENTS.md sections 8, 10, 11 and 12 and the
available brief sections 5 and 8. The checked-in brief ends at section 10;
there is no section 12 to re-read. The review and the user's explicit revisions
are the implementation instructions for this pass. The brief and review were
not edited or staged. No network, build, full developer self-check, material
run, profile, golden gate or figure was run here.

1. Gamma tolerance: `checked_hermitian` uses the largest block norm of the
   operator as its common scale; the tolerance is unchanged and real
   anti-Hermitian defects are still refused.
2. Generalised operators: physical observable blocks are checked before metric
   insertion; the eigen transform is checked/folded for accepted numerical
   round-off. The compiler's existing physical-operator check remains in place;
   no Euclidean-Hermiticity demand is added to an applied generalised matrix.
   Every production and reference word is prepared before any moment run,
   reconstruction, oracle or eigen-reference execution.
3. Exact curve grid: `reconstruction.energies` is an explicit finite increasing
   eV list; `reconstruction.energy_grid` is `[E_min, E_max, n]`. Both survive
   resolution/replay and are mutually exclusive. The default is the two band
   edges adjacent to the requested Fermi level, plus a margin (1001 points
   unless num_points is supplied), rather than the whole spectral window.
4. Progress/timing: metric build, dipole quadrature, operators, moments,
   reconstruction, oracle (when enabled), and eigen reference flush phase
   boundaries and wall seconds to stdout and `run.resolved.json` under
   `receipt.phases`. Failures retain their phase status and elapsed seconds.
5. Cost: the static audit and exact replacement are documented above. Tester
   baseline: 71 min, 264% CPU, 95 threads; short-profile before/after and after
   runtime are unavailable here. No Python per-k loop existed in the audited
   per-moment path; metric Hamiltonian composites are now cached as one C++
   BlockDiagonalOperator application, and the engine releases the GIL.
6. Trace: `EigenMetricBasis.states("exact")` explicitly exposes the C++
   BlockLocalStates trace (W states, all k blocks). All 24 G11 descriptors now
   select basis=eigen and trace=exact, with no stochastic vector count.
   Every eigen reference is exact. Production random phases are in the
   orthonormal orbital Bloch representation: generalised seeds are mapped by
   the cached inverse square root; eigen execution transforms them once by
   the adjoint cached Lowdin eigenvector operator in C++ TransformedStates.
   Rank-2 self-check assertions use exact traces; stochastic comparisons are
   diagnostics. The self-check measures the production plateau standard
   deviation across independent seeded samples at the gate mesh. Formula
   docstrings retain their SymPy source citations; the unchanged SHC units
   now explicitly cite `metric_forms.spin_hall_normalisation`.
7. Report/transport: this "Review 1 applied" entry records the requested changes
   without certifying a gate. `.git` is read-only (index.lock creation refused),
   so seven commits are recorded with detached metadata under
   `/tmp/lsquant-stage6b-r1.git`, based on b51a8c6. No branch was created or
   switched; no push was attempted. The requested portable bundle is
   `/tmp/lsquant-stage6b-r1.bundle`, depending only on b51a8c6. The workspace
   branch reference remains at b51a8c6 until the maintainer imports the bundle.

Editing checks performed in this pass (not integration/golden evidence):

```text
Changed Python files: Python 3.10 AST syntax OK
G11 JSON: 24 descriptors, eigen/exact, no states count
Hermiticity source-body probe: Gamma accepted; real defect refused
Source-body probes: explicit/default energy grids and completed/failed receipts OK
Diff whitespace check: no output
```

The tiny probes extract the actual function/class bodies into a NumPy namespace;
they do not import the package or extension. SciPy and SymPy are unavailable.
The changed C++ binding and the full Python/C++ execution require the tester's
offline build on Aron. The Hamiltonian cache changes floating point association;
all exact-route moment residuals and prior goldens must be measured there.

Expected self-check output for the tester (not output obtained in this sandbox):

```text
(a)-(d) covariance, Hermiticity and Hellmann-Feynman within existing assertions
(e) sigma_xy=2.000000000000 spin_chern=2.000000000000 Z2=1
(e) explicit Lowdin oracle plateau=2.000000000000
SHC exact block Bloch-minus-eigen moments=<measured, required below 1e-10>
(f) generalised exact KPM M=256 Nk=12x12 plateau=<measured> oracle=<measured>
(f) Lowdin exact KPM M=256 Nk=12x12 plateau=<measured>
SHC orbital-phase same-sample Bloch-minus-eigen moments=<measured> ... (diagnostic only)
anticommutator exact plateau=<measured>; equality is NOT required by the derivation
DOS single eigenphase-minus-exact=<measured, required below 1e-12>
SHC orbital-phase noise sample=1/8 plateau=<measured>
...
SHC orbital-phase noise mesh=48x48 M=256 states=16 samples=8 mean=<measured> std=<measured> standard_error=<measured> exact_KPM=<measured> [e^2/h, sigma_z Pauli]
stage_6b_check: developer self-check passed (not a golden gate)
```

The first tester's pre-review exact 12x12 KPM was 2.0055 against its finite-mesh
oracle 1.9999; the retained acceptance is absolute error below 1e-2, not those
rounded numbers. The production standard deviation is deliberately not guessed
from the earlier eigenphase estimator. Defaults use 8 independent samples of
16 orbital-phase vectors at 48x48, M=256; the printed std is the sample standard
deviation of a 16-vector plateau estimate (ddof=1), and standard_error describes
the mean of those independent estimates. To resolve it more precisely and at
the full G11 moment count, run the same one requested self-check with
`--noise-samples 100 --noise-moments 1024`; mesh/vector/moment/sample counts are
printed so a short diagnostic cannot be mistaken for the full production case.
No stochastic tolerance replaces an exact comparison.

## Review 2: item 1 source localisation, before changing the sampler

Read review-2.md completely, AGENTS.md's word/MirrorToBra instruction and
stage-6 theory equations (2a)-(2c). Commit 78e0431 is absent locally, so the
review's quoted tester evidence is the available final-run evidence. No fetch,
build or numerical execution is attempted. review-1.md is absent here;
review-2.md remains untracked and neither review is staged.

The suggested bra/ket mismatch is **not found in the rank-2 source**.
At the review-1 HEAD, moment_engine.cpp:535-540 prepares `ketSeed=K r`, then
starts the bra from that very same `ketSeed`. The block adapter at
bases.py:1371-1383 inserts S before the compiler's bra operations. With
`X=S^-1`, the generalised rank-2 vectors are therefore

```
right_n = X B T_n K r
left_m  = S T_m X A X S K r = S T_m X A K r
<left_m|right_n> = r^+ K^+ S X A T_m X B T_n K r.
```

Here K=C for exact block-local states, K=Y=S^-1/2 for orbital random
phases; on the Lowdin route K=U or I and S is absent. These are consistent
metric duals of a single prepared random vector, and in the eigen
representation the result is `z^+ A_e T_m(eps) B_e T_n(eps) z`.
RandomPhaseStates::Next (state_generator.cpp:34-45) draws N independent
phases in one stream, with amplitude 1/sqrt(N); it does not regenerate
phases for the bra or repeat them across blocks. ApplyKpmWeights at
stage6b_siesta.py:76 applies `(2-delta_m0)(2-delta_n0)/R` once. Combined
with the reconstruction's N factor, its trace normalisation is correct.

Assertion (g) now reads the actual C++ generator vector through the small
copy-returning StateGenerator.next binding, then engine.run resets and
replays that same stream. On a 4x4 mesh, M=8, batch=3, it compares EVERY raw
moment with the NumPy quadratic form (absolute residual <1e-10), for both
fixtures, both routes and both execution adapters. It separately checks the
unsymmetrised KPM weights. This is an assertion for Aron, **not a claimed
local numerical pass**. It is placed before the expensive self-check work.

There is a distinct difference from stage 3: 6b explicitly disables the
trace-cyclic Hermitian projection of the moment tensor. A single stochastic
word need not satisfy `mu_mn=conj(mu_nm)`, although its exact trace does for
Hermitian physical A and B. The noise comparison currently compares two
different estimators, so the review's Jensen bound does not establish a
seed mismatch. Item 2 will restore the common estimator and make the
fresh-state bra semantics explicit; its numerical effect remains for Aron.

Review 2 item 2: static MirrorToBra now starts its bra program from the
original r for ranks 0, 1 and 2. bind_eigen_trace explicitly applies K on
BOTH copies, then S on the generalised bra, rather than depending on the
old implicit prepared-ket copy. Thus equations (2a)-(2c) hold with the dual
metric of the common K r; no phases are redrawn. This is a semantic repair,
not a claim that the earlier rank-2 adapter used two different samples.
The pre-existing TimeAxis register-copy extension is documented separately
and retained for the stage-4 projected/MSD words.

The actual estimator change is ApplyKpmWeights(..., True): stage 3 and 6b
now use the same exact trace-cyclic identity for Hermitian physical A,B,
`mu <- (mu+mu^+)/2`, with weights and state count applied once. In general
the raw per-vector word does not obey that identity. Keeping its unwanted
part through the Bastin contraction can expose large cancelling noise;
the tester must quantify how much this change removes. Assertion (g) still
checks the raw quadratic form, and additionally compares compute_moments
with the explicitly projected NumPy tensor. Exact block traces remain the
gate/reference; no stochastic threshold or measured result is invented.
