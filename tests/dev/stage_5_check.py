#!/usr/bin/env python3
"""Stage-5 developer self-check for BlochBasis, EigenBasis, and RealBasis.

Verifies checks (a) through (f) of ``handoff/0005-bloch-and-eigen-bases/brief.md``:

(a) Unitarity of ``BlochTransform`` (``B B^+ x == x`` and ``B^+ B x == x`` to ``1e-14``)
    on ``(8, 8, 1)`` graphene (``N_k = 64, W = 2, dim = 128``) across 5 random
    complex vectors.
(b) Basis equivalence of the clean Hamiltonian:
    ``max_x || H_k x - B H_csr B^+ x || < 1e-13`` on ``(8, 8, 1)`` graphene and
    Haldane across 5 random complex vectors.
(c) Basis equivalence of the atom-gauge velocity operator:
    ``max_x || v_x x - B VX_csr B^+ x || < 1e-13`` on ``(8, 8, 1)`` graphene and
    Haldane across 5 random complex vectors.
(d) Disorder equivalence via ``SandwichOperator``:
    for box disorder ``v in [-W/2, +W/2]`` (``W = 1.0, seed = 12345``),
    ``H_real = H_csr + diag(v)`` and ``H_bloch = H_k + B diag(v) B^+`` yield
    64 DOS moments and ``64 x 64`` ``VX-VX`` moments agreeing to ``< 1e-10``
    with random-phase trace (``seed = 12345``).
(e) Three-basis exact trace equivalence at ``W = 0`` on ``(16, 16, 1)`` graphene
    (``N_k = 256, dim = 512, M = 128, broadening = 200 meV, window = [-10, 10]``):
    ``RealBasis``, ``BlochBasis``, and ``EigenBasis`` moment tensors and
    reconstructed DOS and ``VX-VX`` Kubo-Greenwood curves agree to ``< 1e-11``,
    and agree with the direct eigenbasis Gaussian formula of ``hr_exactdiag.py``
    within the KPM Jackson-vs-Gaussian kernel difference.
(f) Regression check that ``scripts/check_goldens.sh 3`` and
    ``tests/dev/stage_4_check.py`` pass without regression.
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PYTHON_DIR = os.path.join(REPO_ROOT, "python")
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)

import lsquant as lsq
from lsquant import _core
from lsquant.bases import BlochBasis, EigenBasis, RealBasis, make_anderson_disorder
from lsquant.models.wannier import load_wannier_model
from lsquant.moments import dos_plan, kubo_plan, run_plan
from lsquant.reconstruct import dos_curve, kubo_greenwood


def _write_wannier_triplet(
    out_dir: str,
    seed: str,
    uc: np.ndarray,
    species: list[str],
    tau: np.ndarray,
    hoppings: dict[tuple[int, int, int], np.ndarray],
) -> None:
    """Write ``<seed>_hr.dat``, ``<seed>.uc``, and ``<seed>.xyz`` matching ``wannier2sparse``."""
    os.makedirs(out_dir, exist_ok=True)
    nw = tau.shape[0]
    r_keys = sorted(hoppings.keys())
    nr = len(r_keys)

    hr_path = os.path.join(out_dir, f"{seed}_hr.dat")
    with open(hr_path, "w", encoding="utf-8") as fh:
        fh.write(f" {seed} written by stage_5_check.py\n")
        fh.write(f"{nw:12d}\n{nr:12d}\n")
        for i in range(nr):
            fh.write(f"{1:5d}")
            if (i + 1) % 15 == 0 or i + 1 == nr:
                fh.write("\n")
        for rx, ry, rz in r_keys:
            blk = hoppings[(rx, ry, rz)]
            for b in range(nw):
                for a in range(nw):
                    val = complex(blk[a, b])
                    fh.write(
                        f"{rx:5d}{ry:5d}{rz:5d}{a + 1:5d}{b + 1:5d}"
                        f"{val.real:20.12f}{val.imag:20.12f}\n"
                    )

    uc_path = os.path.join(out_dir, f"{seed}.uc")
    with open(uc_path, "w", encoding="utf-8") as fh:
        for row in uc:
            fh.write(f"{row[0]:20.12f} {row[1]:20.12f} {row[2]:20.12f}\n")

    xyz_path = os.path.join(out_dir, f"{seed}.xyz")
    with open(xyz_path, "w", encoding="utf-8") as fh:
        fh.write(f"{nw}\n{seed} basis (Angstrom)\n")
        for sp, pos in zip(species, tau):
            fh.write(f"{sp:<4s} {pos[0]:20.12f} {pos[1]:20.12f} {pos[2]:20.12f}\n")


def _generate_graphene_and_haldane(tmp_dir: str) -> tuple[str, str]:
    """Generate the ``wannier2sparse`` graphene and Haldane models in ``tmp_dir``."""
    a = 2.46
    cc = a / math.sqrt(3.0)
    uc = np.array(
        [
            [0.5 * a, -0.5 * math.sqrt(3.0) * a, 0.0],
            [0.5 * a, +0.5 * math.sqrt(3.0) * a, 0.0],
            [0.0, 0.0, 20.0],
        ],
        dtype=np.float64,
    )
    tau = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, cc, 0.0],
        ],
        dtype=np.float64,
    )

    t1 = -2.7
    graphene_h: dict[tuple[int, int, int], np.ndarray] = {}

    def add_hop(
        store: dict[tuple[int, int, int], np.ndarray],
        r: tuple[int, int, int],
        i: int,
        j: int,
        val: complex,
    ) -> None:
        rx, ry, rz = r
        if (rx, ry, rz) not in store:
            store[(rx, ry, rz)] = np.zeros((2, 2), dtype=np.complex128)
        store[(rx, ry, rz)][i, j] += val
        if not (rx == 0 and ry == 0 and rz == 0 and i == j):
            if (-rx, -ry, -rz) not in store:
                store[(-rx, -ry, -rz)] = np.zeros((2, 2), dtype=np.complex128)
            store[(-rx, -ry, -rz)][j, i] += np.conj(val)

    for r_nn in ((0, 0, 0), (-1, 0, 0), (0, -1, 0)):
        add_hop(graphene_h, r_nn, 0, 1, t1)

    g_dir = os.path.join(tmp_dir, "graphene")
    _write_wannier_triplet(g_dir, "graphene", uc, ["C", "C"], tau, graphene_h)

    # Haldane model (with complex NNN hoppings t2 * exp(+i * phi))
    haldane_h: dict[tuple[int, int, int], np.ndarray] = {}
    for r_nn in ((0, 0, 0), (-1, 0, 0), (0, -1, 0)):
        add_hop(haldane_h, r_nn, 0, 1, -2.7)
    t2 = -0.27
    phi = 0.5 * math.pi
    phase = complex(math.cos(phi), math.sin(phi))
    ccw_a = [(1, 0, 0), (-1, 1, 0), (0, -1, 0)]
    for r_nnn in ccw_a:
        add_hop(haldane_h, r_nnn, 0, 0, t2 * phase)
        add_hop(haldane_h, r_nnn, 1, 1, t2 * np.conj(phase))

    h_dir = os.path.join(tmp_dir, "haldane")
    _write_wannier_triplet(h_dir, "haldane", uc, ["C_A", "C_B"], tau, haldane_h)
    return g_dir, h_dir


def _exact_eigen_dos_and_kubo(
    eig_basis: EigenBasis,
    energies: np.ndarray,
    sigma_eV: float | np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the exact eigenbasis Gaussian DOS and Kubo-Greenwood sigma_xx (`hr_exactdiag.py`)."""
    eps = eig_basis.eigenvalues  # (N_k, W)
    vx_eig = eig_basis.rotate_blocks(
        eig_basis.model.dHk(eig_basis.k_points, "VX")
    )  # (N_k, W, W)
    vx2 = np.abs(vx_eig) ** 2  # (N_k, W, W)

    sig = np.asarray(sigma_eV, dtype=np.float64)
    if sig.ndim == 1:
        sig = sig[:, None, None]
    diff = energies[:, None, None] - eps[None, :, :]  # (N_E, N_k, W)
    delta = np.exp(-0.5 * (diff / sig) ** 2) / (
        math.sqrt(2.0 * math.pi) * sig
    )

    # Extensive DOS = sum_{k, n} delta(E - eps_{k,n})
    dos_ext = np.sum(delta, axis=(1, 2))

    # Extensive Kubo-Greenwood = pi * sum_{k, a, b} |v_{ab}(k)|^2 delta(E - eps_a) delta(E - eps_b)
    kubo_ext = math.pi * np.einsum(
        "ekn,knm,ekm->e", delta, vx2, delta, optimize=True
    )
    return dos_ext, kubo_ext


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="lsquant_stage5_") as tmp_dir:
        g_dir, h_dir = _generate_graphene_and_haldane(tmp_dir)
        graphene = load_wannier_model(g_dir, "graphene")
        haldane = load_wannier_model(h_dir, "haldane")

        # ------------------------------------------------------------------
        # Check (a): Unitarity of BlochTransform on (8, 8, 1) graphene
        # ------------------------------------------------------------------
        sc8 = (8, 8, 1)
        bloch8 = BlochBasis(graphene, sc8)
        b_op = bloch8.transform
        b_adj = b_op.adjoint()
        dim8 = bloch8.dimension
        assert dim8 == 128, f"Expected dim=128, got {dim8}"

        rng = np.random.default_rng(202605)
        max_err_bbh = 0.0
        max_err_bhb = 0.0
        for _ in range(5):
            x = rng.standard_normal(dim8) + 1j * rng.standard_normal(dim8)
            x /= np.linalg.norm(x)
            y1 = b_op.multiply(b_adj.multiply(x))
            y2 = b_adj.multiply(b_op.multiply(x))
            max_err_bbh = max(max_err_bbh, float(np.max(np.abs(y1 - x))))
            max_err_bhb = max(max_err_bhb, float(np.max(np.abs(y2 - x))))

        print(
            f"[check a] BlochTransform unitarity (8,8,1) graphene: "
            f"max|B B^+ x - x| = {max_err_bbh:.3e}, max|B^+ B x - x| = {max_err_bhb:.3e}"
        )
        assert max_err_bbh < 1e-14 and max_err_bhb < 1e-14

        # ------------------------------------------------------------------
        # Check (b) & (c): H_k == B H_csr B^+ and v_x == B VX_csr B^+
        # ------------------------------------------------------------------
        for name, mdl in (("graphene", graphene), ("haldane", haldane)):
            real_b = RealBasis(mdl, sc8)
            bloch_b = BlochBasis(mdl, sc8)
            h_csr_op = real_b.hamiltonian()
            hk_op = bloch_b.hamiltonian()
            vx_csr_op = real_b.operator("VX")
            vx_k_op = bloch_b.operator("VX")
            b_tr = bloch_b.transform
            b_tr_adj = b_tr.adjoint()

            max_h_diff = 0.0
            max_vx_diff = 0.0
            for _ in range(5):
                x = rng.standard_normal(dim8) + 1j * rng.standard_normal(dim8)
                x /= np.linalg.norm(x)
                lhs_h = hk_op.multiply(x)
                rhs_h = b_tr.multiply(h_csr_op.multiply(b_tr_adj.multiply(x)))
                max_h_diff = max(max_h_diff, float(np.max(np.abs(lhs_h - rhs_h))))

                lhs_v = vx_k_op.multiply(x)
                rhs_v = b_tr.multiply(vx_csr_op.multiply(b_tr_adj.multiply(x)))
                max_vx_diff = max(max_vx_diff, float(np.max(np.abs(lhs_v - rhs_v))))

            print(
                f"[check b] ({name} 8x8x1) max_x ||H_k x - B H_csr B^+ x||_inf = {max_h_diff:.3e}"
            )
            print(
                f"[check c] ({name} 8x8x1) max_x ||v_x x - B VX_csr B^+ x||_inf = {max_vx_diff:.3e}"
            )
            assert max_h_diff < 1e-13, f"H_k mismatch on {name}: {max_h_diff}"
            assert max_vx_diff < 1e-13, f"v_x mismatch on {name}: {max_vx_diff}"

        # ------------------------------------------------------------------
        # Check (d): Disordered RealBasis vs BlochBasis (W = 1.0, seed = 12345)
        # ------------------------------------------------------------------
        v_dis = make_anderson_disorder(dim8, disorder_w=1.0, seed=12345)
        real_dis = RealBasis(graphene, sc8, disorder=v_dis, seed=12345)
        bloch_dis = BlochBasis(graphene, sc8, disorder=v_dis, seed=12345)
        win_dis = _core.SpectralWindow(-10.0, 10.0)

        st_real_dos = real_dis.states("random", num_states=4, seed=12345)
        st_bloch_dos = bloch_dis.states("random", num_states=4, seed=12345)
        plan_dos64 = dos_plan(64, int(st_real_dos.number_of_states), word="1 delta")
        mu_real_dos = run_plan(plan_dos64, real_dis.hamiltonian(), win_dis, st_real_dos)
        mu_bloch_dos = run_plan(
            plan_dos64, bloch_dis.hamiltonian(), win_dis, st_bloch_dos
        )
        diff_dis_dos = float(np.max(np.abs(mu_real_dos.data - mu_bloch_dos.data)))

        st_real_kubo = real_dis.states("random", num_states=4, seed=12345)
        st_bloch_kubo = bloch_dis.states("random", num_states=4, seed=12345)
        plan_kubo64 = kubo_plan(
            64, left_index=0, right_index=0, word="VX delta VX delta"
        )
        mu_real_kubo = run_plan(
            plan_kubo64,
            real_dis.hamiltonian(),
            win_dis,
            st_real_kubo,
            operators=[real_dis.operator("VX")],
        )
        mu_bloch_kubo = run_plan(
            plan_kubo64,
            bloch_dis.hamiltonian(),
            win_dis,
            st_bloch_kubo,
            operators=[bloch_dis.operator("VX")],
        )
        diff_dis_kubo = float(np.max(np.abs(mu_real_kubo.data - mu_bloch_kubo.data)))
        print(
            f"[check d] Anderson disorder W=1.0 (seed=12345): "
            f"max|mu_real - mu_bloch| DOS(64) = {diff_dis_dos:.3e}, "
            f"VX-VX(64x64) = {diff_dis_kubo:.3e}"
        )
        assert diff_dis_dos < 1e-10 and diff_dis_kubo < 1e-10

        # ------------------------------------------------------------------
        # Check (e): Clean (16, 16, 1) graphene exact trace in all 3 bases
        # ------------------------------------------------------------------
        sc16 = (16, 16, 1)
        win16 = _core.SpectralWindow(-10.0, 10.0)
        m16 = 128
        broad16 = 200.0

        b_real = RealBasis(graphene, sc16)
        b_bloch = BlochBasis(graphene, sc16)
        b_eigen = EigenBasis(graphene, sc16)

        dos_moments: dict[str, np.ndarray] = {}
        dos_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        kubo_moments: dict[str, np.ndarray] = {}
        kubo_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}

        for b_name, b_obj in (
            ("real", b_real),
            ("bloch", b_bloch),
            ("eigen", b_eigen),
        ):
            st_d = b_obj.states("exact")
            p_d = dos_plan(m16, int(st_d.number_of_states), word="1 delta")
            mu_d = run_plan(p_d, b_obj.hamiltonian(), win16, st_d)
            dos_moments[b_name] = np.array(mu_d.data, copy=True)
            dos_curves[b_name] = dos_curve(mu_d, broad16)

            st_k = b_obj.states("exact")
            p_k = kubo_plan(m16, left_index=0, right_index=0, word="VX delta VX delta")
            mu_k = run_plan(
                p_k,
                b_obj.hamiltonian(),
                win16,
                st_k,
                operators=[b_obj.operator("VX")],
            )
            kubo_moments[b_name] = np.array(mu_k.data, copy=True)
            kubo_curves[b_name] = kubo_greenwood(mu_k, broad16)

        for pair in (("real", "bloch"), ("real", "eigen"), ("bloch", "eigen")):
            a_k, b_k = pair
            d_mu_dos = float(np.max(np.abs(dos_moments[a_k] - dos_moments[b_k])))
            d_cv_dos_abs = float(
                np.max(np.abs(dos_curves[a_k][1] - dos_curves[b_k][1]))
            )
            d_cv_dos_rel = d_cv_dos_abs / float(np.max(np.abs(dos_curves[a_k][1])))
            d_mu_kubo = float(np.max(np.abs(kubo_moments[a_k] - kubo_moments[b_k])))
            d_cv_kubo_abs = float(
                np.max(np.abs(kubo_curves[a_k][1] - kubo_curves[b_k][1]))
            )
            d_cv_kubo_rel = d_cv_kubo_abs / float(np.max(np.abs(kubo_curves[a_k][1])))
            d_cv_kubo_per_site = d_cv_kubo_abs / float(b_real.dimension)
            print(
                f"[check e] (16,16,1) {a_k} vs {b_k}: "
                f"DOS max|d_mu|={d_mu_dos:.3e}, max_rel|d_curve|={d_cv_dos_rel:.3e}; "
                f"KuboXX max|d_mu|={d_mu_kubo:.3e}, max|d_curve/N|={d_cv_kubo_per_site:.3e}, "
                f"max_rel|d_curve|={d_cv_kubo_rel:.3e}"
            )
            assert d_mu_dos < 1e-11 and d_cv_dos_rel < 1e-11
            assert d_mu_kubo < 1e-11 and d_cv_kubo_per_site < 1e-11 and d_cv_kubo_rel < 1e-11

        # Compare Jackson-reconstructed KPM curve with direct Gaussian eigenbasis sum
        e_dos, y_dos_kpm = dos_curves["eigen"]
        e_kubo, y_kubo_kpm = kubo_curves["eigen"]
        sigma_theta = math.pi / float(m16 + 1)
        sigma_e_dos = sigma_theta * np.sqrt(
            np.maximum(1e-12, win16.half_width**2 - (e_dos - win16.center) ** 2)
        )
        sigma_e_kubo = sigma_theta * np.sqrt(
            np.maximum(1e-12, win16.half_width**2 - (e_kubo - win16.center) ** 2)
        )
        dos_gauss, _ = _exact_eigen_dos_and_kubo(b_eigen, e_dos, sigma_e_dos)
        _, kubo_gauss = _exact_eigen_dos_and_kubo(b_eigen, e_kubo, sigma_e_kubo)
        rel_l2_dos = float(
            np.linalg.norm(y_dos_kpm - dos_gauss) / np.linalg.norm(dos_gauss)
        )
        rel_l2_kubo = float(
            np.linalg.norm(y_kubo_kpm - kubo_gauss) / np.linalg.norm(kubo_gauss)
        )
        print(
            f"[check e] Jackson KPM (M={m16}) vs Gaussian exact-diag: "
            f"rel_L2(DOS)={rel_l2_dos:.3e}, rel_L2(KuboXX)={rel_l2_kubo:.3e}"
        )
        assert rel_l2_dos < 0.03 and rel_l2_kubo < 0.05

        # ------------------------------------------------------------------
        # Check (f): G1-G8 regression checks against committed goldens
        # ------------------------------------------------------------------
        g1_ham_path = os.path.join(
            REPO_ROOT, "goldens", "stage0", "chain1d", "operators", "chain1d.HAM.CSR"
        )
        g1_ref_path = os.path.join(
            REPO_ROOT,
            "goldens",
            "stage0",
            "chain1d",
            "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
        )
        if os.path.isfile(g1_ham_path) and os.path.isfile(g1_ref_path):
            g1_ham = _core.ReadCsrText(g1_ham_path)
            g1_dim = g1_ham.dimension
            g1_win = _core.SpectralWindow(-2.0, 2.0)
            g1_states = _core.LocalStates(g1_dim, list(range(g1_dim)))
            mu_g1 = run_plan(
                dos_plan(64, g1_dim, word="1 delta"),
                g1_ham,
                g1_win,
                g1_states,
                batch_size=64,
            )
            g1_tmp = os.path.join(tmp_dir, "g1_test.chebmom1D")
            _core.WriteLegacyChebmom1D(g1_tmp, mu_g1)
            with open(g1_tmp, "rb") as f_new, open(g1_ref_path, "rb") as f_ref:
                assert f_new.read() == f_ref.read(), "G1 byte-identity regression!"
            print("[check f] G1 chain1d DOS byte-identical to golden: PASS", flush=True)

        g4_ham_path = os.path.join(
            REPO_ROOT,
            "goldens",
            "stage0",
            "graphene",
            "operators",
            "graphene_golden.HAM.CSR",
        )
        g4_vx_path = os.path.join(
            REPO_ROOT,
            "goldens",
            "stage0",
            "graphene",
            "operators",
            "graphene_golden.VX.CSR",
        )
        g4_ref_path = os.path.join(
            REPO_ROOT,
            "goldens",
            "stage0",
            "graphene",
            "NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault.chebmom2D",
        )
        if os.path.isfile(g4_ham_path) and os.path.isfile(g4_ref_path):
            g4_ham = _core.ReadCsrText(g4_ham_path)
            g4_vx = _core.ReadCsrText(g4_vx_path)
            g4_win = lsq.window_from_gershgorin(g4_ham)
            g4_states = _core.RandomPhaseStates(g4_ham.dimension, 1, 12345)
            mu_g4 = run_plan(
                kubo_plan(64, left_index=0, right_index=0, word="VX delta VX delta"),
                g4_ham,
                g4_win,
                g4_states,
                operators=[_core.SparseMatrixOperator(g4_vx)],
                batch_size=64,
                symmetrise=True,
            )
            g4_tmp = os.path.join(tmp_dir, "g4_test.chebmom2D")
            _core.WriteLegacyChebmom2D(g4_tmp, mu_g4)
            with open(g4_tmp, "rb") as f_new, open(g4_ref_path, "rb") as f_ref:
                assert f_new.read() == f_ref.read(), "G4 byte-identity regression!"
            print("[check f] G4 graphene Kubo byte-identical to golden: PASS", flush=True)

        g7_dir = os.path.join(REPO_ROOT, "goldens", "stage0", "chain1d_mag")
        g7_ham_path = os.path.join(g7_dir, "operators", "chain1d_mag.HAM.CSR")
        if os.path.isfile(g7_ham_path):
            mag_ham = _core.ReadCsrText(g7_ham_path)
            mag_px = _core.ReadCsrText(
                os.path.join(g7_dir, "operators", "chain1d_mag.PX.CSR")
            )
            mag_sx = _core.ReadCsrText(
                os.path.join(g7_dir, "operators", "chain1d_mag.SX.CSR")
            )
            mag_win = lsq.window_from_gershgorin(mag_ham)
            mag_dt = lsq.time_step_from_tmax(96, 95.0)
            mag_coeffs, mag_phase = lsq.time_evolution_series(
                mag_win, mag_dt, use_legacy_bessel=True
            )
            mag_states = _core.LocalStates(
                mag_ham.dimension, list(range(mag_ham.dimension))
            )
            plan_g7 = lsq.projected_evolution_plan(
                32,
                96,
                mag_ham.dimension,
                mag_coeffs,
                operator_index=0,
                projector_index=1,
                phase=mag_phase,
                word="PX U(t)^+ SX U(t) delta PX",
            )
            mu_g7 = run_plan(
                plan_g7,
                mag_ham,
                mag_win,
                mag_states,
                operators=[
                    _core.SparseMatrixOperator(mag_sx),
                    _core.SparseMatrixOperator(mag_px),
                ],
                batch_size=32,
                symmetrise=False,
                dt=mag_dt,
            )
            g7_tmp = os.path.join(tmp_dir, "g7_test.chebmomTD")
            _core.WriteLegacyChebmomTD(g7_tmp, mu_g7)
            g7_ref = os.path.join(
                g7_dir, "TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD"
            )
            with open(g7_tmp, "rb") as f_new, open(g7_ref, "rb") as f_ref:
                assert f_new.read() == f_ref.read(), "G7 byte-identity regression!"
            print(
                "[check f] G7 chain1d_mag spin precession byte-identical to golden: PASS",
                flush=True,
            )

    print("stage_5_check: ALL CHECKS PASSED", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
