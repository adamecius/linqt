"""SIESTA LCAO reader and Bloch adapter, independent of Wannier inputs.

Conventions: eV, Angstrom, hbar=1; interleaved (orbital up, orbital down)
spinors; exp(+ik.(R+tau_j-tau_i)) atom gauge. Stored TSHS energies already
have E_F=0: the printed Fermi level is metadata and is never subtracted.
Only this module imports sisl; array fixtures need neither sisl nor SIESTA.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

from .bloch import bloch_sum, k_mesh, reciprocal_lattice
from .nonorthogonal import adjoint, checked_hermitian


def spinor_scalar(matrix):
    """Lift ``S_orb`` to interleaved spinors as ``S_orb tensor I_2``."""
    return np.kron(matrix, np.eye(2, dtype=complex))


def midpoint_connection(position, overlap, centres, r_cart):
    """Convert full position Q(R) to the Hermitian connection used with eq. (4).

    ``A_a,ij(R)=Q_a,ij(R)-(tau_i+tau_j+R)_a S_ij(R)/2``.
    Translation gives ``Q_a(-R)^+ = Q_a(R)-R_a S(R)``, hence
    ``A_a(-R)^+=A_a(R)``. The raw row-centred dipole
    ``D=Q-tau_i S`` obeys ``D(R)-D(-R)^+=(R+tau_j-tau_i) S(R)``;
    its Bloch anti-Hermitian part is ``-i dS``. Using D unmodified in
    the symmetric-dS velocity would count this part a second time.
    For the rotated intra-cell fixture Q'=T^+ diag(tau) T is the input;
    this conversion, not an independently guessed position, supplies A'.
    """
    midpoint = (np.asarray(centres)[None, :, None, :] +
                np.asarray(centres)[None, None, :, :] + np.asarray(r_cart)[:, None, None, :]) / 2
    return np.asarray(position) - np.moveaxis(midpoint, -1, 0) * overlap[None, ...]


@dataclass
class SiestaModel:
    """Plain real-space matrices: H(R), S(R), optional midpoint A_a(R).

    Blocks include spin (W=2*N_orb). ``dipoles`` has shape (3,N_R,W,W)
    and is NOT a full position matrix; use midpoint_connection for that.
    Non-SIESTA fixtures may provide physical_spin to encode a general basis
    rotation, in which case no tensor-product assumption is made.
    """
    uc: np.ndarray
    tau: np.ndarray
    r_int: np.ndarray
    hoppings: np.ndarray
    overlaps: np.ndarray
    dipoles: np.ndarray | None = None
    seed_name: str = "model"
    metadata: dict = field(default_factory=dict)
    input_tolerance: float = 1e-10
    physical_spin: dict | None = None

    def __post_init__(self):
        self.uc = np.asarray(self.uc, float)
        self.tau = np.asarray(self.tau, float)
        self.r_int = np.asarray(self.r_int, int)
        self.r_cart = self.r_int @ self.uc
        self.num_wann = len(self.tau)
        shape = (len(self.r_int), self.num_wann, self.num_wann)
        self.hoppings = np.asarray(self.hoppings, complex)
        self.overlaps = np.asarray(self.overlaps, complex)
        if self.uc.shape != (3, 3) or self.tau.shape != (self.num_wann, 3):
            raise ValueError("invalid lattice or orbital centres")
        if self.hoppings.shape != shape or self.overlaps.shape != shape:
            raise ValueError("H(R) and S(R) must have shape (N_R,W,W)")
        if self.dipoles is not None:
            self.dipoles = np.asarray(self.dipoles, complex)
            if self.dipoles.shape != (3, *shape):
                raise ValueError("dipoles must have shape (3,N_R,W,W)")
        self.metadata.update(spin_ordering="interleaved", gauge="atom", energy_unit="eV",
                             length_unit="Angstrom", stored_fermi_level=0.0)

    def _sum(self, matrices, k, derivative=None, cell_derivative=False):
        """Evaluate M(k), or ``sum i d_a exp(ik.d) M(R)`` in the atom gauge.

        For a cell derivative d_a=R_a but the result is still atom-dressed,
        so H, S and bare velocities always share the same representation.
        """
        if cell_derivative:
            weighted = 1j * self.r_cart[:, derivative, None, None] * matrices
            out = bloch_sum(self.r_cart, self.tau, weighted, k)
        else:
            out = bloch_sum(self.r_cart, self.tau, matrices, k, derivative=derivative)
        return checked_hermitian(out, name="SIESTA Bloch sum", tolerance=self.input_tolerance)

    def Hk(self, k):
        """Return ``sum_R exp(+ik.(R+tau_j-tau_i)) H(R)`` in eV."""
        return self._sum(self.hoppings, k)

    def Sk(self, k):
        """Return the dimensionless spinor overlap with the same Fourier phase."""
        return self._sum(self.overlaps, k)

    def derivatives(self, k, axis, *, cell=False):
        """Return ``(dH_a,dS_a)`` in eV Angstrom and Angstrom respectively."""
        return (self._sum(self.hoppings, k, axis, cell),
                self._sum(self.overlaps, k, axis, cell))

    def Ak(self, k, axis):
        """Return ``sum exp(ik.d) A_a(R)`` in Angstrom; refuse missing dipoles."""
        if self.dipoles is None:
            raise ValueError("covariant velocity requires SIESTA radial-orbital dipole quadrature")
        return self._sum(self.dipoles[axis], k)

    def k_mesh(self, supercell):
        """Return the n1-fastest mesh ``k=(i/N1,j/N2,l/N3) B``."""
        return k_mesh(self.uc, supercell)


def _pair_quadrature(orbital_i, orbital_j, displacement, order):
    """Integrate overlap and midpoint dipole with tensor Gauss-Legendre nodes.

    ``S_ij=integral phi_i(u)^* phi_j(u-d) du``;
    ``A_ij=integral phi_i(u)^* (u-d/2) phi_j(u-d) du``. The domain is the
    intersection of the two cutoff bounding boxes; psi vanishes outside each
    radial cutoff sphere. No overlap-derived approximation replaces a dipole.
    """
    ri, rj = float(orbital_i.R), float(orbital_j.R)
    if ri <= 0 or rj <= 0:
        raise ValueError("missing radial functions/cutoffs in .ion.xml")
    d = np.asarray(displacement, float)
    lower, upper = np.maximum(-ri, d-rj), np.minimum(ri, d+rj)
    if np.any(lower >= upper) or np.linalg.norm(d) >= ri+rj:
        return np.zeros(4, complex)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    axes = [0.5*(upper[a]+lower[a]) + 0.5*(upper[a]-lower[a])*nodes for a in range(3)]
    scales = [0.5*(upper[a]-lower[a])*weights for a in range(3)]
    result = np.zeros(4, complex)
    # One plane at a time bounds scratch memory by O(order^2).
    yz = np.array(np.meshgrid(axes[1], axes[2], indexing="ij")).reshape(2, -1).T
    wyz = (scales[1][:, None] * scales[2][None, :]).ravel()
    for x, wx in zip(axes[0], scales[0]):
        points = np.column_stack((np.full(len(yz), x), yz))
        product = (np.asarray(orbital_i.psi(points)).conj() *
                   np.asarray(orbital_j.psi(points-d))) * wx * wyz
        result[0] += product.sum()
        result[1:] += product @ (points-d/2)
    return result


def dipole_quadrature(orbitals, centres, r_cart, overlaps, support, *, order=24,
                      tolerance=1e-5, max_order=96):
    """Build A(R) once, certifying refinement and agreement with stored S(R).

    Compare orders n and 2n until max absolute change in S and A is below
    tolerance (dimensionless and Angstrom respectively), and require the
    quadrature overlap to match S(R) to that same absolute tolerance. Include
    every stored H/S pair, including zero S entries with nonzero dipoles.
    The estimate is a refinement estimate, not a rigorous quadrature bound.
    """
    if order < 2 or max_order < 2*order or tolerance <= 0:
        raise ValueError("invalid dipole quadrature order/tolerance")
    shape = overlaps.shape
    out = np.zeros((3, *shape), complex)
    worst_change, worst_overlap, largest_order = 0.0, 0.0, order
    for ir, i, j in np.argwhere(support):
        d = r_cart[ir] + centres[j] - centres[i]
        n = order
        previous = _pair_quadrature(orbitals[i], orbitals[j], d, n)
        while True:
            n *= 2
            current = _pair_quadrature(orbitals[i], orbitals[j], d, n)
            change = float(np.max(np.abs(current-previous)))
            if change <= tolerance or n >= max_order:
                break
            previous = current
        overlap_error = float(abs(current[0]-overlaps[ir, i, j]))
        if change > tolerance or overlap_error > tolerance:
            raise ValueError(f"dipole quadrature uncertified at R index {ir}, orbitals {i},{j}: "
                             f"order={n}, change={change:g}, overlap_error={overlap_error:g}")
        out[:, ir, i, j] = current[1:]
        worst_change = max(worst_change, change)
        worst_overlap = max(worst_overlap, overlap_error)
        largest_order = max(largest_order, n)
    return out, dict(method="tensor Gauss-Legendre on cutoff intersection", initial_order=order,
                     largest_order=largest_order, tolerance=tolerance,
                     max_refinement_error=worst_change, max_overlap_error=worst_overlap)


def load_siesta_model(directory, label=None, *, dipoles=False, quadrature_order=24,
                      quadrature_tolerance=1e-5, quadrature_max_order=96):
    """Read .fdf/.TSHS or .HSX and optional .ion.xml radial functions with sisl.

    sisl 0.16.4 SOC storage is (Re Huu,Re Hdd,Re Hud,Im Hud,Im Huu,
    Im Hdd,Re Hdu,Im Hdu). Lift the scalar overlap to spinors. Validate the
    reconstructed H/S against public Hk/Sk in BOTH gauges at a nonzero k;
    fail on an unknown layout rather than permuting or silently folding it.
    Round-off up to 1e-7 relative from TSHS is measured, recorded and folded
    in Bloch sums (notes section 11); no CSR tolerance is relaxed.
    """
    import sisl

    root = Path(directory).expanduser().resolve()
    if label is None:
        candidates = sorted(root.glob("*.fdf"))
        if len(candidates) != 1:
            raise ValueError("specify the SIESTA label when the directory has multiple .fdf files")
        label = candidates[0].stem
    fdf_path = root / f"{label}.fdf"
    if not fdf_path.is_file():
        raise FileNotFoundError(fdf_path)
    fdf = sisl.get_sile(str(fdf_path))
    binary = next((root / f"{label}{suffix}" for suffix in (".TSHS", ".HSX")
                   if (root / f"{label}{suffix}").is_file()), None)
    if binary is None:
        raise FileNotFoundError(f"no {label}.TSHS or {label}.HSX in {root}")
    reader = sisl.get_sile(str(binary))
    source = reader.read_hamiltonian()
    if source.orthogonal or not source.spin.is_spinorbit:
        raise ValueError("stage6b requires a nonorthogonal SIESTA SOC Hamiltonian")
    geometry = source.geometry
    no = int(geometry.no)
    cell = np.asarray(geometry.cell, float)
    r_int = np.asarray(geometry.lattice.sc_off, int)
    atom_indices = np.asarray(geometry.o2a(np.arange(no)), int)
    centres = np.asarray(geometry.xyz)[atom_indices]
    components = [source.tocsr(i) for i in range(8)]
    spatial_s = source.tocsr(source.S_idx)
    h_blocks, s_blocks, support = [], [], []
    for ir in range(len(r_int)):
        sl = slice(ir*no, (ir+1)*no)
        data = [c[:, sl].toarray() for c in components]
        block = np.zeros((2*no, 2*no), complex)
        block[0::2, 0::2] = data[0] + 1j*data[4]
        block[1::2, 1::2] = data[1] + 1j*data[5]
        block[0::2, 1::2] = data[2] + 1j*data[3]
        block[1::2, 0::2] = data[6] + 1j*data[7]
        overlap = spatial_s[:, sl].toarray()
        h_blocks.append(block)
        s_blocks.append(overlap)
        # Stored sparse pattern, not numerical nonzero S: parity may make S=0.
        pattern = spatial_s[:, sl].copy()
        pattern.data[:] = 1
        for component in components:
            extra = component[:, sl].copy()
            extra.data[:] = 1
            pattern = pattern + extra
        support.append(pattern.toarray() != 0)
    s_blocks = np.asarray(s_blocks)
    metadata = dict(sisl_version=sisl.__version__, input_file=str(binary), fdf=str(fdf_path),
                    energy_reference="stored H, E_F already zero", orbital_metadata=[])
    try:
        metadata["printed_fermi_level_eV"] = float(reader.read_fermi_level())
    except (AttributeError, NotImplementedError):
        metadata["printed_fermi_level_eV"] = None
    model = SiestaModel(cell, np.repeat(centres, 2, axis=0), r_int,
                        np.asarray(h_blocks), spinor_scalar(s_blocks), seed_name=label,
                        metadata=metadata, input_tolerance=1e-7)
    # A public-API comparison pins component order, phase sign and spin duplication.
    k_fractional = np.array([0.137, 0.219, 0.071])
    k_cart = k_fractional @ reciprocal_lattice(cell)
    checks = {}
    for gauge, cell_gauge in (("R", True), ("r", False)):
        for name, matrices in (("H", model.hoppings), ("S", model.overlaps)):
            reference = np.asarray(getattr(source, name+"k")(k_fractional, gauge=gauge,
                                                               format="array", dtype=np.complex128))
            if name == "S" and reference.shape == (no, no):
                reference = spinor_scalar(reference)
            rebuilt = bloch_sum(model.r_cart, model.tau, matrices, k_cart, cell=cell_gauge)
            error = float(np.max(np.abs(rebuilt-reference)))
            checks[name+gauge] = error
            if error > 1e-7 * max(1.0, float(np.max(np.abs(reference)))):
                raise ValueError(f"sisl {name}k layout/gauge validation failed: {error:g}")
            checks[name+gauge+"_anti_relative"] = float(
                np.linalg.norm(rebuilt-adjoint(rebuilt)) / (2*max(np.linalg.norm(rebuilt), 1e-300)))
    metadata["sisl_public_api_checks"] = checks
    orbitals = []
    for ia, atom in enumerate(geometry.atoms):
        tag = str(atom.tag)
        ion_path = root / f"{tag}.ion.xml"
        radial_atom = sisl.get_sile(str(ion_path)).read_basis() if dipoles else atom
        if radial_atom.no != atom.no:
            raise ValueError(f"radial orbital count differs for atom {ia} ({tag})")
        for orbital in radial_atom.orbitals:
            metadata["orbital_metadata"].append(dict(atom=ia, species=tag, orbital=str(orbital)))
            orbitals.append(orbital)
    if dipoles:
        A, certificate = dipole_quadrature(orbitals, centres, model.r_cart, s_blocks,
                                          np.asarray(support), order=quadrature_order,
                                          tolerance=quadrature_tolerance,
                                          max_order=quadrature_max_order)
        # Opposite-image pairs are integrated independently. Fold only after
        # the quadrature error has been certified; record the pre-fold residual.
        lookup = {tuple(r): i for i, r in enumerate(r_int)}
        paired = np.stack([adjoint(A[:, lookup[tuple(-r)]] ) for r in r_int], axis=1)
        certificate["pair_hermiticity_error"] = float(np.max(np.abs(A-paired)))
        if certificate["pair_hermiticity_error"] > 2*quadrature_tolerance:
            raise ValueError("dipole quadrature violates opposite-image Hermiticity")
        model.dipoles = spinor_scalar(0.5*(A+paired))
        metadata["dipole_quadrature"] = certificate
    return model
