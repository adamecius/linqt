"""Basis-independent Bloch sums; k in inverse Angstrom, lattice rows in Angstrom."""
from __future__ import annotations

import numpy as np


def reciprocal_lattice(lattice):
    """Return reciprocal rows ``B=2 pi (L^-1)^T``, so ``L B^T=2 pi I``."""
    return 2 * np.pi * np.linalg.inv(lattice).T


def k_mesh(lattice, supercell):
    """Return ``k=(i/N1,j/N2,l/N3) B``, with i varying fastest."""
    if len(supercell) != 3 or any(int(n) != n or n < 1 for n in supercell):
        raise ValueError("supercell must contain three positive integers")
    n1, n2, n3 = map(int, supercell)
    return np.array([[i/n1, j/n2, l/n3] for l in range(n3)
                     for j in range(n2) for i in range(n1)]) @ reciprocal_lattice(lattice)


def bloch_sum(r_cart, centres, matrices, k_points, *, derivative=None, cell=False):
    """Sum ``exp(+ik.(R+tau_j-tau_i)) M_ij(R)`` and optional Cartesian derivative.

    For ``cell=True`` omit tau from the phase and derivative. A derivative
    multiplies by ``+i (R+tau_j-tau_i)_a``. Work is O(N_k N_R W^2), keeping
    only a single k block of phases at a time (no N_k N_R W^2 temporary).
    """
    points = np.asarray(k_points, dtype=float)
    displacement = np.asarray(r_cart)[:, None, None, :] + np.zeros(
        (1, len(centres), len(centres), 3))
    if not cell:
        displacement += np.asarray(centres)[None, None, :, :] - np.asarray(centres)[None, :, None, :]
    data = np.asarray(matrices)
    if derivative is not None:
        data = 1j * displacement[..., derivative] * data
    blocks = [np.sum(np.exp(1j * (displacement @ k)) * data, axis=0)
              for k in np.atleast_2d(points)]
    out = np.ascontiguousarray(blocks, dtype=np.complex128)
    return out[0] if points.ndim == 1 else out
