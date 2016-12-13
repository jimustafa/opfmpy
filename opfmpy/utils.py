from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.spatial import KDTree

from opfmpy.common.constants import eps


def crystal2cartesian(coords, lattice_vectors):
    """
    Convert vectors expressed in crystal coordinates to Cartesian coordinates

    Parameters
    ----------
    coords: ndarray, shape (..., npts, 3)
    lattice_vectors: ndarray, shape (3, 3)

    Returns
    -------
    ndarray

    """

    return np.dot(coords, lattice_vectors)


def cartesian2crystal(coords, lattice_vectors):
    """
    Convert vectors expressed in Cartesian coordinates to crystal coordinates

    Parameters
    ----------
    coords: ndarray, shape (..., npts, 3)
    lattice_vectors: ndarray, shape (3, 3)

    Returns
    -------
    ndarray

    """

    return np.dot(coords, np.linalg.inv(lattice_vectors))


def find_vector(pt_target, pts_in, idx=None, tol=eps):
    """
    Find the index of a point in an array of points

    Parameters
    ----------
    pt_target: array_like, shape (N,)
        point to find
    pts_in: array_like, shape (M,N)
        array of points in which to find target point
    idx: ndarray, optional
        precomputed indices, useful for when this function is called many times
    tol: float, optional
        tolerance for matching the point

    """
    pt_target = np.asarray(pt_target)
    pts_in = np.asarray(pts_in)

    if pt_target.ndim > 1:
        raise TypeError('expected 1D vector for pt_target')
    if pts_in.ndim != 2:
        raise TypeError('expected 2D array for pts_in')
    if len(pt_target) != pts_in.shape[-1]:
        raise TypeError('expected pts_in to be an array of vectors with shape of pt_target')

    diffs = pts_in[:, 0] - pt_target[0]
    diffs2 = diffs**2
    if idx is not None:
        idx = idx[diffs2 < tol]
    else:
        idx = (diffs2 < tol).nonzero()[0]

    diffs = pts_in[idx, 1] - pt_target[1]
    diffs2 = diffs**2
    idx = idx[diffs2 < tol]

    diffs = pts_in[idx, 2] - pt_target[2]
    diffs2 = diffs**2
    idx = idx[diffs2 < tol]

    return idx


def supercell(lattice_vectors, basis_vectors, dims):
    basis_vectors = np.dot(basis_vectors, lattice_vectors)
    supercell_lattice_vectors = lattice_vectors * np.asarray(dims)[:, np.newaxis]

    lattice_supercell = np.mgrid[:dims[0], :dims[1], :dims[2]].reshape((3, -1)).transpose()
    lattice_supercell = np.dot(lattice_supercell, lattice_vectors)

    supercell_basis_vectors = basis_vectors[:, np.newaxis, :] + lattice_supercell[np.newaxis, :, :]
    supercell_basis_vectors = supercell_basis_vectors.reshape((-1, 3))
    supercell_basis_vectors = np.dot(supercell_basis_vectors, np.linalg.inv(supercell_lattice_vectors))

    return supercell_lattice_vectors, supercell_basis_vectors


def coordinate(lattice_vectors, basis_vectors, coordination_numbers):
    """
    Coordinate each atom in the basis with a specfic number of neighbors

    Parameters
    ----------
    lattice_vectors : array_like, shape (3, 3)
    basis_vectors : array_like, shape (natms, 3)
        basis vectors in crystal coordinates
    coordination_numbers: array_like, shape (natms,)

    """
    assert len(coordination_numbers) == len(basis_vectors)
    natms = len(basis_vectors)

    idx = []
    Rvectors = []
    atomic_positions = []

    basis_vectors = crystal2cartesian(basis_vectors, lattice_vectors)

    lattice_supercell = np.mgrid[-2:3, -2:3, -2:3].reshape((3, -1)).transpose()
    lattice_supercell_kdtree = KDTree(lattice_supercell)

    crystal_supercell = basis_vectors[:, np.newaxis, :] + np.dot(lattice_supercell, lattice_vectors)[np.newaxis, :, :]
    crystal_supercell = crystal_supercell.reshape((-1, 3))
    crystal_supercell_kdtree = KDTree(crystal_supercell)

    # cluster basis atoms around the origin
    # minimize (u+A^T.R)^2 for real R
    # -----------------------------------------------------
    A = 2*np.dot(lattice_vectors, lattice_vectors.T)
    for (iatm, tau) in enumerate(basis_vectors):
        b = 2*np.dot(tau, lattice_vectors.T)
        Rvector = np.dot(np.linalg.inv(A), -b)
        (d, i) = lattice_supercell_kdtree.query(Rvector)
        Rvector = lattice_supercell[i]

        atomic_positions.append(basis_vectors[iatm]+np.dot(Rvector, lattice_vectors))
        idx.append(iatm)
        Rvectors.append(Rvector)
    # -----------------------------------------------------
    # for (iatm, tau) in enumerate(basis_vectors):
    #     Rvector = np.zeros((3,))

    #     atomic_positions.append(basis_vectors[iatm]+np.dot(Rvector, lattice_vectors))
    #     idx.append(iatm)
    #     Rvectors.append(Rvector)
    # -----------------------------------------------------

    for (iatm, tau) in enumerate(basis_vectors):
        if coordination_numbers[iatm] == 0:
            continue
        tau = tau + np.dot(Rvectors[iatm], lattice_vectors)
        d, i = crystal_supercell_kdtree.query(tau, coordination_numbers[iatm]+1)
        i = i[d > 0]
        atm_idx, T_idx = np.unravel_index(i, (natms, len(lattice_supercell)))

        for iatm_bond, Tv in zip(atm_idx, lattice_supercell[T_idx]):
            pos = np.dot(Tv, lattice_vectors) + basis_vectors[iatm_bond]
            if not find_vector(pos, atomic_positions).size:
                atomic_positions.append(pos)
                idx.append(iatm_bond)
                Rvectors.append(Tv)

    return np.array(idx), np.array(Rvectors)


def random_unitary(n):
    z = (np.random.randn(n, n) + 1j*np.random.randn(n, n))/np.sqrt(2.0)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r)
    q = np.dot(q, np.diag(d/np.abs(d)))
    return q
