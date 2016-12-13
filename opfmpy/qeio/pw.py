from __future__ import absolute_import, division, print_function
import os
import re
from xml.etree import cElementTree as ET

import numpy as np

from opfmpy.common.constants import Ha2Ry, Ha2eV, Bohr2Angstrom


class PWscfDataIO(ET.ElementTree):
    def __init__(self, data_file_xml):
        self.parse(data_file_xml)
        self.data_file_xml = data_file_xml

    def get_ecutwfc(self, units='Ry'):
        ecutwfc = float(self.find('./PLANE_WAVES/WFC_CUTOFF').text)

        if units == 'Ha':
            pass
        if units == 'Ry':
            ecutwfc *= Ha2Ry
        if units == 'eV':
            ecutwfc *= Ha2eV

        return ecutwfc

    def read_dlv(self):
        return read_dlv(self.data_file_xml)

    def read_atomic_positions(self, units='crystal'):
        return read_atomic_positions(self, units=units)

    def read_crystal_basis(self, units='crystal'):
        return read_crystal_basis(self, units=units)

    def read_kpoints(self, units='crystal'):
        return read_kpoints(self.data_file_xml, units=units)

    def read_eigenvalues(self, units='eV'):
        return read_eigenvalues(self.data_file_xml, units=units)


def read_crystal_basis(data_file_xml, units='crystal'):
    if type(data_file_xml) == str:
        data_file_xml = ET.parse(data_file_xml)

    dlv = read_dlv(data_file_xml)

    basis = []

    atom_tag_pattern = re.compile('ATOM.[\d]+')
    for elem in data_file_xml.findall('IONS/*'):
        if atom_tag_pattern.match(elem.tag) is not None:
            symbol = elem.attrib['SPECIES'].strip()
            tau = np.fromstring(elem.attrib['tau'], sep=' ')
            if units == 'cartesian':
                basis.append((symbol, tau))
            elif units == 'crystal':
                tau = np.dot(tau, np.linalg.inv(dlv))
                tau[np.abs(tau) < np.finfo(np.float).eps] = 0
                basis.append((symbol, tau))
            else:
                raise Exception

    return basis


def read_atomic_positions(data_file_xml, units='cartesian'):
    if type(data_file_xml) == str:
        data_file_xml = ET.parse(data_file_xml)

    dlv = read_dlv(data_file_xml)

    atomic_positions = []
    atom_tag_pattern = re.compile('ATOM.[\d]+')
    for elem in data_file_xml.findall('IONS/*'):
        if atom_tag_pattern.match(elem.tag) is not None:
            symbol = elem.attrib['SPECIES'].strip()
            tau = np.fromstring(elem.attrib['tau'], sep=' ')
            atomic_positions.append((symbol, tau))

    if units == 'crystal':
        atomic_positions = [(symbol, np.dot(tau, np.linalg.inv(dlv))) for (symbol, tau) in atomic_positions]
    elif units == 'cartesian':
        pass
    else:
        raise Exception

    return atomic_positions


def read_lattice_parameter(data_file_path, units='Bohr'):
    data_file_xml = ET.parse(data_file_path)

    try:
        units = str(units)
    except:
        raise TypeError('units must be able to be stringified')

    units = units.upper()
    angstrom_unit_labels = ['A', 'ANG', 'ANGSTROM']
    bohr_unit_labels = ['B', 'BOHR']
    possible_units = angstrom_unit_labels + bohr_unit_labels
    if units not in possible_units:
        raise ValueError('units must be one of the following %s' % repr(possible_units))

    lattice_parameter = float(data_file_xml.find('./CELL/LATTICE_PARAMETER').text)
    if units in angstrom_unit_labels:
        return lattice_parameter * Bohr2Angstrom
    elif units in bohr_unit_labels:
        return lattice_parameter
    else:
        raise ValueError('units must be one of the following %s' % repr(possible_units))


def read_fft_grid_dims(data_file_path):
    data_file_xml = ET.parse(data_file_path)
    fft_grid_dims = data_file_xml.find('./*/FFT_GRID').attrib
    return np.array(map(int, [fft_grid_dims['nr1'], fft_grid_dims['nr2'], fft_grid_dims['nr3']]))


def read_dlv(data_file_xml, units='bohr'):
    if type(data_file_xml) == str:
        data_file_xml = ET.parse(data_file_xml)

    a1 = np.fromstring(data_file_xml.find('./CELL/DIRECT_LATTICE_VECTORS/a1').text.strip(), sep=' ')
    a2 = np.fromstring(data_file_xml.find('./CELL/DIRECT_LATTICE_VECTORS/a2').text.strip(), sep=' ')
    a3 = np.fromstring(data_file_xml.find('./CELL/DIRECT_LATTICE_VECTORS/a3').text.strip(), sep=' ')

    dlv = np.array([a1, a2, a3])

    if units.upper() == 'BOHR':
        pass
    elif units.upper() == 'ANG':
        raise NotImplementedError
    else:
        raise ValueError('units must be "bohr" or "ang"')

    return dlv


def read_rlv(data_file_path):
    data_file_xml = ET.parse(data_file_path)
    alat = read_lattice_parameter(data_file_path, 'Bohr')
    b1 = np.fromstring(data_file_xml.find('./CELL/RECIPROCAL_LATTICE_VECTORS/b1').text.strip(), sep=' ')
    b2 = np.fromstring(data_file_xml.find('./CELL/RECIPROCAL_LATTICE_VECTORS/b2').text.strip(), sep=' ')
    b3 = np.fromstring(data_file_xml.find('./CELL/RECIPROCAL_LATTICE_VECTORS/b3').text.strip(), sep=' ')

    rlv = np.array([b1, b2, b3]) * 2 * np.pi / alat

    return rlv


def read_kpoints(data_file_path, units='crystal'):
    r"""
    Read k-points from data-file.xml

    Parameters
    ----------
    data_file_path: str
        path of data-file.xml

    Returns
    -------
    kpoints: ndarray


    Notes
    -----
    The k-points in data-file.xml are given in Cartesian coordinates in units of :math:`\frac{2\pi}{a}`.
    To convert to Cartesian requires multiplication by a factor of :math:`\frac{2\pi}{a}`.  This is
    done for units='cartesian'.  To convert to crystal coordinates:

    .. math::

        \frac{2\pi}{a} \mathbf{k}^T
        \left(
        \begin{array}{c}
        \mathbf{b}_1 \\ \mathbf{b}_2 \\ \mathbf{b}_3 \\
        \end{array}
        \right)

    The array of reciprocal lattice vectors is exactly that which is returned by :py:func:`read_rlv <opfmpy.qe.io.qexml.read_rlv>`

    """
    data_file_xml = ET.parse(data_file_path)
    alat = read_lattice_parameter(data_file_path)
    rlv = read_rlv(data_file_path)
    nkpts = int(data_file_xml.find('./BRILLOUIN_ZONE/NUMBER_OF_K-POINTS').text)

    kpoints = np.zeros((nkpts, 3))
    for ikpt in range(nkpts):
        kpoints[ikpt] = np.fromstring(data_file_xml.find('./BRILLOUIN_ZONE/K-POINT.%d' % (ikpt+1)).attrib['XYZ'], sep=' ')

    if units == 'crystal':
        kpoints = 2 * np.pi / alat * np.dot(kpoints, np.linalg.inv(rlv))
    elif units == 'cartesian':
        kpoints *= 2 * np.pi / alat
    else:
        raise ValueError('units must be "crystal" or "cartesian"')

    return kpoints


def read_eigenvalues(data_file_path, units='Ha'):
    data_file_path = os.path.abspath(data_file_path)
    data_file_dir = os.path.dirname(data_file_path)

    data_file_xml = ET.parse(data_file_path)

    nkpts = int(data_file_xml.find('./BAND_STRUCTURE_INFO/NUMBER_OF_K-POINTS').text)
    nbnds = int(data_file_xml.find('./BAND_STRUCTURE_INFO/NUMBER_OF_BANDS').text)
    nspns = int(data_file_xml.find('./BAND_STRUCTURE_INFO/NUMBER_OF_SPIN_COMPONENTS').text)

    if nspns == 2:
        eigenvalues = np.zeros((nkpts, nbnds, nspns))
    else:
        eigenvalues = np.zeros((nkpts, nbnds, 1))

    for elem in data_file_xml.findall('./EIGENVALUES/'):
        ikpt = int(elem.tag.split('.')[1]) - 1
        if nspns == 2:
            for ispn in range(nspns):
                eig_file_path = os.path.join(data_file_dir, elem.find('./DATAFILE.%d' % (ispn+1)).attrib['iotk_link'])
                eig_xml = ET.parse(eig_file_path)
                eigenvalues[ikpt, :, ispn] = np.fromstring(eig_xml.find('./EIGENVALUES').text.strip(), sep='\n')
        else:
            eig_file_path = os.path.join(data_file_dir, elem.find('./DATAFILE').attrib['iotk_link'])
            eig_xml = ET.parse(eig_file_path)
            eigenvalues[ikpt, :, 0] = np.fromstring(eig_xml.find('./EIGENVALUES').text.strip(), sep='\n')

    if units == 'Ha':
        pass
    elif units == 'Ry':
        eigenvalues *= Ha2Ry
    elif units == 'eV':
        eigenvalues *= Ha2eV
    else:
        raise ValueError('units must be "Ha", "Ry", or "eV"')

    return eigenvalues


def read_tetrahedra(data_file_path):
    """
    Read the tetrahedra from data-file.xml

    Parameters
    ----------
    data_file_path: str
        Path of data-file.xml

    Returns
    -------
    tetrahedra: ndarray

    """
    data_file_xml = ET.parse(data_file_path)
    ntetra = int(data_file_xml.find('./OCCUPATIONS/NUMBER_OF_TETRAHEDRA').text)

    tetrahedra = np.zeros((ntetra, 4), dtype=int)
    for elem in data_file_xml.findall('./OCCUPATIONS/'):
        if re.compile(r'TETRAHEDRON[.]\d+').match(elem.tag):
            itetra = int(elem.tag.split('.')[1]) - 1
            tetrahedra[itetra] = np.fromstring(elem.text.strip(), sep='\n')

    # make indices 0-based
    tetrahedra -= 1

    return tetrahedra


def read_fermi_energy(data_file_path, units='Ha'):
    data_file_xml = ET.parse(data_file_path)

    if units == 'Ha':
        energy_unit = 1.0
    elif units == 'Ry':
        energy_unit = Ha2Ry
    elif units == 'eV':
        energy_unit = Ha2eV
    else:
        raise ValueError('units must be "Ha", "Ry", or "eV"')

    elem = data_file_xml.find('./BAND_STRUCTURE_INFO/FERMI_ENERGY')

    if elem is not None:
        return float(elem.text) * energy_unit
    else:
        return None


def read_symmetries(data_file_xml):
    if type(data_file_xml) == str:
        data_file_xml = ET.parse(data_file_xml)

    nsyms = int(data_file_xml.find('./SYMMETRIES/NUMBER_OF_SYMMETRIES').text)

    symmetries = {}
    symmetries['translations'] = np.zeros((nsyms, 3))
    symmetries['rotations'] = np.zeros((nsyms, 3, 3))
    for elem in data_file_xml.findall('./SYMMETRIES/'):
        if re.compile(r'SYMM[.]\d+').match(elem.tag):
            i = int(elem.tag.split('.')[1]) - 1
            symmetries['translations'][i] = np.fromstring(elem.find('FRACTIONAL_TRANSLATION').text, sep='\n')
            symmetries['rotations'][i] = np.fromstring(elem.find('ROTATION').text, sep='\n').reshape((3, 3), order='F')

    return symmetries
