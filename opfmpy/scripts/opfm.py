#!/usr/bin/env python
"""Optimized projection functions method"""
from __future__ import absolute_import, division, print_function
import os
import errno
import sys
import re
import pickle
import argparse
import pprint
import contextlib
import shutil
from collections import OrderedDict
from ConfigParser import ConfigParser
from cStringIO import StringIO

import numpy as np
import scipy.linalg
from w90utils import io as w90io
import w90utils
import codiag
import f90nml
try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import lines as mlines
except ImportError:
    matplotlib = None

from opfmpy.common import cli, CrystalConfigParser, KpointsConfigParser
from opfmpy.common.regex import regex_int
from opfmpy import qeio
from opfmpy import utils
#
if codiag.flib is not None:
    givens = codiag.flib.givens
else:
    givens = codiag.givens


def force_symlink(file1, file2):
    try:
        os.symlink(file1, file2)
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.remove(file2)
            os.symlink(file1, file2)


class MyConfigParser(CrystalConfigParser, KpointsConfigParser):
    atomic_orbtial_pattern = re.compile(r'^[A-Z][a-z]?[0-9]*[:]\s*([spdf]|(([spdf];)+[spdf]))$')
    coordination_number_pattern = re.compile(r'^[A-Z][a-z]?[0-9]*[:]\s*[0-9]+$')

    def __init__(self, *args, **kwargs):
        ConfigParser.__init__(self, *args, **kwargs)

    def _validate(self):
        self.precursor

        if self.precursor == 'QE':
            self.qe_outdir
            self.qe_prefix

        self.restart
        self.max_iter
        self.conv_tol
        self.lagrange_multiplier
        self.kgrid
        self.band_ranges
        self.kpoints
        self.atomic_orbitals
        self.exclude_bands
        self.spinor_orbitals
        assert np.prod(self.kgrid) == len(self.kpoints)

    @property
    def precursor(self):
        return self.get('opfm', 'precursor').upper()

    @property
    def qe_outdir(self):
        return self.get('opfm:qe', 'outdir')

    @property
    def qe_prefix(self):
        return self.get('opfm:qe', 'prefix')

    @property
    def projwfc_out(self):
        return self.get('opfm:qe', 'projwfc_out')

    @property
    def restart(self):
        return self.getboolean('opfm', 'restart')

    @property
    def max_iter(self):
        return self.getint('opfm', 'max_iter')

    @property
    def conv_tol(self):
        if self.has_option('opfm', 'conv_tol'):
            return self.getfloat('opfm', 'conv_tol')
        else:
            return None

    @property
    def lagrange_multiplier(self):
        return self.getfloat('opfm', 'lagrange_multiplier')

    @property
    def include_bweights(self):
        return self.getboolean('opfm', 'include_bweights')

    @property
    def include_offdiags(self):
        return self.getboolean('opfm', 'include_offdiags')

    @property
    def kgrid(self):
        return map(int, self.get('opfm', 'kgrid').strip().split())

    @property
    def band_ranges(self):
        if not self.has_option('opfm', 'band_ranges'):
            return None

        band_ranges_tmp = self.get('opfm', 'band_ranges')
        band_ranges = []
        for band_range in band_ranges_tmp.split():
            if not re.compile('{0}-{0}'.format(regex_int.pattern)).match(band_range):
                raise ValueError('bad value "%s", "band_ranges" must be in the form ibnd1-ibnd2' % band_range)

            ibnd1, ibnd2 = map(int, band_range.split('-'))
            if ibnd1 < 1 or ibnd2 < 1:
                raise ValueError('bad value "%s", band indices must be >= 1' % band_range)

            band_ranges.append((ibnd1, ibnd2))

        return band_ranges

    @property
    def _bnd_idx(self):
        if self.band_ranges is not None:
            bnd_idx = np.array([])
            for band_range in self.band_ranges:
                bnd_idx = np.concatenate((bnd_idx, np.arange(band_range[0], band_range[1]+1)))
            bnd_idx = bnd_idx.astype(int) - 1
        else:
            bnd_idx = None

        return bnd_idx

    @property
    def exclude_bands(self):
        return self.getboolean('opfm', 'exclude_bands')

    @property
    def spinor_orbitals(self):
        return self.getboolean('opfm', 'spinor_orbitals')

    @property
    def site_basis(self):
        try:
            site_basis = self.get('opfm', 'site_basis').upper()
            assert site_basis in ['DEFAULT', 'SPECIFY', 'COORDINATE']
        except AssertionError:
            raise

        return site_basis

    @property
    def orbital_basis(self):
        try:
            orbital_basis = self.get('opfm', 'orbital_basis').upper()
            assert orbital_basis in ['AUTO', 'SPECIFY', 'PSP']
        except AssertionError:
            raise

        return orbital_basis

    @property
    def atomic_orbitals(self):
        contents = self.get('opfm', 'atomic_orbitals')
        contents = contents.split('\n')
        symbols = contents[0].split()

        assert set(symbols) <= set(self.basis_symbols)

        atomic_orbitals = {}
        for line in contents[1:]:
            if not self.atomic_orbtial_pattern.match(line):
                raise ValueError(
                    '"{}" does not match "{}"'.format(
                        line, self.atomic_orbtial_pattern.pattern
                        )
                    )

            symbol = line.split(':')[0].strip()
            assert symbol in symbols
            orbitals = line.split(':')[1].split(';')

            if symbol in atomic_orbitals:
                raise ValueError('duplicate atomic symbol "%s"' % symbol)

            atomic_orbitals[symbol] = orbitals

        missing_symbols = list(set(symbols) - set(atomic_orbitals.keys()))
        if len(missing_symbols) > 0:
            raise ValueError('coordination number for "%s" not specified' % missing_symbols[0])

        return atomic_orbitals

    @property
    def nproj(self):
        atomic_orbitals = self.atomic_orbitals
        nspns = self.nspns

        nproj = 0
        for symbol in self.basis_symbols:
            if 's' in atomic_orbitals[symbol]:
                nproj += 1 * nspns
            if 'p' in atomic_orbitals[symbol]:
                nproj += 3 * nspns
            if 'd' in atomic_orbitals[symbol]:
                nproj += 5 * nspns

        return nproj

    @property
    def nproj_atom(self):
        atomic_orbitals = self.atomic_orbitals
        nspns = self.nspns

        nproj_atom = [[] for iatom in range(len(self.basis_symbols))]
        for iatom, symbol in enumerate(self.basis_symbols):
            nproj_atom[iatom] = 0
            if 's' in atomic_orbitals[symbol]:
                nproj_atom[iatom] += 1 * nspns
            if 'p' in atomic_orbitals[symbol]:
                nproj_atom[iatom] += 3 * nspns
            if 'd' in atomic_orbitals[symbol]:
                nproj_atom[iatom] += 5 * nspns

        return np.array(nproj_atom)

    @property
    def nspns(self):
        if self.spinor_orbitals:
            nspns = 2
        else:
            nspns = 1

        return nspns

    @property
    def coordination_numbers(self):
        if self.has_option('opfm', 'coordination_numbers'):
            contents = self.get('opfm', 'coordination_numbers')
            contents = contents.split('\n')
            symbols = contents[0].split()

            assert set(symbols) <= set(self.basis_symbols)

            coordination_numbers = {}
            for line in contents[1:]:
                if not self.coordination_number_pattern.match(line):
                    raise ValueError(
                        '"{}" does not match "{}"'.format(
                            line, self.coordination_number_pattern.pattern
                            )
                        )

                symbol = line.split(':')[0].strip()
                assert symbol in symbols
                coordination = int(line.split(':')[1])

                if symbol in coordination_numbers:
                    raise ValueError('duplicate atomic symbol "%s"' % symbol)

                coordination_numbers[symbol] = coordination

            missing_symbols = list(set(symbols) - set(coordination_numbers.keys()))
            if len(missing_symbols) > 0:
                raise ValueError('atomic orbitals for "%s" not specified' % missing_symbols[0])
        else:
            coordination_numbers = None

        return coordination_numbers

    @property
    def sites(self):
        if self.has_option('opfm', 'sites'):
            contents = self.get('opfm', 'sites')
            contents = contents.split('\n')
            n = int(contents[0])
            if len(contents) - 1 != n:
                raise Exception

            atom_idx = np.zeros((n,), dtype=int)
            Rvectors = np.zeros((n, 3), dtype=int)
            for (i, line) in enumerate(contents[1:]):
                atom_idx[i] = int(line.split()[0]) - 1
                Rvectors[i] = np.array(map(int, line.split()[1:]))
        else:
            return None

        return atom_idx, Rvectors

    @property
    def supercell(self):
        if self.has_option('opfm', 'supercell'):
            nk1, nk2, nk3 = map(int, self.get('opfm', 'supercell').split())
            Rvectors = np.mgrid[-nk1:nk1, -nk2:nk2, -nk3:nk3].reshape((3, -1)).transpose()
            atom_idx = np.tile(np.arange(len(self.basis)), len(Rvectors))
            Rvectors = reduce(lambda x, y: np.column_stack((x, y)), [Rvectors for i in range(len(self.basis))]).reshape((-1, 3))
        else:
            return None

        return atom_idx, Rvectors

    @property
    def initialize_wmat(self):
        if self.has_option('opfm', 'initialize_wmat'):
            return self.get('opfm', 'initialize_wmat').upper()
        else:
            return 'IDENTITY'

    def __str__(self):
        with contextlib.closing(StringIO()) as sio:
            self.write(sio)
            return sio.getvalue()


def read_config(cfg_file='opfm.cfg'):
    config = MyConfigParser(
        defaults={
            'restart': 'False',
            'exclude_bands': 'False',
            'spinor_orbitals': 'False',
            'include_bweights': 'True',
            'include_offdiags': 'True',
        })
    config.read(cfg_file)
    # config._validate()

    return config


def create_cfg(args):
    config = read_config()

    if config.precursor == 'QE':
        _create_cfg_qe(args)
    else:
        raise Exception


def _create_cfg_qe(args):
    config = read_config()

    outdir = config.qe_outdir
    prefix = config.qe_prefix

    data_file_xml = os.path.join(outdir, prefix+'.save', 'data-file.xml')
    pw_dat = qeio.PWscfDataIO(data_file_xml)
    dlv = pw_dat.read_dlv()
    crystal_basis = pw_dat.read_crystal_basis(units='crystal')
    kpoints = pw_dat.read_kpoints(units='crystal')

    config = MyConfigParser()
    if os.path.isfile('opfm.cfg'):
        if not args.force:
            print('WARNING: "opfm.cfg" already exists, use --force to overwrite')
            print('WARNING: the "[crystal]" and "[kpoints]" will be overwritten')
            exit()
        else:
            config.read('opfm.cfg')

    config.remove_section('Crystal')
    config.remove_section('Kpoints')

    config.set_dlv(dlv)
    config.set_basis(crystal_basis)
    config.set_kpoints(kpoints)

    with open('opfm.cfg', 'w') as f:
        config.write(f)

    if config.orbital_basis == 'PSP':
        write_amn = False
    else:
        write_amn = True

    # ------------------------------------------------------------------
    pw2opfm = {}
    pw2opfm['inputpp'] = {}
    pw2opfm['inputpp']['outdir'] = outdir
    pw2opfm['inputpp']['prefix'] = prefix
    pw2opfm['inputpp']['seedname'] = 'opfm'
    pw2opfm['inputpp']['write_amn'] = write_amn
    pw2opfm['inputpp']['write_mmn'] = True
    f90nml.write(pw2opfm, 'pw2opfm.in', force=True)
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    projwfc = {}
    projwfc['projwfc'] = {}
    projwfc['projwfc']['outdir'] = outdir
    projwfc['projwfc']['prefix'] = prefix
    projwfc['projwfc']['filproj'] = 'atomic_proj'
    f90nml.write(projwfc, 'projwfc.in', force=True)
    # ------------------------------------------------------------------


def projwfc2amn(args):
    config = read_config()

    outdir = config.qe_outdir
    prefix = config.qe_prefix

    proj = qeio.projwfc.read_atomic_proj(os.path.join(outdir, prefix+'.save', 'atomic_proj.xml'))
    w90io.write_amn('opfm.amn', proj.conj())


def print_header():
    print('='*80)
    print()
    print('OPFm')
    print()
    print('Run parameters:')
    print('  flib: %s' % str(codiag.flib is not None))
    print()
    print('='*80)


def read_input(args):
    config = read_config()

    print('Reading input arrays')

    print('Reading AMN ...', end=' '); sys.stdout.flush()
    a_bloch = w90io.read_amn('opfm' + '.amn')
    print('done')

    print('Reading MMN ...', end=' '); sys.stdout.flush()
    m_bloch = w90io.read_mmn('opfm' + '.mmn')
    print('done')

    print()

    if config._bnd_idx is not None:
        bnd_idx = config._bnd_idx
        a_bloch = np.copy(a_bloch[:, bnd_idx, :])
        m_bloch = np.copy(m_bloch[:, :, bnd_idx[:, np.newaxis], bnd_idx])

    return a_bloch, m_bloch


# http://stackoverflow.com/questions/5889142/python-numpy-scipy-finding-the-null-space-of-a-matrix
def null(A, eps=1e-15):
    u, s, vh = scipy.linalg.svd(A)
    null_mask = (s <= eps)
    null_space = scipy.compress(null_mask, vh, axis=0)
    return scipy.transpose(null_space)


def construct_site_basis():
    config = read_config()

    dlv = config.dlv
    crystal_basis = config.basis
    basis_symbols = config.basis_symbols
    basis_vectors = config.basis_vectors

    if config.site_basis == 'COORDINATE':
        if config.coordination_numbers is None:
            raise Exception

        coordination_numbers = [config.coordination_numbers[symbol] for symbol, tau in crystal_basis]
        atom_idx, Rvectors = utils.coordinate(dlv, basis_vectors, coordination_numbers)
    elif config.site_basis == 'SPECIFY':
        atom_idx, Rvectors = config.sites
    elif config.site_basis == 'DEFAULT':
        atom_idx = np.arange(len(basis_symbols))
        Rvectors = np.zeros((len(basis_symbols), 3))
    else:
        raise Exception

    return atom_idx, Rvectors


def write_sites(atom_idx, Rvectors):
    config = read_config()

    dlv = config.dlv
    basis_symbols = config.basis_symbols
    basis_vectors = config.basis_vectors

    with open('opfm.sites.xyz', 'w') as f:
        print(len(atom_idx), file=f)
        print('STRUCTURE', file=f)
        for i in range(len(atom_idx)):
            R = np.dot(basis_vectors[atom_idx[i]] + Rvectors[i], dlv)
            print(basis_symbols[atom_idx[i]], '%12.6f %12.6f %12.6f' % tuple(R*.529177), file=f)


def check_sites(args):
    atom_idx, Rvectors = construct_site_basis()
    write_sites(atom_idx, Rvectors)


def construct_orbital_basis():
    config = read_config()

    atom_idx, Rvectors = construct_site_basis()

    if config.orbital_basis == 'SPECIFY':
        projections = w90utils.io.nnkp.read_projections('opfm.nnkp')
    elif config.orbital_basis == 'PSP':
        if config.precursor == 'QE':
            projections = qeio.projwfc.read_atomic_states(config.projwfc_out)
    else:
        raise Exception

    projections = np.split(np.array(projections), np.cumsum(config.nproj_atom))
    projections = reduce(lambda a, b: np.concatenate((a, b)), [projections[i] for i in atom_idx])
    orbitals = [proj['orbital'] for proj in projections]

    site_idx = [(i+1)*np.ones(n, dtype=int) for (i, n) in zip(np.arange(len(atom_idx)), config.nproj_atom[atom_idx])]
    site_idx = reduce(lambda a, b: np.concatenate((a, b)), site_idx)

    atom_idx = [(i+1)*np.ones(n, dtype=int) for (i, n) in zip(atom_idx, config.nproj_atom[atom_idx])]
    atom_idx = reduce(lambda a, b: np.concatenate((a, b)), atom_idx)

    Rvectors = Rvectors[site_idx-1]

    return site_idx, atom_idx, Rvectors, orbitals


def write_basis(site_idx, atom_idx, Rvectors, orbitals):
    data = np.column_stack((np.arange(len(site_idx))+1, site_idx, atom_idx, orbitals, Rvectors))
    np.savetxt('opfm.orbitals.dat', data, fmt='%5s')


def check_basis(args):
    site_idx, atom_idx, Rvectors, orbitals = construct_orbital_basis()
    write_basis(site_idx, atom_idx, Rvectors, orbitals)


def expand_projections(a):
    config = read_config()

    atom_idx, Rvectors = construct_site_basis()

    kpoints = config.kpoints
    nproj_atom = config.nproj_atom

    a1 = w90utils.expand_amn(a, kpoints, atom_idx, Rvectors, nproj_atom[atom_idx])

    return a1


def unitarize_aw(a_bloch, w_codiag):
    aw = np.einsum('...il,...lj->...ij', a_bloch, w_codiag)
    U, _, V = np.linalg.svd(aw, full_matrices=False)
    return np.einsum('...ik,...kj->...ij', U, V)


def update_iter_info(w90dat, m_proj, w_codiag, iter_info, initial=False):
    nkpts, nntot, nwann, nwann = w90dat.mmn.shape
    nkpts, nbnds, nproj = w90dat.amn.shape
    bvectors, bweights = w90dat.bv, w90dat.bw

    w_codiag = w_codiag[:, :nwann]

    m = w90utils.rotate_mmn(w90dat.mmn, unitarize_aw(w90dat.amn, w_codiag), w90dat.kpb_kidx)
    # aw = np.einsum('...ik,...kj->...ij', w90dat.amn, w_codiag)
    aw = np.dot(w90dat.amn, w_codiag)
    s = np.einsum('...ki,...kj->...ij', aw.conj(), aw)
    sI = s - np.eye(nwann)[np.newaxis]
    sI_diag2 = np.sum(np.abs(np.diagonal(sI, axis1=-2, axis2=-1))**2, axis=1)
    sigma = np.linalg.svd(s)[1]

    l1 = codiag.diag2(m_proj, nwann)
    l2 = -np.sum(iter_info['lm'] * sI_diag2)
    l3 = -np.sum(iter_info['lm'] * (np.linalg.norm(sI, axis=(-2, -1))**2 - sI_diag2))

    if initial:
        def _update(iter_info, key, value):
            iter_info[key] = [value]
        _update(iter_info, 'i', 0)
    else:
        def _update(iter_info, key, value):
            iter_info[key] += [value]
        _update(iter_info, 'i', iter_info['i'][-1] + 1)

    _update(iter_info, 'spread_d', w90utils.sprd.omega_d(m, bvectors, bweights))
    _update(iter_info, 'spread_od', w90utils.sprd.omega_od(m, bweights))
    _update(iter_info, 'spread_d+od', w90utils.sprd.omega_dod(m, bvectors, bweights))
    _update(iter_info, 'spread_i+od', w90utils.sprd.omega_iod(m, bweights))
    _update(iter_info, 'spread_tot', w90utils.sprd.omega(m, bvectors, bweights))
    _update(iter_info, 'lagrangian', l1+l2+l3)
    _update(iter_info, 'lagrangian_1', l1)
    _update(iter_info, 'lagrangian_2', l2)
    _update(iter_info, 'lagrangian_3', l3)
    _update(iter_info, 'sigma_min', np.min(sigma))
    # iter_info['elements_d'] = [codiag.diag2(m_proj, nwann)]
    # iter_info['elements_od'] = [codiag.offdiag2(m_proj, nwann)]
    # iter_info['overlap_d'] = [codiag.diag2(s, nwann)]
    # iter_info['overlap_od'] = [codiag.offdiag2(s, nwann)]


def print_iter_info(iter_info):
    def print_header(*args):
        format_string = '{:4s}'+'    {:>13s}'*len(args[1:])
        print(format_string.format(*args))
        print('-'*80)

    def print_row(*args):
        format_string = '{:4s}'+'    {:>13.8f}'*len(args[1:])
        if args[0] == '' or args[0] == 'new':
            data = [iter_info[s][-1] for s in args[1:]]
        if args[0] == 'old':
            data = [iter_info[s][-2] for s in args[1:]]
        if args[0] == 'diff':
            data = [iter_info[s][-1]-iter_info[s][-2] for s in args[1:]]

        print(format_string.format(args[0], *data))

    print()
    print('iter = %5d' % iter_info['i'][-1])
    print('='*80)

    print_header(
        '',
        'spread (D)',
        'spread (I+OD)',
        'spread (TOT)',
        'lagrangian',
        # 'lagrangian_1',
        # 'lagrangian_2',
        # 'lagrangian_3',
        )

    data_keys = [
        'spread_d',
        'spread_i+od',
        'spread_tot',
        'lagrangian',
        # 'lagrangian_1',
        # 'lagrangian_2',
        # 'lagrangian_3',
        ]

    if iter_info['i'][-1] == 0:
        print_row('', *data_keys)
        print('='*80)
        print()
        return

    print_row('old', *data_keys)
    print_row('new', *data_keys)
    print_row('diff', *data_keys)

    print('='*80)
    print()


def setup_opfm(args):
    config = read_config()
    dlv = config.dlv
    crystal_basis = config.basis
    kpoints = config.kpoints
    atomic_orbitals = config.atomic_orbitals
    nproj = config.nproj

    if args.verbose > 0:
        print(config)

    if args.verbose > 0:
        print('Crystal')
        print('-------')
        print('dlv=\n', dlv)
        print('crystal_basis=\n', crystal_basis)
        print('atomic_orbitals=\n', atomic_orbitals)
        print()

    # write Wannier90 input file in order to construct NNKP file
    # --------------------------------------------------------------------------
    with open('opfm.win', 'w') as f:
        print('NUM_BANDS = %d' % nproj, file=f)
        print('NUM_WANN  = %d' % nproj, file=f)
        print('', file=f)

        if config.spinor_orbitals:
            print('SPINORS = TRUE', file=f)
        print('BEGIN PROJECTIONS', file=f)
        for symbol, tau in crystal_basis:
            print(('f=%f,%f,%f:'+';'.join(atomic_orbitals[symbol])) % (tau[0], tau[1], tau[2]), file=f)
        print('END PROJECTIONS', file=f)
        print('', file=f)

        w90io.win.print_unit_cell(dlv, file=f)
        w90io.win.print_atoms(crystal_basis, units='crystal', file=f)
        w90io.win.print_kpoints(kpoints, mp_grid=config.kgrid, file=f)
    # --------------------------------------------------------------------------


def diagnostics(w90dat, w_codiag):
    a_opt = np.einsum('...il,...lj->...ij', w90dat.amn, w_codiag)
    s = np.einsum('...ki,...kj->...ij', a_opt.conj(), a_opt)
    _, sigma, _ = np.linalg.svd(s)

    print('Diagnostics')
    print('='*80)
    with PdfPages('diagnostics.pdf') as diagnostics_pdf:
        # histogram of minimal singular value
        # ----------------------------------------------------------------------
        fig = pyplot.figure()
        ax = fig.gca()

        y, bins, patches = ax.hist(
            np.min(sigma, axis=1),
            bins=np.linspace(0, np.max(sigma), 20),
            histtype='stepfilled',
            weights=np.ones(len(s))/len(s)
            )
        ax.set_title(r'$\min\sigma_i=$'+('%E' % np.min(sigma)))
        ax.set_xlabel('smallest singluar value')
        ax.set_ylabel(r'fraction of $\mathbf{k}$ points')
        diagnostics_pdf.savefig()
        # ----------------------------------------------------------------------
        print('minimal singular value histogram')
        print('-'*40)
        np.savetxt(sys.stdout, np.column_stack((bins[:-1], bins[1:], y)), fmt='%f - %f  |  %f')
        print('-'*40)
        print('minimal singular value =', np.min(sigma))
        print()

        print('-'*80)
        data = (np.abs(s - np.eye(s.shape[1])[np.newaxis])**2).flatten()
        print('mean data=', np.mean(data))
        fig = pyplot.figure()
        ax = fig.gca()

        y, bins, _ = ax.hist(
            data,
            bins=np.linspace(0, 1, 21),
            histtype='stepfilled',
            )
        ax.set_xlim([0, 1])
        ax.set_yscale('log')
        ax.set_ylim([1e0, 10**np.ceil(np.log10(np.max(y)))])
        ax.set_xlabel(r'$\left|\widetilde{S}^{(\mathbf{k})}_{ij}\right|^2$')
        ax.set_ylabel('number of elements')
        diagnostics_pdf.savefig()
        print()
        print('='*80)


def optimize(args):
    """Construct optimal projection functions (OPFs)"""

    print_header()

    config = read_config()

    # read the input arrays
    # ==========================================================================
    a_bloch, m_bloch = read_input(args)
    a_bloch = expand_projections(a_bloch)

    w90dat = w90io.read_data(seedname='opfm', amn=a_bloch, mmn=m_bloch)

    kpb_kidx = w90dat.kpb_kidx
    bweights = w90dat.bw
    projections = w90io.nnkp.read_projections('opfm' + '.nnkp')
    # ==========================================================================

    # construct projected overlap matrix
    # ==========================================================================
    print('Constructing M_proj ...', end=' '); sys.stdout.flush()
    u_proj = w90utils.unitarize(a_bloch)
    m_proj = w90utils.rotate_mmn(m_bloch, u_proj, kpb_kidx)
    print('done')
    print()
    nwann = a_bloch.shape[1]
    # ==========================================================================

    nkpts, nntot, nproj = m_proj.shape[:3]

    # initialize w_codiag
    # ==========================================================================
    print('Initializing W_codiag ...', end=' '); sys.stdout.flush()
    if config.restart:
        print('Restarting with "%s"' % 'w_restart.npy')
        w_codiag = np.load('w_restart.npy')
    else:
        if config.initialize_wmat == 'RANDOM':
            # np.random.seed(0)
            w_codiag = utils.random_unitary(nproj)
        elif config.initialize_wmat == 'IDENTITY':
            w_codiag = np.identity(nproj, dtype=complex)
        else:
            raise Exception
    assert np.all(np.abs(np.dot(w_codiag.conj().T, w_codiag) - np.eye(nproj)) < 1e-12)
    print('done')
    print()
    # ==========================================================================

    print('=' * 80)
    print('Constructing the %d optimal projections from a set of %d orbitals' % (nwann, nproj))
    print('Input data:')
    print('  Amn: %s' % repr(a_bloch.shape))
    print('  Mmn: %s' % repr(m_bloch.shape))
    print('Uproj: %s' % repr(u_proj.shape))
    print('Mproj: %s' % repr(m_proj.shape))
    print('=' * 80)

    # other setup
    # ==========================================================================
    print('Other setup ...', end=' '); sys.stdout.flush()
    m = np.copy(m_proj.reshape((nkpts*nntot, nproj, nproj)))
    aa = np.einsum('...ki,...kj->...ij', a_bloch.conj(), a_bloch)
    lfac = config.lagrange_multiplier

    if config.include_bweights:
        m *= np.sqrt(np.tile(bweights, nkpts)[:, np.newaxis, np.newaxis])
        lfac *= np.sum(bweights)

    if config.include_offdiags:
        def construct_Abc(x1, x2, x3): return x1-(x2+x3)
    else:
        def construct_Abc(x1, x2, x3): return x1-x2

    if not np.allclose(w_codiag, np.identity(nproj, dtype=complex)):
        print('Applying initial W0')
        # cannot use ellipses here for numpy<=1.8.2
        # m = np.einsum('kml,ln->kmn', m, w_codiag)
        # m = np.einsum('lm,kln->kmn', w_codiag.conj(), m)
        m = np.dot(m, w_codiag)
        m = np.tensordot(w_codiag.conj(), m, ([0], [1])).transpose(1, 0, 2)
        # aa = np.einsum('kml,ln->kmn', aa, w_codiag)
        # aa = np.einsum('lm,kln->kmn', w_codiag.conj(), aa)
        aa = np.dot(aa, w_codiag)
        aa = np.tensordot(w_codiag.conj(), aa, ([0], [1])).transpose(1, 0, 2)

    # aaI = aa - np.eye(nproj)[np.newaxis]
    aaI = np.sqrt(lfac)*(aa - np.eye(nproj)[np.newaxis])

    m = np.asfortranarray(m)
    aaI = np.asfortranarray(aaI)
    print('done')
    print()
    # ==========================================================================

    print('Simultaneous diagonalization of %d matrices of shape %dx%d' % (len(m), m.shape[1], m.shape[2]))

    # perform codiagonalization
    # ==========================================================================
    tol = config.conv_tol

    iter_info = OrderedDict()
    iter_info['lm'] = lfac
    update_iter_info(w90dat, m, w_codiag, iter_info, initial=True)
    print_iter_info(iter_info)

    iiter = 1
    while iiter <= config.max_iter:
        for i in range(nwann):
            kvals = np.delete(np.arange(nwann), i)
            for j in range(i+1, nproj):
                if j < nwann:
                    A1, b1, c1 = codiag.diag2_iijj_Abc(m, i, j)
                    A2, b2, c2 = codiag.diag2_iijj_Abc(aaI, i, j)
                    A3, b3, c3 = codiag.offdiag2_ijji_Abc(aaI, i, j)
                    A = construct_Abc(A1, A2, A3)
                    v = codiag.qpqc.solve_qp_s2(A, None)
                    if v[0] < 0:
                        v *= -1
                else:
                    A1, b1, c1 = codiag.diag2_ii_Abc(m, i, j)
                    A2, b2, c2 = codiag.diag2_ii_Abc(aaI, i, j)
                    _, b3, _ = codiag.offdiag2_ikki_Abc(aaI, i, j, kvals)
                    A = A1-A2
                    b = construct_Abc(b1, b2, b3)
                    v = codiag.qpqc.solve_qp_s2(A, b)

                    if v is None:
                        print('INFO: no solution')
                        continue

                c, s = codiag.v2cs(v)

                try:
                    assert abs(abs(c)**2 + abs(s)**2 - 1) < 1e-3
                except AssertionError:
                    print('WARNING: nonunitary Givens rotation (|c|^2+|s|^2=%f)' % (abs(c)**2 + abs(s)**2))
                    continue

                codiag.givens.right_multiply(w_codiag, i, j, c, s)

                givens.left_multiply(aaI, i, j, c, s)
                givens.right_multiply(aaI, i, j, c, s)

                givens.left_multiply(m, i, j, c, s)
                givens.right_multiply(m, i, j, c, s)

        if not args.quiet:
            update_iter_info(w90dat, m, w_codiag, iter_info)
            print_iter_info(iter_info)

        lagrangian_diff = iter_info['lagrangian'][-1] - iter_info['lagrangian'][-2]
        if abs(lagrangian_diff) < tol:
            break

        with open('opfm.out.pkl', 'wb') as f:
            pickle.dump(iter_info, f)

        # write restart file
        np.save('w_restart.npy', w_codiag)

        iiter += 1

    if iiter > config.max_iter:
        print('WARNING: maximum number of iterations reached')
    # ==========================================================================

    # # write restart file
    # np.save('w_restart.npy', w_codiag)

    w_codiag = w_codiag[:, :nwann]

    # check/fix unitarity of codiagonalization matrix
    try:
        assert np.allclose(np.dot(w_codiag.conj().T, w_codiag), np.eye(nwann, dtype=complex), rtol=0, atol=1e-3)
    except:
        print('')
        print('WARNING: w_codiag not unitary')
        print(np.max(np.abs(np.dot(w_codiag.conj().T, w_codiag) - np.eye(nwann, dtype=complex))))

        # unitarize w_codiag
        u, s, v = np.linalg.svd(w_codiag)
        w_codiag = np.dot(np.dot(u, np.eye(nproj, nwann)), v)
        print('')

    # plot/write w_codiag
    # --------------------------------------------------------------------------
    # # if not args.skip_plot_wmat:
    # print('Plotting/writing w_codiag ...', end=' '); sys.stdout.flush()
    # _plot_wmat('w_codiag.svg', w_codiag, projections)
    # _plot_wmat('w_codiag.pdf', w_codiag, projections)
    # print('done')
    # --------------------------------------------------------------------------

    # diagnostics
    # ==========================================================================
    diagnostics(w90dat, w_codiag)
    # ==========================================================================

    with open('opfm.out.pkl', 'wb') as f:
        pickle.dump(iter_info, f)

    np.save('w_codiag.npy', w_codiag)


def setup_lambda_scan(args):
    for l in args.lambda_values:
        lambda_dirname = 'lambda_{}'.format(l)

        if os.path.isdir(lambda_dirname):
            if not args.force:
                print('WARNING: skipping "{}", use --force t overwrite'.format(lambda_dirname))
                continue
            else:
                shutil.rmtree(lambda_dirname)

        os.mkdir(lambda_dirname)

        shutil.copy('opfm.cfg', os.path.join(lambda_dirname, 'opfm.cfg'))
        for fname in ['opfm.nnkp', 'opfm.amn', 'opfm.mmn', 'opfm.eig']:
            force_symlink('../%s' % fname, os.path.join(lambda_dirname, fname))
            # shutil.copy(fname, os.path.join(lambda_dirname, fname))

        config = read_config(os.path.join(lambda_dirname, 'opfm.cfg'))
        config.set('opfm', 'lagrange_multiplier', l)
        with open(os.path.join(lambda_dirname, 'opfm.cfg'), 'w') as f:
            config.write(f)


def plot_iter_info(args):
    with open(args.pkl, 'rb') as f:
        data = pickle.load(f)

    data_keys = ['spread_tot', 'lagrangian', 'sigma_min']

    with PdfPages('./iter.pdf') as pdf:
        for key in data_keys:
            fig = pyplot.figure()
            ax = fig.gca()
            ax.set_title(key.replace('_', '\_'))
            ax.set_xlabel('iter')
            ax.set_ylabel(key.replace('_', '\_'))

            ax.plot(
                np.arange(len(data[key])), data[key],
                linestyle='-', color='b',
                )

            pdf.savefig()


def run_diagnostics(args):
    print_header()

    config = read_config()

    # read the input arrays
    # ==========================================================================
    a_bloch, m_bloch = read_input(args)
    if config.expand_amn:
        a_bloch = expand_projections(args, a_bloch)

    w90dat = w90io.read_data(seedname='opfm', amn=a_bloch, mmn=m_bloch)

    kpb_kidx = w90dat.kpb_kidx
    bweights = w90dat.bw
    projections = w90io.nnkp.read_projections('opfm' + '.nnkp')
    # ==========================================================================

    w_codiag = np.load('w_codiag.npy')

    diagnostics(w90dat, w_codiag)


def opfm2w90(args):
    # if any(map(os.path.isfile, ['wannier.amn', 'wannier.mmn', 'wannier.eig'])):
    if os.path.isfile('wannier.amn') or os.path.isfile('wannier.mmn') or os.path.isfile('wannier.eig'):
        if not args.force:
            print('WARNING: Wannier90 files already exist, use --force to overwrite')
            exit()

    config = read_config()

    dlv = config.dlv
    crystal_basis = config.basis
    kpoints = config.kpoints

    kpb_kidx, kpb_g = w90io.nnkp.read_nnkpts('opfm.nnkp')

    w_codiag = np.load('w_codiag.npy')
    nwann = w_codiag.shape[1]

    # read the input arrays
    # ==========================================================================
    a_bloch, m_bloch = read_input(args)
    e_bloch = w90io.read_eig('opfm.eig')
    if config._bnd_idx is not None:
        e_bloch = e_bloch[:, config._bnd_idx]

    a_bloch = expand_projections(a_bloch)
    a_opt = np.einsum('...il,...lj->...ij', a_bloch, w_codiag)

    print('Output arrays')
    print('  Amn: %s' % repr(a_opt.shape))
    print('  Mmn: %s' % repr(m_bloch.shape))
    print('  Emn: %s' % repr(e_bloch.shape))

    print('Writing AMN and MMN and EIG ...', end=' '); sys.stdout.flush()
    w90io.write_amn('wannier'+'.amn', a_opt)
    w90io.write_mmn('wannier'+'.mmn', m_bloch, kpb_kidx, kpb_g)
    w90io.write_eig('wannier'+'.eig', e_bloch)
    print('done')

    # write WIN file
    # --------------------------------------------------------------------------
    if os.path.isfile('wannier.win'):
        if not args.force or not args.overwrite_win:
            print('WARNING: preserving "wannier.win", use --force and --overwrite-win to overwrite')
            exit()

    with open('wannier'+'.win', 'w') as f:
        print('NUM_BANDS = %d' % nwann, file=f)
        print('NUM_WANN  = %d' % nwann, file=f)
        print('', file=f)

        print('BEGIN PROJECTIONS', file=f)
        print('random', file=f)
        print('END PROJECTIONS', file=f)
        print('', file=f)

        w90io.win.print_unit_cell(dlv, file=f)
        w90io.win.print_atoms(crystal_basis, units='crystal', file=f)
        w90io.win.print_kpoints(kpoints, mp_grid=config.kgrid, file=f)
    # ==========================================================================


parser = argparse.ArgumentParser(description='Optimized Projection Functions Method')
parser_common = argparse.ArgumentParser(add_help=False, parents=[cli.parser_common])
subparsers = parser.add_subparsers(title='subcommands')
#
parser_create_cfg = subparsers.add_parser('create-cfg', parents=[parser_common])
parser_create_cfg.add_argument('--force', action='store_true')
parser_create_cfg.set_defaults(func=create_cfg)
#
parser_projwfc2amn = subparsers.add_parser('projwfc2amn', parents=[parser_common])
parser_projwfc2amn.set_defaults(func=projwfc2amn)
#
parser_check_sites = subparsers.add_parser('check-sites', parents=[parser_common])
parser_check_sites.set_defaults(func=check_sites)
#
parser_check_basis = subparsers.add_parser('check-basis', parents=[parser_common])
parser_check_basis.set_defaults(func=check_basis)
#
parser_setup = subparsers.add_parser('setup', parents=[parser_common])
parser_setup.set_defaults(func=setup_opfm)
#
parser_setup_lambda_scan = subparsers.add_parser('setup-lambda-scan', parents=[parser_common])
parser_setup_lambda_scan.add_argument('lambda_values', nargs='+', type=float)
parser_setup_lambda_scan.add_argument('--force', action='store_true')
parser_setup_lambda_scan.set_defaults(func=setup_lambda_scan)
#
parser_optimize = subparsers.add_parser('optimize', parents=[parser_common])
parser_optimize.add_argument('--restart', action='store_true')
parser_optimize.add_argument('--random-initial-w0', action='store_true')
parser_optimize.add_argument('--initial-w0')
parser_optimize.set_defaults(func=optimize)
#
parser_run_diagnostics = subparsers.add_parser('run-diagnostics', parents=[parser_common])
parser_run_diagnostics.set_defaults(func=run_diagnostics)
#
parser_plot_iter_info = subparsers.add_parser('plot-iter-info', parents=[parser_common])
parser_plot_iter_info.add_argument('pkl')
parser_plot_iter_info.set_defaults(func=plot_iter_info)
#
parser_opfm2w90 = subparsers.add_parser('opfm2w90', parents=[parser_common])
parser_opfm2w90.add_argument('--force', action='store_true')
parser_opfm2w90.add_argument('--overwrite-win', action='store_true')
parser_opfm2w90.set_defaults(func=opfm2w90)

if __name__ == '__main__':
    args = parser.parse_args()

    if args.verbose > 1:
        print()
        print('INPUT ARGUMENTS')
        print('='*80)
        pprint.pprint(vars(args))
        print('='*80)
        print()

    args.func(args)
