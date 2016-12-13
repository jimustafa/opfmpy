from __future__ import absolute_import, division, print_function
import re
import contextlib
import cStringIO as StringIO
from ConfigParser import ConfigParser

import numpy as np

from opfmpy.common.regex import regex_float


class CrystalConfigParser(ConfigParser):
    atomic_symbol_pattern = re.compile(r'[A-Z][a-z]?[0-9]*')
    atomic_position_pattern = re.compile(r'{0}\s+{1}\s+{1}\s+{1}'.format(atomic_symbol_pattern.pattern, regex_float.pattern))

    @property
    def dlv(self):
        if self.has_option('Crystal', 'lattice_constant'):
            alat = self.getfloat('Crystal', 'lattice_constant')
        else:
            alat = 1.0
        a1 = np.fromstring(self.get('Crystal', 'a1'), sep=' ') * alat
        a2 = np.fromstring(self.get('Crystal', 'a2'), sep=' ') * alat
        a3 = np.fromstring(self.get('Crystal', 'a3'), sep=' ') * alat
        return np.array([a1, a2, a3])

    def set_dlv(self, dlv):
        if not self.has_section('Crystal'):
            self.add_section('Crystal')
        self.set('Crystal', 'a1', '%18.12f %18.12f %18.12f' % tuple(dlv[0]))
        self.set('Crystal', 'a2', '%18.12f %18.12f %18.12f' % tuple(dlv[1]))
        self.set('Crystal', 'a3', '%18.12f %18.12f %18.12f' % tuple(dlv[2]))

    @property
    def basis(self):
        contents = self.get('Crystal', 'basis')
        contents = contents.split('\n')
        n = int(contents[0])
        if len(contents[1:]) != n:
            raise Exception

        basis = []
        for line in contents[1:]:
            if not self.atomic_position_pattern.match(line):
                raise ValueError(
                    '"{}" does not match "{}"'.format(
                        line, self.atomic_position_pattern.pattern
                        )
                    )

            symbol = line.split()[0]
            tau = np.array(map(float, line.split()[1:]))
            basis.append((symbol, tau))

        return basis

    def set_basis(self, basis):
        if not self.has_section('Crystal'):
            self.add_section('Crystal')

        with contextlib.closing(StringIO.StringIO()) as sio:
            print(len(basis), file=sio)
            for symbol, tau in basis:
                print('%5s %18.12f %18.12f %18.12f' % (symbol, tau[0], tau[1], tau[2]), file=sio)
            self.set('Crystal', 'basis', sio.getvalue().rstrip())

    @property
    def basis_symbols(self):
        return [symbol for symbol, tau in self.basis]

    @property
    def basis_vectors(self):
        return np.array([tau for symbol, tau in self.basis])


class KpointsConfigParser(ConfigParser):
    @property
    def kpoints(self):
        contents = self.get('Kpoints', 'kpoints')
        contents = contents.split('\n')
        n = int(contents[0])
        if len(contents) - 1 != n:
            raise Exception
        kpoints = np.zeros((n, 3))
        for (i, line) in enumerate(contents[1:]):
            kpoints[i] = np.array(map(float, line.split()))

        return kpoints

    def set_kpoints(self, kpoints):
        if not self.has_section('Kpoints'):
            self.add_section('Kpoints')

        with contextlib.closing(StringIO.StringIO()) as sio:
            print(len(kpoints), file=sio)
            np.savetxt(sio, kpoints, fmt='%18.12f'*3)
            self.set('Kpoints', 'kpoints', sio.getvalue().rstrip())
