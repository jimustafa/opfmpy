from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.constants import codata


# constants
alpha = codata.value('fine-structure constant')

# energy conversions
Ha2eV = codata.value('Hartree energy in eV')
eV2Ha = 1/Ha2eV
Ha2meV = Ha2eV * 1000
Ry2eV = codata.value('Rydberg constant times hc in eV')
Ha2Ry = 2.0
Ry2Ha = 0.5
Ha2K = codata.value('hartree-kelvin relationship')

# length conversions
Bohr2Angstrom = codata.value('Bohr radius') * 1e10
bohr2angstrom = codata.value('Bohr radius') * 1e10
Angstrom2Bohr = 1 / Bohr2Angstrom
angstrom2bohr = 1 / bohr2angstrom

# energy to inverse length conversions
Ha2cm = codata.value('hartree-inverse meter relationship') / 100
Ry2cm = codata.value('Rydberg constant') / 100

# energy to time conversions
Ry2Hz = codata.value('Rydberg constant times c in Hz')
Ry2THz = Ry2Hz / 1e12

# time to energy conversions
THz2eV = 1 / Ry2THz * Ry2eV
THz2meV = THz2eV * 1000
Ha2Hz = codata.value('hartree-hertz relationship')
Ha2THz = Ha2Hz / 1000
Ha2fs = 1 / Ha2Hz * 1e15
eV2fs = Ha2eV / Ha2Hz * 1e15

# energy/temperature
K2eV = codata.value('Boltzmann constant in eV/K')

# epsilons
eps = np.finfo(np.float).eps
eps3 = 1e-3
eps6 = 1e-6
eps9 = 1e-9
eps12 = 1e-12
