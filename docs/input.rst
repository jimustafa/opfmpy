==========
Input File
==========

Description
===========

The OPFm.py input file is formatted according to the rules dictated by Python's
:py:mod:`ConfigParser`. Briefly, the format is similar to Microsoft Windows INI
files, being composed of sections starting with a ``[section]`` header and
followed by ``option=value`` (or ``option:value``) entries.

OPFm.py requires the input file ``opfm.cfg`` in the working directory. This file
contains runtime parameters for setting up and running the calculations and also
the necessary crystallographic information.

In addition to the parsing option values and converting to types supported using
:py:mod:`ConfigParser` (``int``, ``float``, and ``bool``), we define some
additional parsable types.

.. _dict:

simple_dict
-----------
The following format defines a option format for specifying a simple dictionary.
The strings following the ``simple_dict`` option are the keys for the dictionary, while
the subsequent lines give the key-value pairs that make up the dictionary.

::

   simple_dict:  key1 key2 key3 ...
     key1:  value1
     key2:  value2
     key3:  value3
       .       .
       .       .
       .       .

.. _array:

array
-----

A 2d array with shape (n, 3)

::

   array:  n
      a_11  a_12  a_13
      a_21  a_22  a_23
       .     .     .
       .     .     .
       .     .     .
      a_n1  a_n2  a_n3

.. _basis:

basis
-----

While the format for the options with ``simple_dict`` or ``array`` types are
somewhat general, we define also a particular type for parsing the option that
specifies the crystal basis.

::

   basis:  natms
     symbol1  tau1_b1  tau1_b2  tau1_b3
     symbol2  tau2_b1  tau2_b2  tau2_b3
         .        .        .        .
         .        .        .        .
         .        .        .        .


opfm.cfg
========

The input file is composed of three sections: ``[opfm]``, ``[crystal]``, and
``[kpoints]``. Additionally, a section corresponding to the code used to
construct the Bloch states and band structure should be provided. Currently,
only QuantumESPRESSO is supported, requiring the section ``[opfm:qe]`` to be
specified. The options available in each of the sections are given below. All
options are required, unless explicitly stated otherwise.

.. variable definition block
   name
   ~~~~
   :Type:
   :Required:
   :Default:
   :Description:
   :Examples:

``[opfm]`` options
------------------

precursor
~~~~~~~~~
:Type: str
:Required: Yes
:Default: None
:Description: Code used to generate the Bloch states and band structure
:Choices: QE

site_basis
~~~~~~~~~~
:Type: str
:Required: Yes
:Default: None
:Description: Flag to indicate how the site basis is generated
:Choices: default, specify, coordinate

orbital_basis
~~~~~~~~~~~~~
:Type: str
:Required: Yes
:Default: None
:Description: Flag to determine how the orbital basis is constructed
:Choices: auto, specify, psp

restart
~~~~~~~
:Type: bool
:Required: No
:Default: False
:Description: Flag to restart from a previous *W* matrix. The file
              ``w_restart.npy`` must be in the current directory.

max_iter
~~~~~~~~
:Type: integer
:Required: Yes
:Default: None
:Description: Maximum number of iterations.

conv_tol
~~~~~~~~
:Type: float
:Required: Yes
:Default: None
:Description: Convergence tolerance on the total spread.

lagrange_multiplier
~~~~~~~~~~~~~~~~~~~
:Type: float
:Required: Yes
:Default: None
:Description: Value for the parameter :math:`\lambda`.

band_ranges
~~~~~~~~~~~
:Type: TYPE
:Required: No
:Default: all bands
:Description: Bands for which OPFs are constructed.

spinor_orbitals
~~~~~~~~~~~~~~~
:Type: bool
:Required: No
:Default: False
:Description: Flag to indicate that orbitals carry spin.

.. _atomic_orbitals:

atomic_orbitals
~~~~~~~~~~~~~~~
:Type: simple_dict_
:Required: No
:Default: s;p;d on all atoms
:Description: Atomic orbitals on which to compute projections. The specification
              of the atomic orbitals is similar to specifying the projections in
              the WIN file for Wannier90.

coordination_numbers
~~~~~~~~~~~~~~~~~~~~
:Type: simple_dict_
:Required: No
:Default: do not coordinate
:Description: DESCRIPTION.

kgrid
~~~~~
:Type: int
:Required: Yes
:Default: None
:Description: Dimensions of the **k** point grid. Must be consistent with the
              **k** points specified in ``[kpoints]``.

``[opfm:qe]`` options
---------------------

outdir
~~~~~~
:Type: str
:Required: Yes
:Default: None
:Description: QuantumESPRESSO output directory

prefix
~~~~~~
:Type: str
:Required: Yes
:Default: None
:Description: QuantumESPRESSO run prefix

projwfc_out
~~~~~~~~~~~
:Type: str
:Required: No
:Default: None
:Description: QuantumESPRESSO projwfc.x output file

``[crystal]`` options
---------------------

The ``[crystal]`` section specifies the lattice vectors with the options ``a1``,
``a2``, and ``a3`` and the basis atoms with the ``basis`` option.

a1, a2, a3
~~~~~~~~~~
:Type: vector
:Required: Yes
:Default: None
:Description: Cartesian components of the lattice vectors.

basis
~~~~~
:Type: basis_
:Required: Yes
:Default: None
:Description: Atomic species and positions in crystal coordinates.

``[kpoints]`` options
---------------------

The ``[kpoints]`` section specifies the **k** points using the ``kpoints``
option.

kpoints
~~~~~~~
:Type: array_
:Required: Yes
:Default: None
:Description: The **k** points in crystal coordinates.


Example
=======

The format of these options can best be understood by examining the example
below. The ellipses for the ``basis`` and ``kpoints`` options indicate that
these options must have ``natms`` and ``nkpts`` lines, respectively.

::

   [opfm]
   site_basis = coordinate
   orbital_basis = psp
   kgrid:  nk1  nk2  nk3      # dimensions of the k point grid
   atomic_orbitals:  Si O     # atoms for which atomic orbitals are specified
     Si:  s;p                 # s and p (px, py, and pz) orbitals on Si
      O:  s;p                 # s and p (px, py, and pz) orbitals on O

   [crystal]
   a1:  a1_x  a1_y  a1_z      # Cartesian components of a1
   a2:  a2_x  a2_y  a2_z      # Cartesian components of a2
   a3:  a3_x  a3_y  a3_z      # Cartesian components of a3
   basis:  natms              # the number of basis atoms
     symbol1  tau1_b1  tau1_b2  tau1_b3      # crystal coordinates tau1 of atom with symbol1
     symbol2  tau2_b1  tau2_b2  tau2_b3      # crystal coordinates tau2 of atom with symbol2
         .        .        .        .
         .        .        .        .
         .        .        .        .

   [kpoints]
   kpoints:  nkpts                  # the number of k points in the grid
     kpt1_b1  kpt1_b2  kpt1_b3      # crystal coordinates of k point #1
     kpt2_b1  kpt2_b2  kpt2_b3      # crystal coordinates of k point #2
         .        .        .
         .        .        .
         .        .        .

   [opfm:qe]
   outdir: ../_scratch
   prefix: pwscf
   projwfc_out: projwfc.out
