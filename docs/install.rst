============
Installation
============

OPFm.py requires a working Python and the following additional packages:

- numpy
- scipy
- wannier90-utils
- codiag
- f90nml

OPFm.py has been tested with the following versions of numpy and scipy:

===== ======= =======
----- ------- -------
numpy 1.8.2   1.9.3
scipy 0.15.0  0.16.0
===== ======= =======

Assuming all the requirements have been met, the package can be installed as
follows:

Download the source and unpack in a convenient directory::

   tar -zxvf opfm-|version|.tar.gz
   cd opfm-|version|

It is best to first build the package locally, in order to run tests and make
sure the the extensions are working properly, before installing. Build the
package in the current directory::

   python setup.py build

The package can be tested using ``pytest``::

   py.test build/

Once the build is verified, installation can be completed by running::

   python setup.py install

Python
======

Unix/Linux or OS X
------------------
Python is available on almost all Unix/Linux systems. To check if Python is
already installed on your system, execute the following in the terminal::

   python -V

The package has been tested Python version 2.7.x, while support for
Python 3.x is planned for the near future.

Windows
-------

Pyhton can be obtained for Windows through a variety releases or as part of a
`SciPy Stack`_, which will include the required NumPy and SciPy packages.


NumPy and SciPy
===============


Additional software
===================

Optional, but useful, packages needed for features of the package are listed
below:

- IPython
- matplotlib


.. _SciPy Stack: http://www.scipy.org/stackspec.html
