=====
Usage
=====

Workflow
========

#. Perform DFT calculation (only QuantumESPRESSO supported at this point):

   i. compute the self-consistent charge density::

         pw.x < scf.in > scf.out

   ii. compute wavefunctions and energies on a uniform **k** point grid::

         pw.x < bands.in > bands.out

#. Create an ``opfm.cfg`` file from the output::

      opfm.py create-cfg

#. Setup the calculation of the projections for the specified :ref:`atomic orbitals <atomic_orbitals>`::

      opfm.py setup

   This creates a "hack" Wannier90 input file, which allows the generation of an
   NNKP file.

#. Run Wannier90 in post-processing mode::

      wannier90.x -pp opfm

#. Compute the AMN, MMN, and EIG files::

      pw2wannier90.x < pw2wan.in > pw2wan.out

#. Finally, construct the optimized projection functions::

      opfm.py optimize

#. Optionally, run Wannier90 to obtain maximally localized Wannier functions::

      wannier90.x


In most cases, computing the projections needs to be done only once, since we
have projected onto a large set of atomic orbitals for each atom in the basis.

.. autoprogram:: opfmpy.scripts.opfm:parser
   :prog: opfm.py
