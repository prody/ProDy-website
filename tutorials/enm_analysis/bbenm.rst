
Bond-Bending Elastic Network Model (bbENM)
===============================================================================

This example shows how to perform bond bending ENM calculations, and retrieve
normal mode data.  An :class:`.bbENM` instance that stores Hessian matrix
(and also Kirchhoff matrix) and normal mode data describing the intrinsic
dynamics of the protein structure will be obtained.  :class:`.bbENM` instances
and individual normal modes (:class:`.Mode`) can be used as input to functions
in :mod:`.dynamics` module.

See [Srivastava12]_ for more information on the theory of bbENM.

.. [Srivastava12] Srivastava A, Halevi RB, Veksler A, Granek R. 
Tensorial elastic network model for protein dynamics: Integration of
the anisotropic network model with bond-bending and twist elasticities. *Proteins* **2012** 80:2692-2700.


Parse structure
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

We start with parsing a PDB file by passing an identifier.
Note that if a file is not found in the current working directory, it will be
downloaded.

.. ipython:: python

   p38 = parsePDB('1p38')
   p38

We want to use only Cα atoms, so we select them:

.. ipython:: python

   calphas = p38.select('protein and name CA')
   calphas

Note that, ProDy atom selector gives the flexibility to select any set of atoms
to be used in ANM  calculations.

Build Hessian
-------------------------------------------------------------------------------

We instantiate an :class:`.bbENM` instance:

.. ipython:: python

   bbenm = bbENM('p38 bond-bending ENM analysis')

Then, build the Hessian matrix by passing selected atoms (351 Cα's)
to :meth:`.bbENM.buildHessian` method:

.. ipython:: python

   bbenm.buildHessian(calphas)

We can get a copy of the Hessian matrix using :meth:`.bbENM.getHessian` method:

.. ipython:: python

   bbenm.getHessian().round(3)


Parameters
-------------------------------------------------------------------------------

We didn't pass any parameters to :meth:`.ANM.buildHessian` method, but it
accepts *cutoff* and *gamma* parameters, for which  default values are
``cutoff=15.0`` and ``gamma=1.0``.


Calculate normal modes
-------------------------------------------------------------------------------

Calculate modes using :meth:`.ANM.calcModes` method since bbENM object is a form of an ANM object:

.. ipython:: python

   bbenm.calcModes()

Note that by default 20 non-zero (or non-trivial) and 6 trivial modes are
calculated. Trivial modes are not retained. To calculate a different number
of non-zero modes or to keep zero modes, try ``anm.calcModes(50, zeros=True)``.

Normal modes data
-------------------------------------------------------------------------------

.. ipython:: python

   bbenm.getEigvals().round(3)
   bbenm.getEigvecs().round(3)


You can get the covariance matrix as follows:

.. ipython:: python

   bbenm.getCovariance().round(2)

Covariance matrices are calculated using the available modes (slowest 20 modes
in this case). If the user calculates M slowest modes, only they will be used
in the calculation of covariances.

Individual modes
-------------------------------------------------------------------------------

Normal mode indices in Python start from 0, so the slowest mode has index 0.
By default, modes with zero eigenvalues are excluded. If they were retained,
the slowest non-trivial mode would have index 6.

Get the slowest mode by indexing :class:`.bbENM` instance as follows:

.. ipython:: python

   slowest_mode = bbenm[0]
   slowest_mode.getEigval().round(3)
   slowest_mode.getEigvec().round(3)


Write NMD file
-------------------------------------------------------------------------------

ANM results in NMD format can be visualized using :ref:`nmwiz` VMD_ plugin.
The following statement writes the slowest 3 ANM modes into an NMD file:

.. ipython:: python

   writeNMD('p38_bbenm_modes.nmd', bbenm[:3], calphas)


Note that slicing an :class:`.bbENM` objects returns a list of modes.
In this case, slowest 3 bbENM modes were written into NMD file.

View modes in VMD
-------------------------------------------------------------------------------

First make sure that the VMD path is correct

.. ipython:: python

   pathVMD()


.. ipython:: python
   :verbatim:

   # if this is incorrect use setVMDpath to correct it
   viewNMDinVMD('p38_bbenm_modes.nmd')

This will show the slowest 3 modes in VMD using NMWiz. This concludes the bbENM
example. 

