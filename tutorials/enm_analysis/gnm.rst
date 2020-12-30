.. _gnm:

Gaussian Network Model (GNM)
===============================================================================

This example shows how to perform GNM calculations using an X-ray structure
of ubiquitin.  A :class:`.GNM` instance that stores the Kirchhoff matrix and
normal mode data describing the intrinsic dynamics of the protein structure
will be obtained.  :class:`.GNM` instances and individual normal modes
(:class:`.Mode`) can be used as input to functions in :mod:`prody.dynamics`
module.

See [Bahar97]_ and [Haliloglu97]_ for more information on the theory of GNM.

.. [Bahar97] Bahar I, Atilgan AR, Erman B. Direct evaluation of thermal
   fluctuations in protein using a single parameter harmonic potential.
   *Folding & Design* **1997** 2:173-181.

.. [Haliloglu97] Haliloglu T, Bahar I, Erman B. Gaussian dynamics of folded
   proteins. *Phys. Rev. Lett.* **1997** 79:3090-3093.


Parse structure
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion() # turn interactive mode on

First we parse a PDB file by passing its identifier to
:func:`.parsePDB` function. Note that if file is not found in
the current working directory, it will be downloaded.


.. ipython:: python

   ubi = parsePDB('1aar')
   ubi

This file contains 2 chains, and a flexible C-terminal (residues 71-76).
We only want to use Cα atoms of first 70 residues from chain A,
so we select them:

.. ipython:: python

   calphas = ubi.select('calpha and chain A and resnum < 71')
   calphas

See definition of "calpha", "chain", and other selection
keywords in :ref:`selections`.

Note that, flexible design of classes allows users to select atoms other than
alpha carbons to be used in GNM calculations.


Build Kirchhoff matrix
-------------------------------------------------------------------------------


Instantiate a :class:`.GNM` instance:

.. ipython:: python

   gnm = GNM('Ubiquitin')

We can build Kirchhoff matrix using selected atoms and
:meth:`.GNM.buildKirchhoff` method:

.. ipython:: python

   gnm.buildKirchhoff(calphas)


We can get a copy of the Kirchhoff matrix using :meth:`.GNM.getKirchhoff`
method:

.. ipython:: python

   gnm.getKirchhoff()


Parameters
-------------------------------------------------------------------------------

We didn't pass any parameters, but :meth:`.GNM.buildKirchhoff` method accepts
two of them, which by default are ``cutoff=10.0`` and ``gamma=1.0``, i.e.
``buildKirchhoff(calphas, cutoff=10., gamma=1.)``


.. ipython:: python

   gnm.getCutoff()
   gnm.getGamma()

Note that it is also possible to use an externally calculated Kirchhoff matrix.
Just pass it to the GNM instance using :meth:`.GNM.setKirchhoff` method.


Calculate normal modes
-------------------------------------------------------------------------------

We now calculate normal modes from the Kirchhoff matrix.

.. ipython:: python

   gnm.calcModes()

Note that by default 20 non-zero (or non-trivial) modes and 1 trivial mode are
calculated. Trivial modes are not retained. To calculate different numbers
of non-zero modes or to keep zero modes, try ``gnm.calcModes(50, zeros=True)``.


Normal mode data
-------------------------------------------------------------------------------

Get eigenvalues and eigenvectors:

.. ipython:: python

   gnm.getEigvals().round(3)
   gnm.getEigvecs().round(3)

Get covariance matrix:

.. ipython:: python

   gnm.getCovariance().round(2)

Note that covariance matrices are calculated using the available modes in the
model, which is the slowest 20 modes in this case. If the user calculates M
modes, these M modes will be used in calculating the covariance matrix.


Individual modes
-------------------------------------------------------------------------------

Normal mode indices start from 0, so slowest mode has index 0.

.. ipython:: python

   slowest_mode = gnm[0]
   slowest_mode.getEigval().round(3)
   slowest_mode.getEigvec().round(3)

By default, modes with 0 eigenvalue are excluded. If they were retained,
slowest non-trivial mode would have index 6.


Access hinge sites
-------------------------------------------------------------------------------
Hinge sites identified from all calculated modes (``n_modes`` defined when
calling ``gnm.calcModes()``) can be obtain by using the following command.

.. ipython:: python

   hinges = calcHinges(gnm)
   hinges[:5]

Hinge sites in the slowest mode can be obtained by:

.. ipython:: python

   calcHinges(gnm[0])

Hinge sites identified from multiple modes (e.g. 2 modes) can be accessed by:

.. ipython:: python

   calcHinges(gnm[:2])

These numbers correspond to node indices in the GNM object, which does not know 
anything about the original atoms. In order to get the residue numbers corresponding 
to these hinges, we can index the resum array with the hinges list as follows:

.. ipython:: python

   resnums = calphas.getResnums()
   mode2_hinges = calcHinges(gnm[1])
   resnums[mode2_hinges]


Plot results
-------------------------------------------------------------------------------


ProDy plotting functions are prefixed with ``show``. Let's use some of them
to plot data:

Contact Map
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. ipython:: python

   @savefig enm_analysis_gnm_contact_map.png width=4in
   showContactMap(gnm);

Cross-correlations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. ipython:: python

   @savefig enm_analysis_gnm_cross_corr.png width=4in
   showCrossCorr(gnm);

Slow mode shape
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
By default, hinge sites will be shown in mode shape plot indicated by
red stars, and it can be turned off by setting ``hinges=False``.
The option ``zero=True`` is to turn on the reference line of zero.

.. ipython:: python

   showMode(gnm[0], hinges=True, zero=True);
   @savefig enm_analysis_gnm_mode.png width=4in
   grid();

Square fluctuations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. ipython:: python

   @savefig enm_analysis_gnm_sqflucts.png width=4in
   showSqFlucts(gnm[0], hinges=True);

Protein structure bipartition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Given a GNM mode, protein structure can be partitioned into two parts that move 
with respect to each other. The function ``showProtein()`` can take a GNM mode 
as input and visualize the bipartition. 

.. ipython:: python

   @savefig enm_analysis_gnm_show_protein.png width=6in
   showProtein(calphas, mode=gnm[0]);
   plt.close('all')
   