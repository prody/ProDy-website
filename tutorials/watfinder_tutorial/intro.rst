Introduction
===============================================================================

This tutorial shows how to detect water molecules that might form hydrogen 
bonds with protein structures (called water bridges). The prediction method 
introduced here helps evaluate the significance of water molecules on the 
stability and dynamics of protein structure. 


Required Programs
-------------------------------------------------------------------------------

The latest version of ProDy_ is required.


Recommended Programs
-------------------------------------------------------------------------------

Besides ProDy_, the Matplotlib_ library and VMD_ program are required for
some steps in the tutorial. IPython_ is highly recommended for interactive usage.

Moreover, in the case of the lack of hydrogen atoms in protein structures,
additional packages such as Openbabel_ or PDBfixer_ are required for
predicting hydrogen bonds.

.. _Openbabel: https://github.com/openbabel
.. _PDBfixer: https://github.com/openmm/pdbfixer

They can be installed using conda or pip.

Getting Started
-------------------------------------------------------------------------------

To follow this tutorial, you will need the following files:

.. files.txt will be automatically generated

.. literalinclude:: files.txt


We recommend that you will follow this tutorial by typing commands in an
IPython session, e.g.::

  $ ipython

or with pylab environment::

  $ ipython --pylab


First, we will make necessary imports from ProDy and Matplotlib packages.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.
