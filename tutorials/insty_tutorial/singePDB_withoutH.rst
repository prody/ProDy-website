.. _insty_tutorial:
=======

PDB structure without hydrogens
===============================================================================

Very often, PDB structures downloaded directly from the PDB database will not
have determined hydrogen atoms that are required, for example, for predicting
hydrogen bonds. In such a case, we can use the :func:`.addHydrogens` function.
It will allow us to use one of two available methods (*openbabel* or *pdbfixer*)
to predict the position of hydrogen atoms in protein structure.

To use one of those functions, we need to install additional Python package(s).
For Anaconda users, the installation will be the following:

Installation of Openbabel_:

:: conda install -c conda-forge openbabel   

Installation of PDBfixer_:

:: conda install -c conda-forge pdbfixer


.. _Openbabel: https://github.com/openbabel
.. _PDBfixer: https://github.com/openmm/pdbfixer

Add missing hydrogen atoms to the structure
-------------------------------------------------------------------------------

We start by fetching the PDB file with **5KQM** code (*5kqm.pdb*).
Openbabel_ requires having the PDB file in the same folder. Therefore, it
needs to be downloaded and saved to successfully perform the operation with
adding missing hydrogens. A new file will be saved with the same name with
the additional prefix 'addH_'.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on

Openbabel_ or PDBfixer_ require PDB file saved in the direcory. Therefore
first it needs to be downloaded.

.. ipython:: python
   :verbatim:

   fetchPDB('5kqm', compressed=False)

.. parsed-literal::

   '5kqm.pdb'

When PDB file is already in the local directory, we can choose between
Openbabel_ and PDBfixer_ to add missing hydrogen bonds to the protein
structure:

Openbabel:

.. ipython:: python
   :verbatim:

   PDBname = '5kqm.pdb'
   addMissingAtoms(PDBname, method='openbabel')

.. parsed-literal::

   @> Hydrogens were added to the structure. Structure addH_5kqm.pdb is saved in the local directry.

PDBfixer:

.. ipython:: python
   :verbatim:

   addMissingAtoms(PDBname, method='pdbfixer')

.. parsed-literal::

   @> Hydrogens were added to the structure. New structure is saved as addH_5kqm.pdb.

Next, we can parse the saved structure with hydrogen atoms to ProDy and analyze
it in the same way as in the previous paragraph.

.. ipython:: python
   :verbatim:

   atoms = parsePDB('addH_'+str(PDBname)).select('protein')

.. parsed-literal::

   @> 2800 atoms and 1 coordinate set(s) were parsed in 0.03s.


