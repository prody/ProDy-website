<<<<<<< HEAD:tutorials/insty_tutorial/singePDB_withoutH.rst
.. _insty_tutorial:
=======
.. _esty_tutorial:
>>>>>>> 808bf0122eac0c85aa7179d5eaf18b85b7638e3f:tutorials/esty_tutorial/singePDB_withoutH.rst

PDB structure without hydrogens
===============================================================================

Very often PDB structures downloaded directly from the PDB database will not
have determined hydrogen atoms that are required, for example, for predicting
hydrogen bonds. In such case, we need can use :func:`.addHydrogens`
function. It will allow us to use one of two available methods (*openbabel*
or *pdbfixer*) to predict the position of hydrogen atoms in protein/ligand structure.

To use one of those functions, we need to install additional Python package(s).
For Anaconda users the installation will be the following:

Installation of Openbabel:
:: conda install -c conda-forge openbabel   

Installation of PDBfixer:
:: conda install -c conda-forge pdbfixer


Add missing hydrogen atoms to the structure
-------------------------------------------------------------------------------

We start by fetching the PDB file with **5KQM** code (*5kqm.pdb*). Openbabel
requires having the PDB file in the same folder. Therefore, it needs to be 
downloaded and saved to successfully perform the operation with adding 
missing hydrogens. A new file will be saved with the same name with the
additional prefix 'addH_'.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on

Openbable or PDBFixer require PDB file saved in the direcory. Therefore
first it needs to be downloaded.

.. ipython:: python

   fetchPDB('5kqm', compressed=False)

When PDB file is already in the local directory we can choose between
Openbabel and PDBFixer to add missing hydrogen bonds to the protein
structure:

Openbabel:

.. ipython:: python

   PDBname = '5kqm.pdb'
   addMissingAtoms(PDBname, method='openbabel')

PDBfixer:

.. ipython:: python
   addMissingAtoms(PDBname, method='pdbfixer')


Next, we can parse the saved structure with hydrogen atoms to ProDy and analyze
it in the same way as in the previous paragraph.

.. ipython:: python

   pdb2 = parsePDB('addH_'+str(PDBname)).select('protein')




