.. _esty_tutorial:

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

Before that import everything from the ProDy packages unless you already did that.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on


.. ipython:: python

   PDBname2 = '5kqm.pdb'
   addHydrogens(PDBname2, method='openbabel')

Instead Openbabel we can use PDBfixer:

.. ipython:: python
   addHydrogens(PDBname2, method='pdbfixer')


Next, we can parse the saved structure with hydrogen atoms to ProDy and analyze
it in the same way as in the previous paragraph.

.. ipython:: python

   pdb2 = parsePDB('addH_'+str(PDBname2)).select('protein')

 
