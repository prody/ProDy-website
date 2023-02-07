.. _interactions_analysis:

Analysis of interactions for a single PDB structure without hydrogens
===============================================================================

Very offten PDB structures downloaded directly from PDB database will not have
determined hydrogen atoms which are required, for example, for predicting
hydrogen bonds. In such case, we need can use :func:`.addHydrogens`
function. It will allow us to use one of two available methods (*openbabel*
or *pdbfixer*) to predict the position of hydrogen atoms in protein/ligand structure.

To use this function we need to install additional Python package(s).

Instalation of Openbabel:

:: conda install -c conda-forge openbabel   

Instalation of PDBfixer:

:: conda install -c conda-forge pdbfixer


Add missing hydrogen atoms to protein structure
-------------------------------------------------------------------------------

We start as usual by importing everything from the ProDy package. It is not
required when you already did that in the previous part of the tutorial.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on

We start by parsing the PDB file with LMW-PTP **5kqm.pdb**. Openbabel
requires having the PDB file in the same folder. Therefore, it needs to be 
downloaded and saved to successfully perform the operation with adding 
missing hydrogens. New file will be saved with the same name with additional
prefix 'addH_'.

.. ipython:: python

   PDBname2 = '5kqm.pdb'
   addHydrogens(PDBname2, method='openbabel')

To used PDBfixer:

.. ipython:: python
   addHydrogens(PDBname2, method='pdbfixer')


Next, we can parse saved structure with hydrogen atoms to ProDy and analyze
it in the same way as in the previous paragraph.

.. ipython:: python

   pdb2 = parsePDB('addH_'+str(PDBname2)).select('protein')

 
