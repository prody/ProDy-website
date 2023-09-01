.. _wbfinder_tutorial:

Protein Preparation
===============================================================================

Because in PDB structures we will have water molecules without 
hydrogens we would have to add them. We can use :func:`.addMissingAtoms` 
function. This function use either *Openbabel* or *PDBFixer* and both 
packages required PDB file saved in the local directory. Those are 
external packages therefore they should be installed 
independently. 

How to install them can be found in :class:`.Interactions`.
We can also directly provide PDB structure with hydrogens added by other 
software.

Here we will fetch PDB structure with LMW-PTP **5kqm** from 
Protein Data Bank (PDB) in incompressed form using *compressed=False* 
and add missing atoms using :func:`.addMissingAtoms`:


.. ipython:: python

   PDB = '5kqm'
   fetchPDB(PDB, compressed=False)


.. ipython:: python

   addMissingAtoms(PDB+'.pdb')


A new file with hydrogens is now created. It will contain *'addH_'* as 
a prefix. For protein structures that don't contain hydrogens, 
results will also be computed without applying angle criteria.


