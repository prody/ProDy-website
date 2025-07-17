.. _esty_tutorial:

Protein-ligand interactions
===============================================================================

With the additional installation of the PLIP package, the analysis of
protein-ligand interactions is also possible in ProDy.

Example installation of PLIP (for Anaconda users) is the following:

:: conda install -c conda-forge plip


Start by importing everything from the ProDy packages unless you
already did that.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on


Example 1
-------------------------------------------------------------------------------

Parse protein and ligand structure from *5kqm_all_sci.pdb* which we used
before. Selected elements, such as protein and ligand MES, will be saved in
a separate PDB file in the local directory.

.. ipython:: python

   atoms2 = parsePDB('5kqm_all_sci.pdb').select('protein or resname MES')


To compute the interactions use :func:`.calcLigandInteractions` function:

.. ipython:: python

   ligands_interactions, ligands = calcLigandInteractions(atoms2)


The type of the ligand can be checked as follows:

.. ipython:: python

   ligands


To display the interactions between protein and ligand MES use
func:`.listLigandInteractions` function:

.. ipython:: python

   ligandMES_interactions = ligands_interactions[0]
   protein_ligand_interactions = listLigandInteractions(ligandMES_interactions)

We can display them: 

.. ipython:: python

   for nr_k,k in enumerate(protein_ligand_interactions):
       print ("%3i%24s%10s%28s%4s    <---> %8s%12s%4s%6.1f" % (nr_k,k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]))


Interactions can be saved for visualization in VMD_. We can use
func:`.showLigandInteraction_VMD` function in the following way:

.. ipython:: python

   showLigandInteraction_VMD(atoms2,protein_ligand_interactions)


A TCL file will be saved and can be used in VMD_ after uploading PDB file
with protein and ligand structure **5kqm_all_sci_sele.pdb** and by running
in the command line instruction in the VMD_ *TKConsole* (*VMD Main*): 

::  play 5kqm_all_sci_interaction.tcl


The tcl file contains a method for drawing lines between protein and ligand
(hydrogen bonds - *blue*, salt bridges - *yellow*, pi-stacking - *green*, cation-pi -
*orange*, hydrophobic - *silver*, and water bridges - *cyan*).

.. figure:: images/lig1.png



Example 2
-------------------------------------------------------------------------------

Another example is with protein-ligand structure downloaded from PDB. This
structure doesn't have hydrogen atoms which will be added using *openbabel*
method described in the previous paragraph. 

First, we need to add missing hydrogens. Openbabel will save a new structure
under a similar name with 'addH_' prefix and '_sele.pdb' suffix.

.. ipython:: python

   PDBname3 = '3ugc.pdb'
   addHydrogens(PDBname3, method='openbabel')
   pdb3 = parsePDB('addH_'+str(PDBname3[:-4])+'_sele.pdb')


To select protein and ligand structures for analysis and compute
interactions use the following functions. The procedure is similar to
Example 1.

.. ipython:: python

   atoms3 = pdb3.select('protein or resname 046')
   ligands_interactions3, ligands3 = calcLigandInteractions(atoms3)


.. ipython:: python

   ligand046_interactions = ligands_interactions3[0]
   protein_ligand_interactions3 = listLigandInteractions(ligand046_interactions)


.. ipython:: python

   for nr_k,k in enumerate(protein_ligand_interactions3):
       print ("%3i%24s%10s%28s%4s    <---> %8s%12s%4s%6.1f" % (nr_k,k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]))


Interactions can be saved for visualization in VMD_:

.. ipython:: python

   showLigandInteraction_VMD(atoms3,protein_ligand_interactions3)


A tcl file will be saved and can be used in VMD_ after uploading the PDB file
with protein and ligand structure **addH_3ugc_sele.pdb** and by running
in the command line instruction in the VMD_ *TKConsole* (*VMD Main*): 

   ::  play addH_3ugc_sele_interaction.tcl


The tcl file contains a method for drawing lines between protein and ligand
(hydrogen bonds - *blue*, salt bridges - *yellow*, pi-stacking - *green*, cation-pi -
*orange*, hydrophobic - *silver*, and water bridges - *cyan*).

.. figure:: images/lig2.png

