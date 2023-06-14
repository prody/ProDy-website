.. _insty_tutorial:

Interactions/Stability Evaluation (InSty)
===============================================================================

This example shows how to perform Interactions/Stability Evaluation
(**InSty**) analysis for a small protein (<200 residues) called tyrosine
phosphatase LMW-PTP (**5KQM**) and visualize the results using Matplotlib_
library and VMD_ program. 
In the tutorial, we will use already preapared structure for
simulation (with hydrogens added). The same structure will be later
analyzed with the trajectory file to show how the analysis of interactions 
in the course of simulation can change. 

The tutorial will also include an example of a PDB structure directly
downloaded from Protein Data Bank (PDB) which requires adding the missing hydrogen
atoms to the protein and ligand structure. The analysis will be performed for
protein-ligand interactions.


Analysis of interactions for a single PDB structure
-------------------------------------------------------------------------------

We start by parsing PDB file with LMW-PTP **5kqm_all_sci.pdb** which is avalable
as the tutorial files. PDB file contains protein structures with water and 
counter ions prepared using VMD_ program.

Before that import everything from the ProDy packages.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on


.. ipython:: python

   PDBfile = '5kqm_all_sci.pdb'
   coords = parsePDB(PDBfile)
   coords

For the analysis we will use only protein coordinates (*atoms*):

.. ipython:: python

   atoms = coords.select('protein')
   atoms


Compute all types of interactions
-------------------------------------------------------------------------------

In the next step, we instantiate an :class:`.Interactions` instance:

.. ipython:: python

   interactions = Interactions()


Now we can compute all available types of interactions (seven types: hydrogen
bonds, salt bridges, repulsive ionic bonding, Pi-cation, Pi-stacking,
hydrphobic interactions, and disulfide bonds) for protein structure by passing
selected atoms (*atoms*) to :meth:`.Interactions.calcProteinInteractions` method:

.. ipython:: python

   all_interactions = interactions.calcProteinInteractions(atoms)

All types of interactions will be displayed on the screen with all types of
information such as distance or angle (if applied).

Moreover, we will have access to the details of each interaction type
using the following methods: 

:meth:`.Interactions.getHydrogenBonds` - hydrogen bonds:

.. ipython:: python
   
   interactions.getHydrogenBonds()


:meth:`.Interactions.getSaltBridges` - salt bridges (residues with oposite
charges):

.. ipython:: python
   
   interactions.getSaltBridges()


:meth:`.Interactions.getRepulsiveIonicBonding` - repulsive ionic bonding
(between residues with the same charges):

.. ipython:: python

   interactions.getRepulsiveIonicBonding()


:meth:`.Interactions.getPiStacking` - Pi-stacking interactions:

.. ipython:: python

   interactions.getPiStacking()


:meth:`.Interactions.getPiCation` - Pi-cation:

.. ipython:: python

   interactions.getPiCation()


:meth:`.Interactions.getHydrophohic` - hydrophobic interactions:

.. ipython:: python

   interactions.getHydrophohic()


:meth:`.Interactions.getDisulfideBonds` - disulfide bonds:

.. ipython:: python

   interactions.getDisulfideBonds()


To display residues with the biggest number of potential interactions and their
types, we can use :meth:`.Interactions.getFrequentInteractions` method:

.. ipython:: python

   frequent_interactions = interactions.getFrequentInteractions(contacts_min=3)
   frequent_interactions

The value of *contacts_min* can be modified to display residues with smaller
number of interactions. 


Visualize interactions in VMD
-------------------------------------------------------------------------------

We can generate tcl files for visualizing each type of interaction with VMD_ 
using the :func:`.showProteinInteractions_VMD` function in the following way:

.. ipython:: python

   showProteinInteractions_VMD(atoms, interactions.getHydrogenBonds(), color='blue', filename='HBs.tcl')
   showProteinInteractions_VMD(atoms, interactions.getSaltBridges(), color='yellow',filename='SBs.tcl')
   showProteinInteractions_VMD(atoms, interactions.getRepulsiveIonicBonding(), color='red',filename='RIB.tcl')
   showProteinInteractions_VMD(atoms, interactions.getPiStacking(), color='green',filename='PiStacking.tcl') 
   showProteinInteractions_VMD(atoms, interactions.getPiCation(), color='orange',filename='PiCation.tcl') 
   showProteinInteractions_VMD(atoms, interactions.getHydrophobic(), color='silver',filename='HPh.tcl')
   showProteinInteractions_VMD(atoms, interactions.getDisulfideBonds(), color='black',filename='DiBs.tcl') 

A TCL file will be saved and can be used in VMD_ after uploading the PDB file
with protein structure **5kqm_all_sci.pdb** and by running the following command 
line instruction in the VMD_ *TKConsole* (*VMD Main*) for Linux, Windows and Mac users: 

::  play HBs.tcl

The tcl file contains a method for drawing lines between selected pairs of 
residues. Those residues are also displayed.

.. figure:: images/HBs.tga
   :scale: 60 %


::  play SBs.tcl

.. figure:: images/SBs.tga
   :scale: 60 %


::  play RIB.tcl

.. figure:: images/RIB.tga
   :scale: 60 %


::  play PiStacking.tcl

.. figure:: images/PiStacking.tga
   :scale: 60 %


::  play PiCation.tcl

.. figure:: images/PiCation.tga
   :scale: 60 %


::  play HPh.tcl

.. figure:: images/Hydrophobic.tga
   :scale: 60 %



Additional selections
-------------------------------------------------------------------------------

From the predicted interactions we can select only interactions assigned to
a certain regions, chains or between different chains.

We can compute them by adding additional parameters to the selected
function. See examples below:

.. ipython:: python

   interactions.getSaltBridges(selection='chain P')


.. ipython:: python

   interactions.getRepulsiveIonicBonding(selection='resid 102')


.. ipython:: python

   interactions.getPiStacking(selection='chain P and resid 26')


It can be done for all kinds of interactions as well. The function will
return a list of interactions with following order:

    (1) Hydrogen bonds
    (2) Salt Bridges
    (3) RepulsiveIonicBonding 
    (4) Pi stacking interactions
    (5) Pi-cation interactions
    (6) Hydrophobic interactions
    (7) Disulfide bonds

.. ipython:: python

   allRes_20to50 = interactions.getInteractions(selection='resid 20 to 50')
   allRes_20to50


The list of hydrogen bonds, salt bridges and other types of interactions can
be displayed as follows:

.. ipython:: python

   allRes_20to50[0]


Salt Bridges:

.. ipython:: python

   allRes_20to50[1]


We can also select one particular residue of our interest:

.. ipython:: python

   interactions.getPiCation(selection='resid 85')


.. ipython:: python

   interactions.getHydrophobic(selection='resid 26 to 100')


Change selection criteria for interaction type
-------------------------------------------------------------------------------

The :meth:`.Interactions.buildInteractionMatrix` method computes interactions 
using default parameters for interactions. However, it can be changed
according to our needs. To do that, we need to recalculate the selected type
of interactions. 

We can do it using the following functions: :func:`.calcHydrogenBonds`,
:func:`.calcHydrogenBonds`, :func:`.calcSaltBridges`,
:func:`.calcRepulsiveIonicBonding`, :func:`.calcPiStacking`,
:func:`.calcPiCation`, :func:`.calcHydrophohic`,
:func:`.calcDisulfideBonds`, and use
:meth:`.Interactions.setNewHydrogenBonds`,
:meth:`.Interactions.setNewSaltBridges`,
:meth:`.Interactions.setNewRepulsiveIonicBonding`,
:meth:`.Interactions.setNewPiStacking`,
:meth:`.Interactions.setNewPiCation`,
:meth:`.Interactions.setNewHydrophohic`,
:meth:`.Interactions.setNewDisulfideBonds` method to replace it in the main
Instance. 

For example:

.. ipython:: python

   newHydrogenBonds2 = calcHydrogenBonds(atoms, distA=2.8, angle=30, cutoff_dist=15)
   interactions.setNewHydrogenBonds(newHydrogenBonds2)
   
.. ipython:: python

   interactions.getHydrogenBonds()

.. ipython:: python

   sb2 = calcSaltBridges(atoms, distA=6)
   interactions.setNewSaltBridges(sb2)

.. ipython:: python

   rib2 = calcRepulsiveIonicBonding(atoms, distA=9)
   interactions.setNewRepulsiveIonicBonding(rib2)

.. ipython:: python

   picat2 = calcPiCation(atoms, distA=7)
   interactions.setNewPiCation(picat2)



Assess the functional significance of a residue
-------------------------------------------------------------------------------

For assessing the functional significance of each residue in protein
structure, we counted the number of possible contacts based on:

    (1) Hydrogen bonds (HBs)
    (2) Salt Bridges (SBs)
    (3) Repulsive Ionic Bonding (RIB)  
    (4) Pi stacking interactions (PiStack)
    (5) Pi-cation interactions (PiCat) 
    (6) Hydrophobic interactions (HPh) 
    (7) Disulfide Bonds (DiBs)


To compute the weighted interactions use the 
:meth:`.Interactions.buildInteractionMatrix` method:

.. ipython:: python

   matrix = interactions.buildInteractionMatrix()


The results can be displayed in the following way:

.. ipython:: python

    import matplotlib.pylab as plt
    showAtomicMatrix(matrix, atoms=atoms.ca, cmap='seismic', markersize=8)
    plt.xlabel('Residue')
    plt.ylabel('Residue')
    plt.clim([-3,3])


The total number of interaction for each residue can be displayed on the plot using
:func:`.showCumulativeInteractionTypes()` function.

.. ipython:: python

   interactions.showCumulativeInteractionTypes()


The results with the higest number of possible contacts can be saved in PDB
file. They will be restored in Occupancy column and display in VMD_.

.. ipython:: python

   interactions.saveInteractionsPDB(filename='5kqm_meanMatrix.pdb')



Visualize number of interactions onto 3D structure
-------------------------------------------------------------------------------

The number of the interaction can be saved to a PDB file in the
*Occupancy* column by using :meth:`.Interactions.saveInteractionsPDB`
method. Then the score would be displayed in color in any available graphical
program, for example, in VMD_.
 
.. ipython:: python

   interactions.saveInteractionsPDB(filename='5kqm_meanMatrix.pdb')


A file *5kqm_meanMatrix.pdb* will be saved and can be used in VMD_ by 
uploading PDB structure and displaying it with *Coloring Method*
*Occupancy*. By default blue colors correspond to the highest values but we
can change it in *VMD Main* -> *Graphics* -> *Color Controls* -> *Color
Scale* -> *Method* to *BWR*.

.. figure:: images/fig1.tga
   :scale: 60 %


Exclude some interaction types from calculations
-------------------------------------------------------------------------------

For analysis we can exclude some of the interaction types by assigning zero
to the type of interactions (HBs - hydrogen bonds, SBs - salt bridges, RIB -
repulsive ionic bonding, PiCat - Pi-Cation, PiStack - Pi-Stacking, HPh -
hydrophobic interactions and finally DiBs - disulfide bonds). 

.. ipython:: python

   matrix = interactions.buildInteractionMatrix(RIB=0, HBs=0, HPh=0, DiBs=0)


The results can be displayed in a similar way:

.. ipython:: python

    showAtomicMatrix(matrix, atoms=atoms.ca, cmap='seismic', markersize=8)
    plt.xlabel('Residue')
    plt.ylabel('Residue')
    plt.clim([-3,3])
