.. _stabilee_tutorial:

Stability/Energetic Evaluation
===============================================================================

This example shows how to perform Stability/Energetic Evaluation (StabilEE)
analysis for a small protein (<200 residues) called tyrosine phosphatase 
LMW-PTP (**5KQM**) and visualize the results using Matplotlib_ library and 
VMD_ program. 
In the tutorial, we will use already preapared structure for
simulation (with hydrogens added). The same structure will be later
analyzed with the trajectory file to show the analysis of interactions 
in the course of simulation. 

The tutorial will also include an example of a PDB structure directly
downloaded from Protein Data Bank which requires adding the missing hydrogen
atoms to the protein and ligand structure. The analysis will be performed for
protein-ligand interactions.


Parse structure
-------------------------------------------------------------------------------

We start by parsing PDB file with LMW-PTP **5kqm_all_sci.pdb** which is avalable
as the tutorial files. PDB file contains protein structures with water and 
counter ions prepared using VMD_ program.

Before that import everything from the ProDy packages unless you already did that.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on


.. ipython:: python

   PDBfile = '5kqm_all_sci.pdb'
   coords = parsePDB(PDBfile)
   coords

For the analysis we will use only protein structure (*atoms*):

.. ipython:: python

   atoms = coords.select('protein')
   atoms


Compute all types of interactions
-------------------------------------------------------------------------------

In the next step, we instantiate an :class:`.Interactions` instance:

.. ipython:: python

   interactions = Interactions('5kqm')


Now we can compute all available types of interactions (six types: hydrogen
bonds, salt bridges, repulsive ionic bonding, Pi-cation, Pi-stacking, and
hydrphobic) for protein structure by passing selected atoms (*atoms*) to
:meth:`.Interactions.calcProteinInteractions` method:

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


To display residues with the biggest number of interactions and their type, we
can use :meth:`.Interactions.getFrequentInteractions` method:

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

   showProteinInteractions_VMD(atoms, interactions.getHydrogenBonds(), color='blue', output='HBs.tcl')
   showProteinInteractions_VMD(atoms, interactions.getSaltBridges(), color='yellow',output='SBs.tcl')
   showProteinInteractions_VMD(atoms, interactions.getRepulsiveIonicBonding(), color='red',output='RIB.tcl')
   showProteinInteractions_VMD(atoms, interactions.getPiStacking(), color='green',output='PiStacking.tcl') 
   showProteinInteractions_VMD(atoms, interactions.getPiCation(), color='orange',output='PiCation.tcl') 
   showProteinInteractions_VMD(atoms, interactions.getHydrophohic(), color='silver',output='HPh.tcl')


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



Assess the functional significance of a residue
-------------------------------------------------------------------------------

As a criterion for assessing the functional significance of a residue, we 
included weights of interaction type in the algorithm. The occurrence 
(and relative strength) of different types of interactions might suggest 
the importance of the region in protein structure.

We will use by default the following scoring: 

    (1) Hydrogen bonds (HBs) +2
    (2) Salt Bridges (SBs) +3
    (3) Repulsive Ionic Bonding (RIB) -1 
    (4) Pi stacking interactions (PiStack) +3
    (5) Pi-cation interactions (PiCat) +3
    (6) Hydrophobic interactions (HPh) +1


To compute the weighted interactions use the 
:meth:`.Interactions.buildInteractionMatrix` method:

.. ipython:: python

   matrix = interactions.buildInteractionMatrix()


The results can be displayed in the following way:

.. ipython:: python

   import matplotlib.pylab as plt
   plt.imshow(matrix, interpolation='none', cmap='seismic')
   plt.clim([-3,3])
   plt.xlabel('Residue')
   plt.ylabel('Residue')
   plt.colorbar()
   plt.tight_layout()


Mean value of interaction for each residue can be displayed on the plot using
:func:`.showInteractions` function.

.. ipython:: python

   interactions.showInteractions()


Residues with the higest score can be displayed using 
:meth:`.Interactions.showFrequenctInteractions` method.

.. ipython:: python

   interactions.showFrequenctInteractions()

We can change the minimum value of score using *cutoff* option:

.. ipython:: python

   interactions.showFrequenctInteractions(cutoff=3)


Visualize weighted interactions onto 3D structure
-------------------------------------------------------------------------------

The mean value of the interaction map can be saved to a PDB file in the
*Occupancy* column by using :meth:`.Interactions.saveInteractionsPDB`
method. Then the score would be displayed in color in any available graphical
program, for example, in VMD_.
 
.. ipython:: python

   interactions.saveInteractionsPDB(output='5kqm_meanMatrix.pdb')


A file *5kqm_meanMatrix.pdb* will be saved and can be used in VMD_ by 
uploading PDB structure and displaying it with *Coloring Method*
*Occupancy*. By default blue colors correspond to the highest values but we
can change it in *VMD Main* -> *Graphics* -> *Color Controls* -> *Color
Scale* -> *Method* to *BWR*.

.. figure:: images/fig1.tga
   :scale: 60 %


Change weights for interaction types
-------------------------------------------------------------------------------

The default weights for interaction types can be easily changed by using
keywords: *HBs*, *SBs*, *RIB*, *PiStack*, *PiCat*, *HPh* in the following way: 

.. ipython:: python

   matrix = interactions.buildInteractionMatrix(HBs=3, SBs=4)


In such a way, we can exclude some types of interactions:

.. ipython:: python

   matrix = interactions.buildInteractionMatrix(RIB=0, HBs=0, HPh=0) 


.. ipython:: python

   interactions.showInteractions()


.. ipython:: python

   interactions.showFrequenctInteractions()


The number of interactions for each residue in the protein structure can be
checked by using a score equal to 1 for each interaction type:

.. ipython:: python

   matrix = interactions.buildInteractionMatrix(RIB=1, PiStack=1, PiCat=1, HBs=1, HPh=1, SBs=1)


.. ipython:: python

   interactions.showInteractions()


.. ipython:: python

   interactions.showFrequenctInteractions()



Change selection criteria for interaction type
-------------------------------------------------------------------------------

The :meth:`.Interactions.buildInteractionMatrix` method computes interactions 
using default parameters for interactions. However, it can be changed
according to our needs. To do that, we need to recalculate the selected type
of interactions. 

We can do it using the following functions: :func:`.calcHydrogenBonds`,
:func:`.calcHydrogenBonds`, :func:`.calcSaltBridges`,
:func:`.calcRepulsiveIonicBonding`, :func:`.calcPiStacking`,
:func:`.calcPiCation`, :func:`.calcHydrophohic`, and use
:meth:`.Interactions.setNewHydrogenBonds`,
:meth:`.Interactions.setNewSaltBridges`,
:meth:`.Interactions.setNewRepulsiveIonicBonding`,
:meth:`.Interactions.setNewPiStacking`,
:meth:`.Interactions.setNewPiCation`,
:meth:`.Interactions.setNewHydrophohic` method to replace it in the main
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




