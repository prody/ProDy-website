.. _interactions_analysis:

Analysis of interactions for a trajectory (DCD file)
===============================================================================

This example shows how to compute interactions for a trajectory performed
using NAMD_ software for a small protein tyrosine phosphatase LMW-PTP 
in a complex with inhibitor MES (**5KQM**) and visualize the results using 
Matplotlib_ library and VMD_ program. 

In the tutorial, we will use already preapared files for
simulation (*PDB* and *DCD* file).


Parse structure with trajectory
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on

We start by parsing PDB and DCD files whcih contain LMW-PTP protein
structure (avalable as the tutorial files). PDB file contains the
coordinates of protein structure with water and counter ions. DCD
file is a binary file which contain short simulation commputed in NAMD_
package (20 frames). Commands shown below are explained in *Trajectory
Analysis* tutorial.

.. ipython:: python

   PDBfile = '5kqm_all_sci.pdb'
   DCDfile = 'NAMD_D2_co100.dcd'
   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)


Compute interactions
-------------------------------------------------------------------------------

To compute hydrogen bonds for each frame of the simulation use
:func:`.calcHydrogenBondsDCD` function:

.. ipython:: python

   calcHydrogenBondsDCD(atoms, dcd)


Similarly, it can be done with other interaction types. Salt bdriges
(residues with oposite changes) with :func:`.calcSaltBridgesDCD`:  

.. ipython:: python

   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)
   
   calcSaltBridgesDCD(atoms, dcd)


Repulsive Ionic Bonding using :func:`.calcRepulsiveIonicBondingDCD` for residues with
the same charges:

.. ipython:: python

   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)

   calcRepulsiveIonicBondingDCD(atoms, dcd, distA=7)


Pi-Stacking interactions using :func:`.calcPiStackingDCD`:

.. ipython:: python

   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)

   calcPiStackingDCD(atoms, dcd, distA=5)


Pi-Cation interactions using :func:`.calcPiCationDCD`:

.. ipython:: python

   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)

   calcPiCationDCD(atoms, dcd)

Hydrophobic interactions using :func:`.calcHydrophohicDCD`:

.. ipython:: python

   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)

   calcHydrophohicDCD(atoms, dcd)



Compute all availabe type of interactions at once
-------------------------------------------------------------------------------

First, we need to parse PDB and DCD file:

.. ipython:: python
  
   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)


Next, we instantiate an :class:`.InteractionsDCD` instance which stores all the
information about interactions for protein structure for multiple frames.
With :meth:`.InteractionsDCD.calcProteinInteractionsDCD` we can compute all
type of interactions such as: hydrogen bonds, salt bridges, repulsive ionic bonding, 
Pi-cation, Pi-stacking, and hydrphobic) at once. Be aware that those
computations may take a while depending on the size of the system and number
of frames which are storaged by DCD file, therefore, we recomend to save the
results as an *output* file. *Output* file, here
*calcProteinInteractionsDCD.pkl*, can be reloaded and used with all availabe
functions and methods. 

.. ipython:: python

   interactionsDCD = InteractionsDCD('trajectory')
   interactionsDCD.calcProteinInteractionsDCD(atoms, dcd, output='interactions_data_5kqm')


The results are displayed on the screen but they can display them also
using :meth:`.InteractionsDCD.getInteractions()` method.

.. ipython:: python

   interactionsDCD.getInteractions()


Moreover, we can display the evolution of each interaction type during the
simulation. There are following type of plots: hydrogen bonds (*blue*), salt
bridges (*yellow*), hydrophobic interactions (*silver*), Pi-stacking
(*green*), Pi-cation (*orange*), repulsive ionic bonding (*red*).  

.. ipython:: python

   number_of_counts = interactionsDCD.getTimeInteractions()


Similar to the single PDB analysis, we have an access to each interaction
type by using: :meth:`.InteractionsDCD.getHydrogenBonds` method, etc.

.. ipython:: python
   
   interactionsDCD.getHydrogenBonds()


Change selection criteria for interaction type
-------------------------------------------------------------------------------

The :meth:`.interactionsDCD.calcProteinInteractionsDCD` method compute interactions 
using default paramaters for interactions. However, it can be changed
accoridng to our needs. To do that, we need to recalculate selected type of
interactions. 

We can do it using following functions: :func:`.calcHydrogenBondsDCD`,
:func:`.calcHydrogenBondsDCD`, :func:`.calcSaltBridgesDCD`,
:func:`.calcRepulsiveIonicBondingDCD`, :func:`.calcPiStackingDCD`,
:func:`.calcPiCationDCD`, :func:`.calcHydrophohicDCD`, and use
:meth:`.InteractionsDCD.setNewHydrogenBondsDCD`,
:meth:`.InteractionsDCD.setNewSaltBridgesDCD`,
:meth:`.InteractionsDCD.setNewRepulsiveIonicBondingDCD`,
:meth:`.InteractionsDCD.setNewPiStackingDCD`,
:meth:`.InteractionsDCD.setNewPiCationDCD`,
:meth:`.InteractionsDCD.setNewHydrophohicDCD` method to replace it in the main
Instance. 

For example:

.. ipython:: python

   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)
   
   newRIB = calcRepulsiveIonicBondingDCD(atoms, dcd, distA=8)
   interactionsDCD.setNewRepulsiveIonicBondingDCD(newRIB)
   
.. ipython:: python

   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)
   
   newPiCation = calcPiCationDCD(atoms, dcd, distA=6)
   interactionsDCD.setNewPiCationDCD(newPiCation)


Statistics
-------------------------------------------------------------------------------

Using :func:`.calcStatisticsInteractions` function, we can compute the statistics 
of interaction in the tajectory such as number of counts, avarage distance between 
residues (usualy center of the mass, details are described in function which compute 
specific type of interactions) and standard deviation. For example:


.. ipython:: python

   interactions = interactionsDCD.getPiStacking()
   calcStatisticsInteractions(interactions)


.. ipython:: python

   calcStatisticsInteractions(interactionsDCD.getHydrogenBonds())


Parse previously saved interactions
-------------------------------------------------------------------------------

To upload and further use the inetractions data use
:meth:`.InteractionsDCD.parseInteractions` function:

.. ipython:: python

   interactionsDCD2 = InteractionsDCD('5kqm_import')
   interactionsDCD2.parseInteractions('interactions_data_5kqm.pkl')


After uploading, we have access to all data, for example:

.. ipython:: python

   interactionsDCD2.getHydrophohic()

.. ipython:: python

   calcStatisticsInteractions(interactionsDCD2.getHydrogenBonds())

.. ipython:: python

   interactionsDCD2.getTimeInteractions()


