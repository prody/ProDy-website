.. _insty_tutorial:

Trajectory analysis
===============================================================================

This example shows how to compute interactions for a trajectory performed
using NAMD_ software for a small protein tyrosine phosphatase LMW-PTP 
in a complex with inhibitor MES (**5KQM**) and visualize the results using 
Matplotlib_ library and VMD_ program. 

In the tutorial, we will use already prepared files for
simulation (*PDB* and *DCD* file).


Parse structure with trajectory
-------------------------------------------------------------------------------

We start by parsing PDB and DCD files which contain LMW-PTP protein
structure (available as the tutorial files). PDB file contains the
coordinates of protein structure with water and counter ions. DCD
file is a binary file that contains a short simulation computed in NAMD_
package (20 frames). The commands shown below are explained in *Trajectory
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
:func:`.calcHydrogenBondsTrajectory` function:

.. ipython:: python

   calcHydrogenBondsTrajectory(atoms, dcd)


Similarly, it can be done with other interaction types. Salt bridges
(residues with opposite changes) with :func:`.calcSaltBridgesTrajectory`:  

.. ipython:: python

   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)
   
   calcSaltBridgesTrajectory(atoms, dcd)


Repulsive Ionic Bonding using :func:`.calcRepulsiveIonicBondingTrajectory` for residues with
the same charges:

.. ipython:: python

   calcRepulsiveIonicBondingTrajectory(atoms, dcd, distA=7)


Pi-Stacking interactions using :func:`.calcPiStackingTrajectory`:

.. ipython:: python

   calcPiStackingTrajectory(atoms, dcd, distA=5)


Pi-Cation interactions using :func:`.calcPiCationTrajectory`:

.. ipython:: python

   calcPiCationTrajectory(atoms, dcd)

Hydrophobic interactions using :func:`.calcHydrophohicTrajectory`:

.. ipython:: python

   calcHydrophohicTrajectory(atoms, dcd)



Compute all availabe types of interactions
-------------------------------------------------------------------------------

First, we need to parse PDB and DCD file:

.. ipython:: python
  
   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)


Next, we instantiate an :class:`.InteractionsTrajectory` instance which stores all the
information about interactions for protein structure for multiple frames.
With :meth:`.InteractionsTrajectory.calcProteinInteractionsTrajectory`, we can compute all
types of interactions such as hydrogen bonds, salt bridges, repulsive ionic bonding, 
Pi-cation, Pi-stacking, and hydrophobic) at once. Be aware that those
computations may take a while, depending on the size of the system and the number
of frames that are stored by the DCD file. Therefore, we recommend saving the
results as an *output* file. *Output* file, here
*calcProteinInteractionsTrajectory.pkl*, can be reloaded and used with all availabe
functions and methods. 

.. ipython:: python

   interactionsTrajectory = InteractionsTrajectory('trajectory')
   interactionsTrajectory.calcProteinInteractionsTrajectory(atoms, dcd, output='interactions_data_5kqm')


The results are displayed on the screen but they can display them also
using :meth:`.InteractionsTrajectory.getInteractions()` method.

.. ipython:: python

   interactionsTrajectory.getInteractions()


Moreover, we can display the evolution of each interaction type during the
simulation. There are the following types of plots: hydrogen bonds (*blue*),
salt bridges (*yellow*), hydrophobic interactions (*silver*), Pi-stacking
(*green*), Pi-cation (*orange*), repulsive ionic bonding (*red*).  

.. ipython:: python

   number_of_counts = interactionsTrajectory.getTimeInteractions()


Similar to the single PDB analysis, we have an access to each interaction
type by using: :meth:`.InteractionsTrajectory.getHydrogenBonds` method, etc.

.. ipython:: python
   
   interactionsTrajectory.getHydrogenBonds()


Change selection criteria for interaction type
-------------------------------------------------------------------------------

The :meth:`.interactionsTrajectory.calcProteinInteractionsTrajectory` method computes
interactions using default parameters for interactions. However, it can be
changed according to our needs. To do that, we need to recalculate the
selected types of interactions. 

We can do it using the following functions: :func:`.calcHydrogenBondsTrajectory`,
:func:`.calcHydrogenBondsTrajectory`, :func:`.calcSaltBridgesTrajectory`,
:func:`.calcRepulsiveIonicBondingTrajectory`, :func:`.calcPiStackingTrajectory`,
:func:`.calcPiCationTrajectory`, :func:`.calcHydrophohicTrajectory`, and use
:meth:`.InteractionsTrajectory.setNewHydrogenBondsTrajectory`,
:meth:`.InteractionsTrajectory.setNewSaltBridgesTrajectory`,
:meth:`.InteractionsTrajectory.setNewRepulsiveIonicBondingTrajectory`,
:meth:`.InteractionsTrajectory.setNewPiStackingTrajectory`,
:meth:`.InteractionsTrajectory.setNewPiCationTrajectory`,
:meth:`.InteractionsTrajectory.setNewHydrophohicTrajectory` method to replace it in the main
Instance. 

For example:

.. ipython:: python

   newRIB = calcRepulsiveIonicBondingTrajectory(atoms, dcd, distA=8)
   interactionsTrajectory.setNewRepulsiveIonicBondingTrajectory(newRIB)
   
.. ipython:: python

   newPiCation = calcPiCationTrajectory(atoms, dcd, distA=6)
   interactionsTrajectory.setNewPiCationTrajectory(newPiCation)


Statistics
-------------------------------------------------------------------------------

Using :func:`.calcStatisticsInteractions` function, we can compute the statistics 
of interaction in the trajectory such as the number of counts, average distance
between residues (usually the center of the mass, details are described in
the function which computes the specific type of interactions), and
standard deviation. For example:


.. ipython:: python

   interactions = interactionsTrajectory.getPiStacking()
   calcStatisticsInteractions(interactions)


.. ipython:: python

   calcStatisticsInteractions(interactionsTrajectory.getHydrogenBonds())


Parse previously saved data
-------------------------------------------------------------------------------

To upload and further use the interactions data use
:meth:`.InteractionsTrajectory.parseInteractions` function:

.. ipython:: python

   interactionsTrajectory2 = InteractionsTrajectory('5kqm_import')
   interactionsTrajectory2.parseInteractions('interactions_data_5kqm.pkl')


After uploading, we have access to all data, for example:

.. ipython:: python

   interactionsTrajectory2.getHydrophohic()

.. ipython:: python

   calcStatisticsInteractions(interactionsTrajectory2.getHydrogenBonds())

.. ipython:: python

   interactionsTrajectory2.getTimeInteractions()


