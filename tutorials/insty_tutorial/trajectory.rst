.. _insty_tutorial:
=======

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


To compute hydrogen bonds for each frame of the simulation use
:func:`.calcHydrogenBondsTrajectory` function:

.. ipython:: python

   calcHydrogenBondsTrajectory(atoms, dcd)


Similarly, it can be done with other interaction types. Salt bridges
(residues with opposite changes) with :func:`.calcSaltBridgesTrajectory`:  

.. ipython:: python

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


In this particular example you will not have disulfide bonds but you can
compute it using :func:`.calcDisulfideBondsTrajectory`:

.. ipython:: python

   calcDisulfideBondsTrajectory(atoms, dcd)



Compute all availabe types of interactions
-------------------------------------------------------------------------------

First, we instantiate an :class:`.InteractionsTrajectory` instance which stores all the
information about interactions for protein structure for multiple frames.
With :meth:`.InteractionsTrajectory.calcProteinInteractionsTrajectory`, we can compute all
types of interactions such as hydrogen bonds, salt bridges, repulsive ionic bonding, 
Pi-cation, Pi-stacking, and hydrophobic) at once. Be aware that those
computations may take a while, depending on the size of the system and the number
of frames that are stored by the DCD file. Therefore, we recommend saving the
results as an *filename* file. *filename* file, here
*calcProteinInteractionsTrajectory.pkl*, can be reloaded and used with all availabe
functions and methods. 

.. ipython:: python

   interactionsTrajectory = InteractionsTrajectory('trajectory')
   interactionsTrajectory.calcProteinInteractionsTrajectory(atoms, dcd,
filename='calcProteinInteractionsTrajectory')


The results are displayed on the screen and they can be fetch by
using :meth:`.InteractionsTrajectory.getInteractions()` method.

.. ipython:: python

   interactionsTrajectory.getInteractions()


Moreover, we can display the evolution of each interaction type during the
simulation. There are the following types of plots: hydrogen bonds (*blue*),
salt bridges (*yellow*), hydrophobic interactions (*silver*), Pi-stacking
(*green*), Pi-cation (*orange*), repulsive ionic bonding (*red*), disulfide
bonds(*black*).  

.. ipython:: python

   number_of_counts = interactionsTrajectory.getTimeInteractions()


If the structure is stable we will not observe a lot of changes in protein
structure.


Similar to the single PDB analysis, we have an access to each interaction
type by using: :meth:`.InteractionsTrajectory.getHydrogenBonds` method, etc.

.. ipython:: python
   
   interactionsTrajectory.getHydrogenBonds()


.. ipython:: python
   
   interactionsTrajectory.getSaltBridges()


.. ipython:: python
   
   interactionsTrajectory.getHydrophobic()


.. ipython:: python
   
   interactionsTrajectory.getPiCation()


.. ipython:: python
   
   interactionsTrajectory.getPiStacking()


Once we compute interactions we can also select two which are interested for
us by using *selection* or *selection* and *selection2* if we want to
compare two chains of protein structure.
Here, we can display all interactions with residues with numbers between 100
and 106.

.. ipython:: python

    interactionsTrajectory.getInteractions(selection='resid 100 to 106')


We can apply the same selection to any type of interaction, for example:

.. ipython:: python
    
    interactionsTrajectory.getHydrophobic(selection='chain P and resid 90 to 100')


.. ipython:: python

    interactionsTrajectory.getHydrogenBonds(selection='chain P and resid 112 to 115')


.. ipython:: python

    interactionsTrajectory.getSaltBridges(selection='chain P and resid 100 to 120')


.. ipython:: python

    interactionsTrajectory.getRepulsiveIonicBonding(selection='chain P')



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


We can check whether the interactions were replaced. The repulsive ionic
bonding can be found by using :meth:`.getInteractions` and selecting *2*
(0 - hydrogen bonds, 1 - salt bridges, 2 - repulsive ionic bondnig, 
3 - Pi-Stacking, 4 - Pi-Cation, 5 - hydrophobic, 6- disulfide bonds).

.. ipython:: python

    interactionsTrajectory.getInteractions()[2]


.. ipython:: python

    interactionsTrajectory.getInteractions()[4]



Statistics
-------------------------------------------------------------------------------

Using :func:`.calcStatisticsInteractions` function, we can compute the statistics 
of interaction in the trajectory such as the average distance between residues
(usually the center of the mass, details are described in the function which
computes the specific type of interactions), standard deviation for the
distance value and weight. Weight is the number of counts for the whole
trajectory and devided by the number of frames in dcd file. Weight equal to 1
corresponds to the contact that was main the whole time of the simulation.
Note that weight can be >1 when the multiple contacts are present between the same
residues. 

For example:

.. ipython:: python

   interactions = interactionsTrajectory.getPiCation()
   calcStatisticsInteractions(interactions)


.. ipython:: python

   calcStatisticsInteractions(interactionsTrajectory.getHydrogenBonds())


For better visualization of those results we can use
func:`.showInteractionsGraph` which displayed results as a graph with a
residue-residue pairs of interactions. The intensity of the color of the
lines connecting two residues corresponds to the number of counts. Darker
lines are assigned to the most frequent appearence of interaction. The
distance between pairs corresponds to the average distance accross
all the frames. Moreover, ovals with residue names are color-coded: acidic
residues: *red*, basic: *blue*, polar: *green*, non-polar: *silver*, and
proline: *pink*. The same function can be applied for ensemble PDB (more
example can be found there).

.. ipython:: python

    showInteractionsGraph(statistics, code='1-letter', cutoff=0.5)



Compare two independent frames
-------------------------------------------------------------------------------

With this analysis we can also compare interactions between frames. Below we
will compute hydrogen bonds for frame 0 and frame 18 and we will compare it
using :func:`.compareInteractions` function. That function will be helpful
is checking the difference between interactions. The results will be saved
as *diff_fr0_vsfr18.dat* file.

.. ipython:: python

    frame0 = dcd.getFrame(0)
    at0 = frame0.getAtoms()
    hb0 = calcHydrogenBonds(at0.select('protein')) 


.. ipython:: python

    frame18 = dcd.getFrame(18)
    at18 = frame18.getAtoms()
    hb18 = calcHydrogenBonds(at18.select('protein'))


.. ipython:: python

    compareInteractions(hb0, hb18, filename='diff_fr0_vsfr18.dat')


We can also all tools which are shown for single PDB analysis in this
tutorial. For example, we can compute all interactions for frame0 and
frame18 and display the inetractions:

.. ipython:: python

    interactions0 = Interactions()
    interactions0.calcProteinInteractions(at0)
    matrix0 = interactions0.buildInteractionMatrix()


.. ipython:: python

    interactions0.showCumulativeInteractionTypes()


.. ipython:: python

    interactions18 = Interactions()
    interactions18.calcProteinInteractions(at18)
    matrix18 = interactions18.buildInteractionMatrix()


.. ipython:: python

    interactions18.showCumulativeInteractionTypes()



Parse previously saved data
-------------------------------------------------------------------------------

To upload and further use the interactions data use
:meth:`.InteractionsTrajectory.parseInteractions` function:

.. ipython:: python

   interactionsTrajectory2 = InteractionsTrajectory('5kqm_import')
   interactionsTrajectory2.parseInteractions('calcProteinInteractionsTrajectory.pkl')


After uploading, we have access to all data, for example:

.. ipython:: python

   interactionsTrajectory2.getHydrophohic()


.. ipython:: python

   calcStatisticsInteractions(interactionsTrajectory2.getHydrogenBonds())


.. ipython:: python

   interactionsTrajectory2.getTimeInteractions()


