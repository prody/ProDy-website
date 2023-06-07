.. _insty_tutorial:

Ensemble PDB analysis
===============================================================================

This example shows how to compute interactions for an Ensemble PDB
(e.g. NMR data). The example is prepared for a NMR structure of ubiquitin 
(**2K39**) and visualize the results using Matplotlib_ library and VMD_ program. 


Parse structure
-------------------------------------------------------------------------------

We start by parsing PDB file which contain multiple conformations of
ubiquitin structure.

.. ipython:: python

   atoms = parsePDB('2k39')


Compute interactions for an Ensemble PDB
-------------------------------------------------------------------------------

To compute hydrogen bonds for each frame use :func:`.calcHydrogenBondsTrajectory`
function:

.. ipython:: python

   calcHydrogenBondsTrajectory(atoms)


Similarly, it can be done with other interaction types. Salt bridges with
:func:`.calcSaltBridgesTrajectory`:  

.. ipython:: python

   calcSaltBridgesTrajectory(atoms)


Repulsive Ionic Bonding using :func:`.calcRepulsiveIonicBondingTrajectory`
for residues with the same charges:

.. ipython:: python

   calcRepulsiveIonicBondingTrajectory(atoms)


Pi-Stacking interactions using :func:`.calcPiStackingTrajectory`:

.. ipython:: python

   calcPiStackingTrajectory(atoms)


Pi-Cation interactions using :func:`.calcPiCationTrajectory`:

.. ipython:: python

   calcPiCationTrajectory(atoms)


Hydrophobic interactions using :func:`.calcHydrophohicTrajectory`:

.. ipython:: python

   calcHydrophobicTrajectory(atoms)


.. ipython:: python

   calcDisulfideBondsTrajectory(atoms)



Select particular frames and change default parameters of the interactions
-------------------------------------------------------------------------------

The default parameters which are assigned to the interaction types could be
changed as follows:

.. ipython:: python
  
   calcHydrogenBondsTrajectory(atoms, distA=2.7, angle=35, cutoff_dist=10)

Similarly for other interactions type. Moreover, we can also select frames
which we would like to analyze as well as the selection with the protein
structure. Below you will find such examples:

.. ipython:: python
  
   calcPiCationTrajectory(atoms, distA=7, start_frame=15, stop_frame=20)


.. ipython:: python
  
   calcHydrophobicTrajectory(atoms, start_frame=10, stop_frame=13, selection='resid 50 to 60')



Compute all availabe types of interactions at once
-------------------------------------------------------------------------------

Next, we instantiate an :class:`.InteractionsTrajectory` instance which stores all the
information about interactions for protein structure for multiple frames.
With :meth:`.InteractionsTrajectory.calcProteinInteractionsTrajectory`, we can compute
all types of interactions such as hydrogen bonds, salt bridges, repulsive ionic bonding, 
Pi-cation, Pi-stacking, and hydrophobic) at once. Be aware that those
computations may take a while, depending on the size of the system and the number
of frames that are stored by the Ensemble PDB file. Therefore, we recommend saving the
results as an *output* file. *Output* file, *calcProteinInteractionsEnseblePDB.pkl*,
can be reloaded and used with all availabe functions and methods. 

.. ipython:: python

   interactionsTrajectoryNMR = InteractionsTrajectory('ensambleNMR')
   interactionsTrajectoryNMR.calcProteinInteractionsTrajectory(atoms,
   filename='calcProteinInteractionsEnseblePDB.pkl')


The results can be displayed using :meth:`.getTimeInteractions` where all
the interactions are displayed and could be tracked per each
confrmation (frame in the Ensemble PDB file).

.. ipython:: python

   number_of_counts = interactionsTrajectoryNMR.getTimeInteractions()
   

Each interactions type could be further counted with some additional
quantitative analysis using :func:`.calcStatisticsInteractions`:

.. ipython:: python

   statistics = calcStatisticsInteractions(interactionsTrajectoryNMR.getHydrogenBonds())


To provide a better way for visualization of those results another function
func:`.showInteractionsGraph` could be used which provides graph with a
residue-residue pairs of interactions. The intensity of the color of the
lines connecting two residues corresponds to the number of counts. Draker
lines are assigned to the most frequent appearence of interaction. The
distance between pairs corresponds to the average distance accross
all the frames. Moreover, ovals with residue names are color-coded: acidic
residues: red, basic: blue, polar: green, non-polar: silver, and proline: pink.

Below an example with additional parameters: *1-letter* code of residues
which be used instead of 3-letter code, *cutoff* = 3 for the number of counts
for residue interaction, *font_size* for the residue names displayed on the
graph and *seed* which is a random number which can help to organize the
graph in a nicer way.


.. ipython:: python

   showInteractionsGraph(statistics, code='1-letter', cutoff=3, font_size=8, seed=42)


.. ipython:: python

   showInteractionsGraph(statistics, code='1-letter', cutoff=50,
   font_size=16, node_distance=3, seed=1)


Selection of the residue pairs can be made as needed by choosing pairs which
higher number of counts or by changing the selection to certain region:

.. ipython:: python

   showInteractionsGraph(statistics, code='1-letter', cutoff=50, font_size=16,
   node_distance=3, seed=1)


.. ipython:: python

   hbs_20to30 = interactionsTrajectoryNMR.getHydrogenBonds(selection='resid 20 to 30')
   statistics2 = calcStatisticsInteractions(hbs_20to30)
   showInteractionsGraph(statistics2)


The selection can be made at different stages of analysis. The example below
shows how to analyse only certain frames (from 5th - 10th frame) for
residues number betwee 10 and 30.

.. ipython:: python

   interactionsTrajectoryNMR.calcProteinInteractionsTrajectory(atoms, start_frame=5, stop_frame=10, selection='resid 10 to 30')



