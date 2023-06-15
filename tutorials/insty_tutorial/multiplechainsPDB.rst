.. _insty_tutorial:

PDB structure with multiple chains
===============================================================================

This time we will use protein with two chains, lipoxygenase (**7LAF**) which
contain chain A and chain B. First, we will add missing hydrogens to the
protein structures and then we will perform analysis of interactions between
two chains. 

Add missing hydrogen atoms to the structure
-------------------------------------------------------------------------------

We start by fetching the PDB file and adding missing hydrogens using
Openbabel.

.. ipython:: python

   fetchPDB('7laf', compressed=False)

::

   addMissingAtoms('7laf.pdb', method='openbabel')


.. ipython:: python

   atoms = parsePDB('addH_7laf.pdb').select('protein')


.. ipython:: python

   interactions = Interactions('7laf')


To compute all interactions:

.. ipython:: python

   all_interactions = interactions.calcProteinInteractions(atoms)


.. ipython:: python

   interactions.getHydrogenBonds(selection='chain A', selection2='chain B')


.. ipython:: python

   interactions.getSaltBridges(selection='chain A', selection2='chain B')


.. ipython:: python

   interactions.getHydrophobic(selection='chain A', selection2='chain B')


.. ipython:: python

   interactions.getPiStacking(selection='chain A', selection2='chain B')


.. ipython:: python
   
   interactions.getPiCation(selection='chain A', selection2='chain B')


.. ipython:: python

   interactions.getRepulsiveIonicBonding(selection='chain A', selection2='chain B')


Non-zero interactions could be futher saved and used in VMD_ program to
display them:

.. ipython:: python

   showProteinInteractions_VMD(atoms, interactions.getHydrogenBonds(), color='blue', filename='HBs_7laf.tcl')
   showProteinInteractions_VMD(atoms, interactions.getSaltBridges(), color='yellow',filename='SBs_7laf.tcl')
   showProteinInteractions_VMD(atoms, interactions.getHydrophobic(), color='silver',filename='HPh_7laf.tcl')











