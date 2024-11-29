Introduction
===============================================================================

ClustENMD [KD21]_ is a highly efficient, unbiased conformational search algorithm, 
which is suitable for generating atomistic conformers for ensemble docking purposes and for large, 
supramolecular assemblies like the ribosome [KD16]_. In this hybrid method, the following steps are 
carried out at each generation/cycle: (1) conformer generation by random sampling along global ANM modes,
(2) hierarchical clustering of generated conformers and (3) relaxation of cluster representatives by short 
MD simulations using OpenMM [E17]_. The relaxed conformers are then fed back to Step 1, each being used as 
a starting point for a new generation of conformers. This iterative procedure (Steps 1-3) is repeated for 
several generations to allow for sufficiently large excursions from the initial structure. Thus, conformational 
sampling can be efficiently performed for highly flexible systems composed of proteins, RNA and/or DNA chains. 
Furthermore, the generated ensemble can be analyzed and utilized using the available tools in ProDy.

ANMD [CM22]_ is a simple conformers generation algorithm that explores motions along single ANM modes.
The normal mode calculation is only performed once and a set of structures along each of the first few modes
is subjected to energy minimization, producing physically reasonable conformers near the initial energy well.
A key feature of ANMD is that each successive mode is sampled with a lower amplitude related to its frequency.

The first part of this tutorial demonstrates how to use ClustENMD to perform conformational sampling for the homo-dimeric enzyme 
HIV-1 protease in an open conformation without any ligand (PDB id: 1tw7). Furthermore, we will show the application of
ProDy ensemble analysis tools to study the conformers and generate their population distribution.

The last part demonstrates how to use ANMD to explore the first two modes of a metabotropic glutamate receptor
N-terminal domain (PDB id: 1ewk). 

Required Programs
-------------------------------------------------------------------------------

The latest versions of ProDy_, OpenMM_, and PDBFixer_ are required for ClustENMD.

.. _OpenMM: https://openmm.org/
.. _PDBFixer: https://github.com/openmm/pdbfixer

Getting Started
-------------------------------------------------------------------------------

We recommend that you will follow this tutorial by typing commands in an
IPython session, e.g.::

  $ ipython

or with pylab environment::

  $ ipython --pylab

First, we will make necessary imports from ProDy_, NumPy_, and Matplotlib_ packages.

.. ipython:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from prody import *
    plt.ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.

How to Cite
-------------------------------------------------------------------------------
If you benefited from ClustENMD in your research, please cite the following paper:

.. [KD16] Kurkcuoglu Z., Bahar I., and Doruker P., 
   ClustENM: ENM-Based Sampling of Essential Conformational Space at Full Atomic
   Resolution, *J Chem Theory Comput* **2016** 12: 4549.

.. [KD21] Kaynak B.T., Zhang S., Bahar I., and Doruker P., 
   ClustENMD: Efficient sampling of biomolecular conformational space at atomic resolution,
   *Bioinformatics* **2021** 37(21): 3956–3958. 

.. [CM22] Mary Hongying Cheng, James M Krieger, Anupam Banerjee, Yufei Xiang, 
   Burak Kaynak, Yi Shi, Moshe Arditi, Ivet Bahar. 
   Impact of new variants on SARS-CoV-2 infectivity and neutralization: 
   A molecular assessment of the alterations in the spike-host protein 
   interactions, *iScience* **2022** 25(3):103939.


Additionally, please also cite the following paper for OpenMM:

.. [E17] Eastman P., et al. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics, *PLoS Comput Biol* **2017** 13:e1005659.
