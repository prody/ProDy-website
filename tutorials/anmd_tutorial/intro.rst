Introduction
===============================================================================

ANMD [CM22]_ is a simple conformers generation algorithm that explores motions along single ANM modes.
The normal mode calculation is only performed once and a set of structures along each of the first few modes
is subjected to energy minimization, producing physically reasonable conformers near the initial energy well.
A key feature of ANMD is that each successive mode is sampled with a lower amplitude related to its frequency.

This tutorial demonstrates how to use ANMD to explore the first two modes of a metabotropic glutamate receptor
N-terminal domain (PDB id: 1ewk). 

Required Programs
-------------------------------------------------------------------------------

The latest versions of ProDy_, OpenMM_, and PDBFixer_ are required for ANMD.

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
If you benefited from ANMD in your research, please cite the following paper:

.. [CM22] Mary Hongying Cheng, James M Krieger, Anupam Banerjee, Yufei Xiang, 
   Burak Kaynak, Yi Shi, Moshe Arditi, Ivet Bahar. 
   Impact of new variants on SARS-CoV-2 infectivity and neutralization: 
   A molecular assessment of the alterations in the spike-host protein 
   interactions, *iScience* **2022** 25(3):103939.

Additionally, please also cite the following paper for OpenMM:

.. [E17] Eastman P., et al. OpenMM 7: Rapid development of high performance algorithms for molecular dynamics, *PLoS Comput Biol* **2017** 13:e1005659.
