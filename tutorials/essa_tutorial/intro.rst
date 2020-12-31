Introduction
===============================================================================

This tutorial demonstrates how to use ESSA for determining the essential sites
that would significantly alter the global/functional dynamics of a protein upon
ligand binding ([KB20]_). ESSA emulates ligand binding by increasing the node
density around each scanned residue. This is achieved by adding extra nodes at
the positions of the specific residue’s heavy atoms (other than the C-alpha atoms
that define the reference network in GNM or ANM). Changes in the global mode
spectrum in response to crowding near each scanned residue is measured by the
mean shift in the frequency of selected softest modes after pairwise matching of 
the modes between the reference and perturbed models. For quantifying *essentiality*,
we convert the mean shifts to z-scores. Essential residues correspond to hot
regions in terms of ligand binding and/or allostery.

Integration of ESSA with pocket information has been successful in predicting
the allosteric pockets for apo and holo structures of proteins ([KB20]_).
The pockets obtained from the Fpocket algorithm ([LGV09]_) are rank-ordered using 
the ESSA score calculated based on the residues lining each pocket.
Further screening using local hydrophobic density (LHD, a pocket feature of Fpocket)
improves the predictions. In the second part of this tutorial, we will demonstrate
how to use our automated protocol ([KB20]_) that combines ESSA scores,
pocket geometry and LHD for detecting allosteric pockets.

The example used in this tutorial for ESSA profile generation and the prediction
of allosteric pockets is TEM-1 beta-lactamase (PDB id: 1pzo), which we have
previously studied (see Figure S2 and Table S2 of [KB20]_ and the third image
on the ESSA webpage). 

Required Programs
-------------------------------------------------------------------------------

The latest version of ProDy_ is required for ESSA. Additionally, Fpocket 3.0 and
Pandas are required for the ESSA-based allosteric pocket prediction. 

Getting Started
-------------------------------------------------------------------------------

We recommend that you will follow this tutorial by typing commands in an
IPython session, e.g.::

  $ ipython

or with pylab environment::

  $ ipython --pylab

First, we will make necessary imports from ProDy_, NumPy_ and Matplotlib_
packages.

.. ipython:: python

    from prody import *
    from numpy import *
    from matplotlib.pyplot import *
    ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.

How to Cite
-------------------------------------------------------------------------------
If you benefited from ESSA in your research, please cite the following paper:

    .. [KB20] Kaynak B.T., Bahar I., Doruker P.,
        Essential site scanning analysis: A new approach for detecting sites that 
        modulate the dispersion of protein global motions,
        *Comput. Struct. Biotechnol. J.* **2020** 18:1577-1586.


Additionally, if you performed ESSA-based allosteric site prediction in your 
research, please also cite the following paper for Fpocket:

    .. [LGV09] Le Guilloux, V., Schmidtke P., Tuffery P.,
        Fpocket: An open source platform for ligand pocket detection,
        *BMC Bioinformatics* **2009** 10:168.
