.. _essa_prediction:

Prediction of allosteric pockets
===============================================================================

First, we will make necessary imports from ProDy_, NumPy_ and Matplotlib_
packages if you haven't already done it:

.. ipython:: python

    from prody import *
    from numpy import *
    from matplotlib.pyplot import *
    ion()


Allosteric pocket prediction requires Fpocket 3.0 and Pandas installed in your 
system. The first step is the pocket hunting, which is automatically carried out 
in the background, by calling :func:`.scanPockets`. This method parses the pocket 
features provided by Fpocket, and all identified pockets are stored in a folder 
ending with *_out* in the current working directory. Additionally, maximum/median 
ESSA z-scores are assigned to each pocket based on the residues forming it.

.. ipython:: python

    essa.scanPockets()

Pocket features that are stored in a Pandas DataFrame can be obtained by 
:func:`.getPocketFeatures`, and saved as a Python pickle file by 
:func:`.savePocketZscores`.

Key features of the pockets to be used in the prediction, namely ESSA and local 
hydrophobic density (LHD) z-scores, can be listed by :func:`.getPocketZscores`.

.. ipython:: python

    essa.getPocketZscores()

The prediction protocol ranks the pockets with respect to their ESSA and LHD 
z-scores. At the same time, the pockets with negative LHD z-scores are also 
filtered out as allosteric sites are known to have relatively higher LHD. For 
the details of this protocol, please refer to the original ESSA article ([KB20]_).

Ranking of the pockets can be performed and obtained by :func:`.rankPockets` and 
:func:`.getPocketRanks`, respectively.

.. ipython:: python

    essa.rankPockets()
    essa.getPocketRanks()

Pocket 6 with the highest ESSA_max score has been identified as the only allosteric 
pocket in this structure. Interestingly, other pockets have been filtered 
out due to their negative LHD z-scores. This is a large pocket that includes CBT 
allosteric ligands at A300 and A301, as well as a part of the orthosteric ligand 
(shown as image three on the ESSA webpage).

In order to visualize the pockets, the `.pqr` file, an output of Fpocket needs 
to be opened by PyMOL or VMD together with the original pdb file. 

Pocket z-scores and ranks can be saved by :func:`.savePocketZscores` and 
:func:`.writePocketRankstoCSV`, respectively.

.. ipython:: python

    essa.savePocketZscores()
    essa.writePocketRankstoCSV()