.. _signdy-data:

Data Collection
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures and build a :class:`.PDBEnsemble`. This can be achieved by 
multiple routes: a query search of the PDB using :func:`.blastPDB` or :func:`.searchDali`, 
extraction of PDB IDs from the Pfam or CATH database, or input of a pre-defined list. 

We demonstrate the usage of SignDy with a pre-defined list of transporter proteins sharing 
the common LeuT fold [YS13]_. These proteins cycle through four typical states to transport 
a substrate molecule: outward-facing open (OFo), outward-facing closed (OFc), inward-facing 
open (IFo), inward-facing closed (IFc), but only the first three states have PDB structures 
available. If you know how to prepare an ensemble of structural homologs and wish to skip 
this part, you can download the ensemble file used in [SZ18]_ from here and proceed to the 
next tutorial.

First, make necessary imports from ProDy_ and Matplotlib_ packages if you haven't already.

.. ipython:: python

    from prody import *
    from pylab import *
    ion()

Prepare ensemble
-------------------------------------------------------------------------------

For convinience and clarity, we define LeuT folds in separate lists taxonomically. For example,
the PDB identifiers for bacterial Leucine transporters are defined as follows:

.. ipython:: python

    LeuTs = ['2A65', '2Q6H', '2Q72', '2QB4', '2QEI', '2QJU', '3F3A', '3F3C', 
             '3F3D', '3F3E', '3F48', '3F4I', '3F4J', '3GJC', '3GJD', '3GWU', 
             '3GWV', '3GWW', '3MPN', '3MPQ', '3QS4', '3QS5', '3QS6', '3TT1', 
             '3TU0', '3USG', '3USI', '3USJ', '3USK', '3USL', '3USM', '3USO', 
             '3USP', '4FXZ', '4FY0', '4HMK', '4HOD', '4MM4', '4MM5', '4MM6', 
             '4MM7', '4MM8', '4MM9', '4MMA', '4MMB', '4MMC', '4MMD', '4MME', 
             '4MMF', '3TT3']

Despite the fact that bacterial Leucine transporters can form dimers, we will only take the 
chain A in each structure:

.. ipython:: python

    LeuTs = [protID + 'A' for protID in LeuTs]

Note that in the above line, we use `list comprehension`_ to add a letter 'A' to each PDB 
identifier in the list. We define other LeuT folds similarily:

.. ipython:: python

    DATs = ['4M48', '4XNU', '4XNX', '4XP1', '4XP4', '4XP5', '4XP6', 
            '4XP9', '4XPA', '4XPB', '4XPF', '4XPG', '4XPH', '4XPT']
    DATs = [protID + 'A' for protID in DATs]
    MhsTs = ['4US4A', '4US3A']
    vSGLTs = ['2XQ2A']
    Mhp1s = ['2JLN', '2X79', '4D1A', '4D1B', '4D1C', '4D1D']
    Mhp1s = [protID + 'A' for protID in Mhp1s]
    BetPs = ['2WITA', '2WITB', '2WITC', '3P03A', '3P03B', '3P03C', 
             '4AINA', '4AINB', '4AINC', '4C7RA', '4C7RB', '4C7RC', 
             '4DOJA', '4DOJB', '4DOJC', '4LLHA', '4LLHB', '4LLHC']
    AdiCs = ['3L1L', '3LRB', '3LRC', '3NCY', '3OB6', '5J4I', '5J4N']
    AdiCs = [protID + 'A' for protID in AdiCs]
    CaiTs = ['4M8JA', '2WSXA', '2WSXB', '2WSXC', '2WSWA', '3HFXA']

:func:`.parsePDB` allows us to parse multiple structures all at once, and we can use it to 
load all the PDB structures into ProDy_ in one line. We only need the alpha carbon for our 
purpose, so we set ``subset='ca'``:

.. ipython:: python

    pdb_ids = LeuTs + DATs + MhsTs + vSGLTs + Mhp1s + BetPs + AdiCs + CaiTs
    pdbs = parsePDB(*pdb_ids)
    len(pdbs)

Any element in the list *pdbs* should be an :class:`.AtomGroup` instance. We can conveniently 
feed this list to :func:`.buildPDBEnsemble` and let it build an :class:`.PDBEnsemble` for downstream 
analyses. We use set ``mapping=ce`` to tell the function to use a structure alignment algorithm, 
CEalign [IS98]_, for building the ensemble. We also set ``seqid=0`` and ``overlap=0`` to make sure 
we apply no threshold of sequence identity or coverage/overlap to the building process. 

.. ipython:: python

    ens = buildPDBEnsemble(pdbs, mapping='ce', seqid=0, overlap=0, title='LeuT', subset='ca')
    ens

Finally we save the ensemble for later processing:

.. ipython:: python

    saveEnsemble(ens, 'LeuT')

A refiner alignment procedure was adopted in the [SZ18]_ paper. A representative structure is chosen 
from each subtype of the proteins, e.g. LeuT, DAT, etc., and they are aligned to the LeuT representative 
using CEalign [IS98]_. Then the rest are aligned to the representative structure of their own kind using 
the pairwise alignment algorithm because they are sequentially the same despite small differences. The 
ensemble used in the [SZ18]_ paper is provided in the download files and will be used in the next tutorial, 
but you are also welcome to use the ensemble we created using above code.

.. _`Structure Analysis Tutorial`: http://prody.csb.pitt.edu/tutorials/structure_analysis/blastpdb.html
.. _`list comprehension`: https://docs.python.org/2/tutorial/datastructures.html#list-comprehensions

.. [YS13] Shi Y.
    Common folds and transport mechanisms of secondary active transporters.
    *Annu. Rev. Biophys.* **2013** 42:51-72
