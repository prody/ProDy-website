.. _signdy-data:

Data Collection
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures and build a :class:`.PDBEnsemble`. This can be achieved by 
multiple routes: a query search of the PDB using :func:`.blastPDB` or :func:`.searchDali`, 
extraction of PDB IDs from the Pfam or CATH database, or input of a pre-defined list. 

We demonstrate the usage of SignDy with a pre-defined list of transporter proteins sharing 
the common LeuT fold [Y13]_. These proteins cycle through four typical states to transport 
a substrate molecule: outward-facing open (OFo), outward-facing closed (OFc), inward-facing 
open (IFo), inward-facing closed (IFc), but only the first three states have PDB structures 
available. 

First, make necessary imports from ProDy_ and Matplotlib_ packages if you haven't already.

.. ipython:: python

    from prody import *
    from pylab import *
    ion()

Define the LeuT fold PDB identifiers and building a PDBEnsemble
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

    LeuTs = [leut + 'A' for leut in LeuTs]

In the above line, we use list comprehension to add a letter 'A' to each PDB identifier in the 
list. We define other LeuT folds similarily:

.. ipython:: python

    MhsTs = ['4US4', '4US3']
    DATs = ['4M48', '4XNU', '4XNX', '4XP1', '4XP4', '4XP5', '4XP6', 
            '4XP9', '4XPA', '4XPB', '4XPF', '4XPG', '4XPH', '4XPT']
    vSGLTs = ['2XQ2']
    Mhp1s = ['2JLN', '2X79', '4D1A', '4D1B', '4D1C', '4D1D']
    BetPs = ['2WITA', '2WITB', '2WITC', '3P03A', '3P03B', '3P03C', 
             '4AINA', '4AINB', '4AINC', '4C7RA', '4C7RB', '4C7RC', 
             '4DOJA', '4DOJB', '4DOJC', '4LLHA', '4LLHB', '4LLHC']
    AdiCs = ['3L1L', '3LRB', '3LRC', '3NCY', '3OB6', '5J4I', '5J4N']
    CaiTs = ['4M8JA', '2WSXA', '2WSXB', '2WSXC', '2WSWA', '3HFXA']

.. ipython:: python

    pdbs = parsePDB(*pdb_ids, subset='ca')
    len(pdbs)

.. ipython:: python

    dali_ens = buildPDBEnsemble(pdbs, mapping=mappings, seqid=20)
    dali_ens

Finally we save the ensemble for later processing:

.. ipython:: python

   saveEnsemble(dali_ens, 'dali_ensemble')



.. _`Structure Analysis Tutorial`: http://prody.csb.pitt.edu/tutorials/structure_analysis/blastpdb.html
.. _`list_comprehensions`: https://docs.python.org/2/tutorial/datastructures.html#list-comprehensions

.. [Y13] Shi Y
   Common folds and transport mechanisms of secondary active transporters.
   *Annu. Rev. Biophys.* **2013** 42:51-72
