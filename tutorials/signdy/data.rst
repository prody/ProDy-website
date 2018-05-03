.. _signdy-data:

Data Collection
===============================================================================

The first step in signature dynamics analysis is to collect a set of related 
protein structures. This can be achieved by multiple routes. In this example, belonging to three 
main types: a query search of the PDB using BLAST or Dali, extraction of PDB IDs 
from the Pfam or CATH database, or input of a pre-defined list. We demonstrate the 
Dali method here. The Pfam and CATH methods are still under development and the 
function :func:`blastPDB` is described in the `Structure Analysis Tutorial`_.

We apply these methods to the type-I periplasmic binding protein domains, 
a group of protein structures originally found in bacteria for transport of solutes 
across the periplasmic space and later seen in various eukaryotic receptors including 
ionotropic and metabotropic glutamate receptors. We use the N-terminal domain of AMPA
receptor subunit GluA2 (gene name GRIA2) as a query.

First, make necessary imports from ProDy and Matplotlib packages if you haven't already.

.. ipython:: python

    from prody import *
    from pylab import *
    ion()

Dali
-------------------------------------------------------------------------------

First we use the function searchDali to search the PDB, which returns a DaliRecord 
that contains a list of PDBs and their corresponding alignments to the reference structure. 
This can be used to build a PDB ensemble directly using the method associated to this class.

Dali 

.. ipython:: python

    dali_rec = searchDali('3H5V','A')
    dali_ens = dali_rec.buildDaliEnsemble()


.. _`Structure Analysis Tutorial`: http://prody.csb.pitt.edu/tutorials/structure_analysis/blastpdb.html
.. _`list_comprehensions`: https://docs.python.org/2/tutorial/datastructures.html#list-comprehensions