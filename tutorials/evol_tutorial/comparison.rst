.. _comparison:

Sequence-Structure Comparison
===============================================================================

The part shows how to compare sequence conservation properties with
structural mobility obtained from Gaussian network model (GNM) calculations.

.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()  # turn interactive mode on


Entropy Calculation
-------------------------------------------------------------------------------

First, we retrieve MSA for protein for protein family :pfam:`PF00074`:

.. ipython::
   :verbatim:

   In [1]: fetchPfamMSA('PF00074')
   Out[1]: 'PF00074_full.sth'

We parse the MSA file:

.. ipython:: python

   msa = parseMSA('PF00074_full.sth')

Then, we refine it using :func:`.refineMSA` based on the sequence of
:uniprot:`RNAS1_BOVIN`:

.. ipython:: python

   msa_refine = refineMSA(msa, label='RNAS1_BOVIN', seqid=0.98)


We calculate the entropy for the refined MSA:

.. ipython:: python

   entropy = calcShannonEntropy(msa_refine)


Mobility Calculation
-------------------------------------------------------------------------------

Next, we obtain residue fluctuations or mobility for protein member of the
above family. We will use chain B of :pdb:`2W5I`.

.. ipython:: python

   pdb = parsePDB('2W5I', chain='B')

We can use the function :func:`.alignSequenceToMSA` to identify the part of 
the PDB file that matches with the MSA. In cases where this fails, a label 
keyword argument can be provided to the function as well.

.. ipython:: python

   aln, idx_1, idx_2 = alignSequenceToMSA(pdb, msa_refine, label='RNAS1_BOVIN')
   showAlignment(aln, indices=[idx_1, idx_2])

This tells us that the first two residues are missing as are the last three, ending the 
sequence at residue 121. Hence, we make a selection accordingly::

.. ipython:: python

   chB = pdb.select('resnum 3 to 121')


We can see from the sequence that this gives us the right portion:

.. ipython:: python

   chB.ca.getSequence()

We write this selection to a PDB file for use later, e.g. with evol apps.

.. ipython:: python

   writePDB('2W5IB_3-121.pdb', chB)

We perform GNM as follows:

.. ipython:: python

   gnm = GNM('2W5I')
   gnm.buildKirchhoff(chB.ca)
   gnm.calcModes(n_modes=None)  # calculate all modes

Now, let's obtain residue mobility using the slowest mode, the slowest 8 modes,
and all modes:


.. ipython:: python

   mobility_1 = calcSqFlucts(gnm[0])
   mobility_1to8 = calcSqFlucts(gnm[:8])
   mobility_all = calcSqFlucts(gnm[:])


See :ref:`gnm` for details.

Comparison of mobility and conservation
-------------------------------------------------------------------------------

We use the above data to compare structural mobility and degree of
conservation. We can calculate a correlation coefficient between the two
quantities:

.. ipython:: python

   result = corrcoef(mobility_all, entropy)
   result.round(3)[0,1]

We can plot the two curves simultaneously to visualize the correlation.
We have to scale the values of mobility to display them in the same plot.

Plotting
^^^^^^^^

.. ipython:: python

   indices = range(1,122)
   bar(indices, entropy, width=1.2, color='grey');
   xlim(min(indices)-1, max(indices)+1);
   @savefig entropy_mobility.png width=4in
   plot(indices, mobility_all*(max(entropy)/max(mobility_all)), color='b',
   linewidth=2);


Writing PDB files
-------------------------------------------------------------------------------

We can also write PDB with b-factor column replaced by entropy and mobility
values respectively. We can then load the PDB structure in VMD or PyMol to
see the distribution of entropy and mobility on the structure.

.. ipython:: python

   selprot = chB.copy()
   resindex = selprot.getResindices()
   entropy_prot = [entropy[ind] for ind in resindex]
   mobility_prot = [mobility_all[ind]*10 for ind in resindex]
   selprot.setBetas(entropy_prot)
   writePDB('2W5I_entropy.pdb', selprot)
   selprot.setBetas(mobility_prot)
   writePDB('2W5I_mobility.pdb', selprot)

We can see on the structure just as we could in the bar graph that there is 
some correlation with highly conserved (low entropy) regions having low 
mobility and high entropy regions have higher mobility.