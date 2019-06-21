.. _pca-blast:

Homologous Proteins
===============================================================================

This example shows how to perform PCA of a structural dataset obtained by BLAST
searching the PDB. The protein of interest is :wiki:`cytochrome c` (cyt *c*).
This dataset will contain structures sharing 44% or more sequence identity with
human *cyt c*, i.e. its paralogs and/or orthologs.

A :class:`.PCA` instance that stores the covariance matrix and principal modes that
describe the dominant changes in the dataset will be obtained. :class:`.PCA`
instance and principal modes (:class:`.Mode`) can be used as input to functions
in :mod:`.dynamics` module for further analysis.

Input is amino acid sequence of the protein, a reference PDB identifier,
and some parameters.

Setup
-------------------------------------------------------------------------------

Import ProDy and matplotlib into the current namespace.


.. ipython:: python

   from prody import *
   from pylab import *
   ion()



Name of the protein (a name without a white space is preferred)

.. ipython:: python

   name = 'cyt_c'
   ref_pdb = '1hrc'

In order to perform a BLAST search of the PDB, we will need the amino acid 
sequence of our reference protein.  We could get the FASTA format from the PDB, 
or we could get the sequence from the PDB file itself. A more attractive 
method (to us) is to get the sequence using ProDy.

.. ipython:: python

   ref_prot = parsePDB(ref_pdb)
   ref_hv = ref_prot.getHierView()['A']
   sequence = ref_hv.getSequence()

This is the same as simply using

.. ipython:: python

   sequence = '''GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEE
   TLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE'''

Optionally, a list of PDB files to be excluded from analysis can be provided.
In this case dimeric Cyt c structures are excluded from the analysis. To use
all PDB hits, provide an empty list.

.. ipython:: python

   exclude = ['3nbt', '3nbs']

Parameters
-------------------------------------------------------------------------------

It is sometimes useful to set parameters in variables to use multiple times. 
In this case, we use ``seqid`` for the minimum sequence identity for including 
sequences at both selection of BLAST hits and ensemble building.

.. ipython:: python

   seqid = 44

Blast and download
-------------------------------------------------------------------------------

A list of PDB structures can be obtained using :func:`.blastPDB`
as follows:

.. ipython:: python

   blast_record = blastPDB(sequence)

It is a good practice to save this record on disk, as NCBI may not respond to
repeated searches for the same sequence. We can do this using the Python standard
library :mod:`pickle` as follows:

.. ipython:: python

   import pickle

The record is saved using the :func:`~pickle.dump` function:

.. ipython:: python

   pickle.dump(blast_record, open('cytc_blast_record.pkl', 'wb'))


Then, it can be loaded using the :func:`~pickle.load` function:

.. ipython:: python

   blast_record = pickle.load(open('cytc_blast_record.pkl'))

We then read information from the record to extract a list of 
PDB IDs and chain IDs.

.. ipython:: python

   pdb_hits = []
   for key, item in blast_record.getHits(seqid).iteritems():
       pdb_hits.append((key, item['chain_id']))

Let's parse the PDB files and see how many there are:

.. ipython:: python

   pdbs = parsePDB(*[pdb for pdb, ch in pdb_hits], subset='ca', compressed=False)


.. ipython:: python

   len(pdbs)


Set reference
-------------------------------------------------------------------------------

We first parse the reference structure. Note that we parse only Cα atoms from
chain A. The analysis will be performed for a single chain (monomeric) protein.
For analysis of a dimeric protein see :ref:`pca-dimer`

.. ipython:: python

   reference_structure = parsePDB(ref_pdb, subset='ca', chain='A')
   reference_hierview = reference_structure.getHierView()
   reference_chain = reference_hierview['A']

Prepare ensemble
-------------------------------------------------------------------------------

X-ray structural ensembles are heterogenenous, i.e. different structures
have different sets of unresolved residues. Hence, it is not straightforward
to analyzed them as it would be for NMR models (see :ref:`pca-nmr`).

ProDy has special functions and classes for facilitating efficient analysis
of the PDB X-ray data. In this example we use :func:`.mapOntoChain`
function which returns an :class:`.AtomMap` instance. See :ref:`atommaps` for more details.

The resulting :class:`.AtomMap` instances are used to  prepare a :class:`.PDBEnsemble` 
by mapping each structure against the reference chain and adding a coordinates set 
corresponding to the mapped atoms. The overall procedure is shown in detail below 
so you can understand the process and think about case specific changes such as those in 
the `Multimeric Structures tutorial`_. This process can also be automated using 
:func:`.buildPDBEnsemble` as shown in the `Heterogeneous X-ray Structures tutorial`_.

.. ipython:: python

   startLogfile('pca_blast')
   ensemble = PDBEnsemble(name)
   ensemble.setAtoms(reference_chain)
   ensemble.setCoords(reference_chain.getCoords())

.. ipython:: python

   for structure in pdbs:
       if structure.getTitle()[:4] in exclude:
           continue
       if structure is None:
           plog('Failed to parse ' + pdb_file)
           continue
       mappings = mapOntoChain(structure, reference_chain, seqid=seqid)
       if len(mappings) == 0:
           plog('Failed to map', structure.getTitle()[:4])
           continue
       atommap = mappings[0][0]
       ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
   ensemble.iterpose()
   saveEnsemble(ensemble)


Let's check how many conformations are extracted from PDB files:

.. ipython:: python

   len(ensemble)

Note that the number of conformations is larger than the number of PDB structures
we retrieved. This is because some of the PDB files contained NMR structures
with multiple models. Each model in NMR structures are added to the ensemble
as individual conformations.

Write aligned conformations into a PDB file as follows:

.. ipython:: python

   writePDB(name+'.pdb', ensemble)


This file can be used to visualize the aligned conformations in a modeling
software.



Align PDB files
-------------------------------------------------------------------------------

:func:`.alignPDBEnsemble` function can be used to align PDB structures used
in the analysis and write new PDB files, e.g. ``alignPDBEnsemble(ensemble)``. 
The resulting files will contain intact structures and can be used for 
visualization purposes. In this case, we will align only select PDB files:

.. ipython:: python

   conf1_aligned = alignPDBEnsemble(ensemble[0])
   conf2_aligned = alignPDBEnsemble(ensemble[1])


Let's take a quick look at the aligned structures:

.. ipython:: python


   showProtein(parsePDB(conf1_aligned), parsePDB(conf2_aligned));
   @savefig ensemble_analysis_blast_aligned.png width=4in
   legend();


Perform PCA
-------------------------------------------------------------------------------

Once the ensemble is ready, performing PCA is 3 easy steps:

.. ipython:: python

   pca = PCA(name)
   pca.buildCovariance(ensemble)
   pca.calcModes()

The calculated data can be saved as a compressed file using :func:`.saveModel`
function:

.. ipython:: python

   saveModel(pca)


Plot results
-------------------------------------------------------------------------------


Let's plot RMSDs of all conformations from the average conformation:


.. ipython:: python

   rmsd = calcRMSD(ensemble)
   plot(rmsd);
   xlabel('Conformation index');
   @savefig ensemble_analysis_blast_rmsd.png width=4in
   ylabel('RMSD (A)');


Let's show a projection of the ensemble onto PC1 and PC2:

.. ipython:: python

   @savefig ensemble_analysis_blast_projection.png width=4in
   showProjection(ensemble, pca[:2]);

.. _`Multimeric Structures tutorial`: http://prody.csb.pitt.edu/tutorials/ensemble_analysis/dimer.html
.. _`Heterogeneous X-ray Structures tutorial`: http://prody.csb.pitt.edu/tutorials/ensemble_analysis/xray_calculations.html