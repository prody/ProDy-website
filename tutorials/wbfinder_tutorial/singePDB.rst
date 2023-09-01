.. _wbfinder_tutorial:

Water bridges detection in a single PDB structure
===============================================================================


Water bridges prediction
-------------------------------------------------------------------------------

To analyze the structure we need to parse a structure **addH_5kqm.pdb**
using :func:`.parsePDB`:


.. ipython:: python

   PDBfile = 'addH_'+PDB+'.pdb'
   coords = parsePDB(PDBfile)


Before analysis, we can check how many water molecules are present in our PDB 
structure and later compare how many of them were meaningful for protein structure:


.. ipython:: python

   water_molecules = coords.select('water')
   len(water_molecules)


Next, we can use :func:`.calcWaterBridges` and one out of two methods to detect 
water bridges, *'chain'* or *'cluster'*:


1. **Method 'chain' (default)** which will detect water molecules between pairs of 
hydrophilic residues:


.. ipython:: python

   waterBridges_chain = calcWaterBridges(coords)


2. **Method 'cluster'** which will detect water molecules between multiple hydrophilic 
residues:


.. ipython:: python

   waterBridges_cluster = calcWaterBridges(coords, method='cluster')


*'Chain' method* detected **42** water bridges and *'Cluster' method* second **49**. 
The total number of water molecules in the crystal structure is **363**. As we can 
see, many of them are not meaningful for the protein stability.



Save results in PDB file
-------------------------------------------------------------------------------

We can use :func:`.savePDBWaterBridges` to save the results in PDB file. File will 
contain water molecules that are forming potential hydrogen bridges and protein 
structure. Residues involved in water bridges can be displayed using occupancy column.


.. ipython:: python

   savePDBWaterBridges(waterBridges_cluster, coords, PDBfile[:-4]+'_wb_cluster.pdb')
   savePDBWaterBridges(waterBridges_chain, coords, PDBfile[:-4]+'_wb_chain.pdb')


The results can be displayed in VMD. Below we can see a comparison between
results obtained by 'chain' vs 'cluster' (additional molecules shown in
green) method.


.. figure:: images/Fig1.tga
   :scale: 60 %


Access to the raw data
-------------------------------------------------------------------------------

To have acces to the raw data, we need to include paramater 
additional parameter *ouput='info'* in :func:`.calcWaterBridges`.


.. ipython:: python

   waterBridges_cluster = calcWaterBridges(coords, method='cluster', output='info')
   waterBridges_cluster


.. ipython:: python

   waterBridges_chain = calcWaterBridges(coords, output='info')


We can check which residues are involved in water bridges using the code below. 
First we need to extract residues names and display them without repetitions.


.. ipython:: python

   allresidues = []
   
   for i in waterBridges_chain:
       allresidues.append(i[0])
       allresidues.append(i[3])

   import numpy as np
   allresidues_once = np.unique(allresidues)    
   allresidues_once


We can also count how many times each residue was involved in water bridges 
(with different waters) and display the number of counts as a histogram.


.. ipython:: python

   from collections import Counter
   aa_counter = Counter(allresidues)
   sorted_aa_counter = dict(sorted(aa_counter.items(), key=lambda item: item[1], reverse=True))
   sorted_aa_counter


.. ipython:: python

   import matplotlib.pyplot as plt

   values = list(sorted_aa_counter.values())
   labels = list(sorted_aa_counter.keys())

   plt.figure(figsize=(10, 4))
   plt.bar(labels, values)
   plt.xticks(rotation=90)
   plt.xlabel('Residues')
   plt.ylabel('#')
   plt.tight_layout()
   plt.show()


Based on the results we can see that there is one residue, GLU23, which 
participate often in the interactions with water molecules.

