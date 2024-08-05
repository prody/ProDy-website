.. _watfinder_parameters:

Changes in the default parameters
===============================================================================

There are a lot of parameters that can be changed in the water bridge
analysis, including the distances, angles, and the number of involved water
molecules or residues.

``atoms``: Atomic object from which atoms are considered

``method``: ``'cluster'`` or ``'chain'``, where 'chain' will find the shortest 
    water bridging path between two protein atoms
    default is ``'chain'``

``distDA``: maximal distance between water/protein donor and acceptor
    default is ``3.5``

``distWR``: maximal distance between considered water and any residue
    default is ``4``

``anglePDWA``: angle range where protein is donor and water is acceptor
    default is ``(100, 200)``

``anglePAWD``: angle range where protein is acceptor and water is donor
    default is ``(100, 140)``

``angleWW``: angle between water donor/acceptor
    default is ``(140, 180)``

``maxDepth``: maximum number of waters in chain/depth of residues in cluster
    default is ``2``

``maxNumRes``: maximum number of water+protein residues in cluster
    default is ``None``

``donors``: which atoms to count as donors 
    default is ``['N', 'O', 'S', 'F']``

``acceptors``: which atoms to count as acceptors 
    default is ``['N', 'O', 'S', 'F']``

``output``: return information arrays, (protein atoms, water atoms), 
    or just atom indices per bridge
    Options are ``'info'``, ``'atomic'`` and ``'indices'``
    Default is ``'atomic'``

``isInfoLog``: whether to log information
    default is ``True``


The above criteria are used to find hydrogen bonds in a PDB structure which
are used to predict possible water bridges. 

The angles are calculated so that the hydrogen atom that makes the hydrogen
bond is the vertex of the angle, where the donor and acceptor atoms form the
two ends. The angle is measured between the line connecting the donor to the
hydrogen and the line connecting the hydrogen to the acceptor. 

Most PDB structures do not have hydrogen atoms, therefore we can ignore this
criterium. Hydrogen bonds will be calculated just based on distance and atom types. 


Here's an example of how to apply changes in parameters:

.. ipython:: python
   :verbatim:

   waterBridges_2 = calcWaterBridges(atoms, method='cluster', distWR=3.0, distDA=3.2, maxDepth=3)

.. parsed-literal::

   @> 17 water bridges detected.
   @> SER7 OG_21 A LYS110 NZ_838 A 4.70722699686344 1 ['A_1316']
   @> GLY14 O_72 A SER47 O_328 A GLU50 N_347 A 5.288116583434976 4.544135231262378 5.10489901956934 1 ['A_1274']
   @> PRO20 O_115 A GLU23 OE1_139 A 4.571172934816621 1 ['A_1292']
   @> GLU23 OE2_140 A SER71 O_514 A HIS72 ND1_523 A 4.934310286149422 4.127239392136104 4.463734423103597 1 ['A_1244']
   @> ASP42 OD2_301 A ARG40 NH2_286 A 5.163938516287738 1 ['A_1246']
   @> ASP81 OD1_598 A THR84 OG1_621 A ARG40 NH2_286 A 4.415462942886057 4.365525627000715 3.9229717052255175 1 ['A_1262']
   @> GLY14 O_72 A SER47 O_328 A GLU50 N_347 A 5.288116583434976 4.544135231262378 5.10489901956934 1 ['A_1274']
   @> GLU50 OE2_355 A TYR131 OH_1009 A 5.157987010452818 1 ['A_1299']
   @> ARG65 NH2_473 A ASP135 O_1037 A 4.820906553751068 1 ['A_1267']
   @> GLU23 OE2_140 A SER71 O_514 A HIS72 ND1_523 A 4.934310286149422 4.127239392136104 4.463734423103597 1 ['A_1244']
   @> GLN105 NE2_800 A LYS102 O_772 A LYS79 NZ_582 A 4.2363221076778395 3.8752885053889865 4.929221236666094 1 ['A_1249']
   @> ASP81 OD1_598 A THR84 OG1_621 A ARG40 NH2_286 A 4.415462942886057 4.365525627000715 3.9229717052255175 1 ['A_1262']
   @> GLN105 NE2_800 A LYS102 O_772 A LYS79 NZ_582 A 4.2363221076778395 3.8752885053889865 4.929221236666094 1 ['A_1249']
   @> SER7 OG_21 A LYS110 NZ_838 A 4.70722699686344 1 ['A_1316']
   @> LEU125 O_953 A GLY117 N_886 A 3.8595291163560352 1 ['A_1264']
   @> GLU50 OE2_355 A TYR131 OH_1009 A 5.157987010452818 1 ['A_1299']
   @> ARG65 NH2_473 A ASP135 O_1037 A 4.820906553751068 1 ['A_1267']
