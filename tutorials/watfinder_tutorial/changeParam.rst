.. _watfinder_tutorial:

Changes in the default parameters
===============================================================================

There are a lot of paramaters that can be changed in the water bridges analysis 
including the distances, angles and the number of involved water molecules or 
residues.

:arg **atoms**: Atomic object from which atoms are considered

:arg **method**: cluster or chain, where chain find shortest water bridging 
    path between two protein atoms
    default is 'chain'

:arg **distDA**: maximal distance between water/protein donor and acceptor
    default is 3.5

:arg **distWR**: maximal distance between considered water and any residue
    default is 4

:arg **anglePDWA**: angle range where protein is donor and water is acceptor
    default is (100, 200)

:arg **anglePAWD**: angle range where protein is acceptor and water is donor
    default is (100, 140)

:arg **angleWW**: angle between water donor/acceptor
    default is (140, 180)

:arg **maxDepth**: maximum number of waters in chain/depth of residues in cluster
    default is 2

:arg **maxNumRes**: maximum number of water+protein residues in cluster
    default is None

:arg **donors**: which atoms to count as donors 
    default is ['N', 'O', 'S', 'F']

:arg **acceptors**: which atoms to count as acceptors 
    default is ['N', 'O', 'S', 'F']

:arg **output**: return information arrays, (protein atoms, water atoms), 
    or just atom indices per bridge
    default is 'atomic'
:type output: 'info' | 'atomic' | 'indices'

:arg **isInfoLog**: should log information
    default is True


.. ipython:: python
   :verbatim:

   waterBridges_2 = calcWaterBridges(coords, method='cluster', distWR=3.0, distDA=3.2, maxDepth=3)

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
