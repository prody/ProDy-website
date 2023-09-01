.. _wbfinder_tutorial:

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

   waterBridges_2 = calcWaterBridges(coords, method='cluster', distWR=3.0, distDA=3.2, maxDepth=3)

