.. _watfinder_tutorial:
=======

Water bridges detection in a trajectory
===============================================================================

Now, we will perform calculations for a trajectory file. We will use 
:func:`.calcWaterBridgesTrajectory` for which we need to provide PDB and DCD file.

We will use files that were prepared in NAMD. The system (protein in a water box) 
can be found in *5kqm_all_sci.pdb*. Trajectory, *NAMD_D2_sample.dcd*, has dcd format.


Parse structure with trajectory
-------------------------------------------------------------------------------


.. ipython:: python
   :verbatim:

   PDBtraj_file = "5kqm_all_sci.pdb"
   coords_traj = parsePDB(PDBtraj_file)
   trajectory = parseDCD("NAMD_D2_sample.dcd")

.. parsed-literal::

   @> 19321 atoms and 1 coordinate set(s) were parsed in 0.18s.
   @> DCD file contains 17 coordinate sets for 19321 atoms.
   @> DCD file was parsed in 0.01 seconds.
   @> 3.76 MB parsed at input rate 748.06 MB/s.
   @> 17 coordinate sets parsed at input rate 3382 frame/s.

The analysis od water bridges can be performed on selected frames by using 
*start_frame* or *stop_frame*. 

.. ipython:: python
   :verbatim:

   wb_traj = calcWaterBridgesTrajectory(coords_traj, trajectory, start_frame=5, stop_frame=15, output='info')

.. parsed-literal::

   @> Frame: 5
   @> 101 water bridges detected.
   @> Frame: 6
   @> 107 water bridges detected.
   @> Frame: 7
   @> 90 water bridges detected.
   @> Frame: 8
   @> 97 water bridges detected.
   @> Frame: 9
   @> 122 water bridges detected.
   @> Frame: 10
   @> 101 water bridges detected.
   @> Frame: 11
   @> 130 water bridges detected.
   @> Frame: 12
   @> 132 water bridges detected.
   @> Frame: 13
   @> 126 water bridges detected.
   @> Frame: 14
   @> 88 water bridges detected.
   @> Frame: 15
   @> 105 water bridges detected.

Because of the number of data there results will not be displayed. We can have 
an access to the raw data by using *output='info'*.

.. ipython:: python
   :verbatim:

   wb_traj

.. parsed-literal::

   [[['THR5',
      'OG1_8',
      'P',
      'TRP39',
      'NE1_547',
      'P',
      3.261452627054999,
      1,
      ['3W_13313']],
     ['THR5',
      'O_15',
      'P',
      'ASP86',
      'OD1_1269',
      'P',
      5.986350034086454,
      2,
      ['3W_12974', '3W_18431']],
     ['THR5',
      'O_15',
      'P',
      'LYS110',
      'NZ_1667',
      'P',
      7.375256709599827,
      2,
      ['3W_12974', '3W_18431']],
     ['THR5',
      'O_15',
      'P',
      'LYS6',
      'NZ_32',
      'P',
      6.414308925017051,
      2,
      ['3W_12974', '3W_12152']],
     ['LYS6',
      'NZ_32',
      'P',
      'TYR87',
      'OH_1286',
      'P',
      4.891713264838611,
      1,
      ['3W_9209']],
     ['LYS6',
      'NZ_32',
      'P',
      'ASP86',
      'O_1272',
      'P',
      6.000079664458025,
      2,
      ['3W_7319', '3W_17114']],
     ['GLY14',
      'O_156',
      'P',
      'GLU50',
      'O_711',
      'P',
      4.4703701847403154,
      1,
      ['3W_14210']],
     ['GLY14',
      'O_156',
      'P',
      'ASN53',
      'ND2_747',
      'P',
      5.847016041542153,
      1,
      ['3W_14210']],
     ['ARG18',
      'NH1_217',
      'P',
      'ASP129',
      'N_1970',
      'P',
      4.176516917230124,
      1,
      ['3W_9143']],
     ['ARG18',
      'NH2_220',
      'P',
      'ASN95',
      'OD1_1406',
      'P',
      3.3808444537964073,
      1,
      ['3W_9365']],
     ['ARG18',
      'NH2_220',
      'P',
      'ASN95',
      'ND2_1407',
      'P',
      4.564229284186745,
      1,
      ['3W_9365']],
     ['ALA24',
      'O_303',
      'P',
      'LYS28',
      'N_364',
      'P',
      2.890261036024075,
      1,
      ['3W_6209']],
     ['ARG27',
      'NE_353',
      'P',
      'ARG27',
      'NH2_359',
      'P',
      2.238207059787312,
      1,
      ['3W_18434']],
     ['ARG27',
      'NE_353',
      'P',
      'VAL41',
      'N_585',
      'P',
      6.508044630925681,
      2,
      ['3W_18434', '3W_6350']],
     ['ARG27',
      'NH2_359',
      'P',
      'VAL41',
      'N_585',
      'P',
      5.304460459374652,
      2,
      ['3W_18434', '3W_6350']],
     ['LYS28',
      'NZ_380',
      'P',
      'TYR142',
      'OH_2172',
      'P',
      7.115418987200776,
      2,
      ['3W_5423', '3W_10319']],
     ['SER36',
      'N_497',
      'P',
      'GLU37',
      'N_508',
      'P',
      2.8466821983592507,
      1,
      ['3W_13610']],
     ['GLU37',
      'OE2_520',
      'P',
      'ASN38',
      'ND2_532',
      'P',
      5.9092287655354365,
      2,
      ['3W_9515', '3W_9869']],
     ['GLU37',
      'O_522',
      'P',
      'TRP39',
      'O_560',
      'P',
      5.639727451637652,
      2,
      ['3W_14753', '3W_9479']],
     ['TRP39',
      'O_560',
      'P',
      'VAL41',
      'N_585',
      'P',
      3.686292263336709,
      2,
      ['3W_9479', '3W_6350']],
     ['ARG40',
      'NE_574',
      'P',
      'ARG40',
      'NH2_580',
      'P',
      2.284446592624093,
      1,
      ['3W_13502']],
      ..
      ..


Save the results
-------------------------------------------------------------------------------

The results can be saved using :func:`.saveWaterBridges` in two formats. Txt 
file will contain all the results for analysis and can be visualized in text 
editor, and wb file will restore data for further analysis. It can be uploaded 
using :func:`.parseWaterBridges` as shown below.


.. ipython:: python
   :verbatim:

   waterBridges_save = calcWaterBridgesTrajectory(coords_traj, trajectory, stop_frame=15)
   saveWaterBridges(waterBridges_save,'wb_saved.txt')
   saveWaterBridges(waterBridges_save,'wb_saved.wb')

.. parsed-literal::

   @> Frame: 0
   @> 48 water bridges detected.
   @> Frame: 1
   @> 48 water bridges detected.
   @> Frame: 2
   @> 94 water bridges detected.
   @> Frame: 3
   @> 110 water bridges detected.
   @> Frame: 4
   @> 105 water bridges detected.
   @> Frame: 5
   @> 101 water bridges detected.
   @> Frame: 6
   @> 107 water bridges detected.
   @> Frame: 7
   @> 90 water bridges detected.
   @> Frame: 8
   @> 97 water bridges detected.
   @> Frame: 9
   @> 122 water bridges detected.
   @> Frame: 10
   @> 101 water bridges detected.
   @> Frame: 11
   @> 130 water bridges detected.
   @> Frame: 12
   @> 132 water bridges detected.
   @> Frame: 13
   @> 126 water bridges detected.
   @> Frame: 14
   @> 88 water bridges detected.
   @> Frame: 15
   @> 105 water bridges detected.

To upload wb file use :func:`.parseWaterBridges` and protein coordinates 
as follows:

.. ipython:: python
   :verbatim:

   uploaded_results = parseWaterBridges('wb_saved.wb',coords_traj)

Uploaded results are of type atomic (.wb file), therefore it can be used for 
analysis later. The atomic output can be transformed to 
detailed information using :func:`.getWaterBridgesInfoOutput`.


Analysis of the results
-------------------------------------------------------------------------------

To perform analysis we can not use *output='info'*, therefore we will run 
the calculations again. This time we will run first 15 frames of the simulation.


.. ipython:: python
   :verbatim:

   waterBridges = calcWaterBridgesTrajectory(coords_traj, trajectory, stop_frame=15)

.. parsed-literal::

   @> Frame: 0
   @> 48 water bridges detected.
   @> Frame: 1
   @> 48 water bridges detected.
   @> Frame: 2
   @> 94 water bridges detected.
   @> Frame: 3
   @> 110 water bridges detected.
   @> Frame: 4
   @> 105 water bridges detected.
   @> Frame: 5
   @> 101 water bridges detected.
   @> Frame: 6
   @> 107 water bridges detected.
   @> Frame: 7
   @> 90 water bridges detected.
   @> Frame: 8
   @> 97 water bridges detected.
   @> Frame: 9
   @> 122 water bridges detected.
   @> Frame: 10
   @> 101 water bridges detected.
   @> Frame: 11
   @> 130 water bridges detected.
   @> Frame: 12
   @> 132 water bridges detected.
   @> Frame: 13
   @> 126 water bridges detected.
   @> Frame: 14
   @> 88 water bridges detected.
   @> Frame: 15
   @> 105 water bridges detected.

Information about residues contributiong to water bridges
-------------------------------------------------------------------------------

Analysis of the data can be performed using :func:`.calcWaterBridgesStatistics`.
The analysis presented below gave information about pairs of residues involved 
in water bridges, how often they occur, and the average distance between them. 
Standard deviation provides information on how the distance was changing during 
the simulation.Additionally, the analysis can be saved by using *filename* option.


.. ipython:: python
   :verbatim:
   
   analysisAtomic = calcWaterBridgesStatistics(waterBridges, trajectory, filename='data.txt')

   for item in analysisAtomic.items():
      print(item)

.. parsed-literal::

   @> RES1           RES2           PERC      DIST_AVG  DIST_STD
   @> SER7P          ARG40P         12.500    4.901     0.000
   @> ARG18P         ASP92P         68.750    4.285     1.159
   @> ARG18P         ASN95P         68.750    5.099     1.192
   @> PRO20P         GLU23P         12.500    4.571     0.000
   @> GLU23P         HSE72P         12.500    3.669     0.458
   @> ARG27P         VAL41P         56.250    5.565     0.781
   @> ARG27P         SER71P         75.000    6.116     0.445
   @> ASP32P         ASN34P         25.000    4.218     0.652
   @> SER36P         GLU37P         75.000    3.700     1.154
   @> ARG40P         THR84P         50.000    4.235     0.671
   @> ASP42P         ARG75P         68.750    3.159     0.652
   @> THR46P         ASN95P         62.500    4.067     0.842
   @> SER47P         TYR49P         50.000    4.320     0.757
   @> SER47P         GLU50P         37.500    6.161     1.070
   @> TYR49P         GLU50P         68.750    3.194     0.730
   @> PRO55P         ALA74P         18.750    5.998     0.070
   @> ARG58P         ARG65P         12.500    6.682     0.000
   @> ARG58P         ASP135P        18.750    6.446     0.236
   @> ARG65P         ASP135P        50.000    5.664     0.955
   @> ARG75P         GLN76P         31.250    3.659     0.871
   @> LYS79P         LYS102P        25.000    5.724     1.117
   @> LYS79P         GLN105P        50.000    5.392     1.134
   @> CYS90P         MET91P         25.000    3.745     0.246
   @> CYS90P         GLY117P        37.500    4.888     0.232
   @> CYS90P         LEU125P        18.750    5.984     0.285
   @> MET91P         GLY117P        37.500    5.674     0.242
   @> MET91P         LEU125P        25.000    4.293     0.315
   @> ASP92P         ASN95P         62.500    4.548     1.036
   @> GLU93P         SER94P         43.750    2.784     0.092
   @> ASN100P        ILE113P        25.000    6.672     1.103
   @> LYS102P        GLN105P        50.000    5.948     1.507
   @> GLY117P        LEU125P        25.000    3.841     0.064
   @> LEU125P        ILE126P        81.250    3.264     0.744
   @> ASP137P        THR140P        12.500    5.251     0.000
   @> GLN143P        ARG147P        37.500    5.429     0.843
   @> THR5P          ASN38P         56.250    5.146     1.287
   @> THR5P          LYS6P          43.750    4.728     0.954
   @> LYS6P          ALA156P        18.750    7.223     0.266
   @> GLU23P         ARG27P         43.750    3.505     0.675
   @> GLU23P         SER71P         18.750    7.137     0.716
   @> LYS28P         TYR142P        18.750    7.302     0.262
   @> THR31P         SER36P         18.750    4.212     0.904
   @> ASN34P         ASN38P         6.250     6.497     0.000
   @> SER36P         TRP39P         31.250    4.077     1.199
   @> GLU37P         ASN38P         37.500    4.182     1.106
   @> ARG40P         ARG75P         37.500    6.434     0.768
   @> SER47P         ARG75P         6.250     7.044     0.000
   @> TYR49P         ARG75P         6.250     7.837     0.000
   @> TYR49P         ASN53P         25.000    6.312     1.235
   @> GLU50P         ARG75P         6.250     7.079     0.000
   @> ILE51P         ASN53P         12.500    4.529     0.793
   @> TYR57P         GLN60P         31.250    5.961     0.978
   @> ARG58P         GLY133P        37.500    3.202     0.549
   @> ARG58P         ASN134P        6.250     7.308     0.000
   @> SER61P         ARG65P         62.500    5.253     0.722
   ..
   ..


To have an access to the data we can use :func:`.getWaterBridgeStatInfo`.

.. ipython:: python
   :verbatim:
   
   getWaterBridgeStatInfo(analysisAtomic, coords_traj)

.. parsed-literal::

   {('SER7P', 'ARG40P'): {'percentage': 12.5,
     'distAvg': 4.9006157,
     'distStd': 0.0},
    ('ARG18P', 'ASP92P'): {'percentage': 68.75,
     'distAvg': 4.2853837,
     'distStd': 1.159262},
    ('ARG18P', 'ASN95P'): {'percentage': 68.75,
     'distAvg': 5.0986476,
     'distStd': 1.1916962},
    ('PRO20P', 'GLU23P'): {'percentage': 12.5,
     'distAvg': 4.571081,
     'distStd': 0.0},
    ('GLU23P', 'HSE72P'): {'percentage': 12.5,
     'distAvg': 3.668869,
     'distStd': 0.45773232},
    ('ARG27P', 'VAL41P'): {'percentage': 56.25,
     'distAvg': 5.5646605,
     'distStd': 0.7812058},
    ('ARG27P', 'SER71P'): {'percentage': 75.0,
     'distAvg': 6.1158614,
     'distStd': 0.4451057},
    ('ASP32P', 'ASN34P'): {'percentage': 25.0,
     'distAvg': 4.217685,
     'distStd': 0.65174913},
    ('SER36P', 'GLU37P'): {'percentage': 75.0,
     'distAvg': 3.7002723,
     'distStd': 1.1539057},
    ('ARG40P', 'THR84P'): {'percentage': 50.0,
     'distAvg': 4.2351923,
     'distStd': 0.67144614},
    ('ASP42P', 'ARG75P'): {'percentage': 68.75,
     'distAvg': 3.1586912,
     'distStd': 0.65187186},
    ('THR46P', 'ASN95P'): {'percentage': 62.5,
     'distAvg': 4.067392,
     'distStd': 0.84244806},
    ('SER47P', 'TYR49P'): {'percentage': 50.0,
     'distAvg': 4.3195844,
     'distStd': 0.7574061},
    ('SER47P', 'GLU50P'): {'percentage': 37.5,
     'distAvg': 6.160927,
     'distStd': 1.0696731},
    ('TYR49P', 'GLU50P'): {'percentage': 68.75,
     'distAvg': 3.1939995,
     'distStd': 0.7295457},
    ('PRO55P', 'ALA74P'): {'percentage': 18.75,
     'distAvg': 5.9978576,
     'distStd': 0.0696475},
    ('ARG58P', 'ARG65P'): {'percentage': 12.5,
     'distAvg': 6.6822205,
     'distStd': 0.0},
    ('ARG58P', 'ASP135P'): {'percentage': 18.75,
     'distAvg': 6.4458184,
     'distStd': 0.23627748},
    ('ARG65P', 'ASP135P'): {'percentage': 50.0,
     'distAvg': 5.663641,
     'distStd': 0.9553071},
    ('ARG75P', 'GLN76P'): {'percentage': 31.25,
     'distAvg': 3.6593163,
     'distStd': 0.871043},
    ('LYS79P', 'LYS102P'): {'percentage': 25.0,
     'distAvg': 5.7237897,
     'distStd': 1.1169385},
    ('LYS79P', 'GLN105P'): {'percentage': 50.0,
     'distAvg': 5.3919144,
     'distStd': 1.133555},
    ('CYS90P', 'MET91P'): {'percentage': 25.0,
     'distAvg': 3.744512,
     'distStd': 0.2455233},
    ('CYS90P', 'GLY117P'): {'percentage': 37.5,
     'distAvg': 4.887792,
     'distStd': 0.2316904},
    ('CYS90P', 'LEU125P'): {'percentage': 18.75,
     'distAvg': 5.984488,
     'distStd': 0.28450695},
    ('MET91P', 'GLY117P'): {'percentage': 37.5,
     'distAvg': 5.673981,
     'distStd': 0.24238132},
    ('MET91P', 'LEU125P'): {'percentage': 25.0,
     'distAvg': 4.2932153,
     'distStd': 0.31463936},
    ('ASP92P', 'ASN95P'): {'percentage': 62.5,
     'distAvg': 4.5477376,
     'distStd': 1.0361063},
    ('GLU93P', 'SER94P'): {'percentage': 43.75,
     'distAvg': 2.7842126,
     'distStd': 0.09215464},
    ('ASN100P', 'ILE113P'): {'percentage': 25.0,
     'distAvg': 6.6723633,
     'distStd': 1.1034071},
    ('LYS102P', 'GLN105P'): {'percentage': 50.0,
     'distAvg': 5.948009,
     'distStd': 1.5065572},
    ('GLY117P', 'LEU125P'): {'percentage': 25.0,
     'distAvg': 3.8412733,
     'distStd': 0.06422541},
    ('LEU125P', 'ILE126P'): {'percentage': 81.25,
     'distAvg': 3.2639554,
     'distStd': 0.7435347},
    ('ASP137P', 'THR140P'): {'percentage': 12.5,
     'distAvg': 5.251236,
     'distStd': 0.0},
    ('GLN143P', 'ARG147P'): {'percentage': 37.5,
     'distAvg': 5.4292507,
     'distStd': 0.8427817},
    ('THR5P', 'ASN38P'): {'percentage': 56.25,
     'distAvg': 5.1459327,
     'distStd': 1.287013},
    ('THR5P', 'LYS6P'): {'percentage': 43.75,
     'distAvg': 4.7280626,
     'distStd': 0.954359},
    ('LYS6P', 'ALA156P'): {'percentage': 18.75,
     'distAvg': 7.2233906,
     'distStd': 0.26610205},
    ('GLU23P', 'ARG27P'): {'percentage': 43.75,
     'distAvg': 3.5054235,
     'distStd': 0.6751442},
    ('GLU23P', 'SER71P'): {'percentage': 18.75,
     'distAvg': 7.136764,
     'distStd': 0.7162538},
    ('LYS28P', 'TYR142P'): {'percentage': 18.75,
     'distAvg': 7.3022885,
     'distStd': 0.26167196},
    ('THR31P', 'SER36P'): {'percentage': 18.75,
     'distAvg': 4.211753,
     'distStd': 0.9039756},
    ('ASN34P', 'ASN38P'): {'percentage': 6.25,
     'distAvg': 6.497203,
     'distStd': 0.0},
      ..
      ..

To obtain maps of interactions for protein structure, we can use 
:func:`.showWaterBridgeMatrix` which is equipted in three paramaters: 
*'percentage'* (how often two residues were forming water bridges), 
*'distAvg'* (how close there were), and *'distStd'* (how stable that 
water bridge was).


.. ipython:: python
   :verbatim:
   
   showWaterBridgeMatrix(analysisAtomic, 'percentage')

.. figure:: images/traj_percentage.png
   :scale: 60 %

.. ipython:: python
   :verbatim:
   
   showWaterBridgeMatrix(analysisAtomic, 'distAvg')

.. figure:: images/traj_distAvg.png
   :scale: 60 %   

.. ipython:: python
   :verbatim:   

   showWaterBridgeMatrix(analysisAtomic, 'distStd')

.. figure:: images/traj_distStd.png
   :scale: 60 %

Raw data of the matrices can be obatined with :func:`.calcWaterBridgeMatrix`. 
The type of the matrix can be selected among: *'percentage'*, *'distAvg'*, *'distStd'*.


.. ipython:: python
   :verbatim:

    M1 = calcWaterBridgeMatrix(analysisAtomic, 'percentage')
    M2 = calcWaterBridgeMatrix(analysisAtomic, 'distAvg')
    M3 = calcWaterBridgeMatrix(analysisAtomic, 'distStd')

.. ipython:: python
   :verbatim:

   M1

.. parsed-literal::

   array([[ 0.  ,  0.  ,  0.  , ...,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  , ...,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  , ...,  0.  ,  0.  ,  0.  ],
          ...,
          [ 0.  ,  0.  ,  0.  , ...,  0.  , 12.5 , 31.25],
          [ 0.  ,  0.  ,  0.  , ..., 12.5 ,  0.  , 12.5 ],
          [ 0.  ,  0.  ,  0.  , ..., 31.25, 12.5 ,  0.  ]])

.. ipython:: python
   :verbatim:

   M2

.. parsed-literal::

   array([[0.        , 0.        , 0.        , ..., 0.        , 0.        ,
           0.        ],
          [0.        , 0.        , 0.        , ..., 0.        , 0.        ,
           0.        ],
          [0.        , 0.        , 0.        , ..., 0.        , 0.        ,
           0.        ],
          ...,
          [0.        , 0.        , 0.        , ..., 0.        , 4.58851337,
           5.82083416],
          [0.        , 0.        , 0.        , ..., 4.58851337, 0.        ,
           3.52366138],
          [0.        , 0.        , 0.        , ..., 5.82083416, 3.52366138,
           0.        ]])


.. ipython:: python
   :verbatim:

   M3

.. parsed-literal::

   array([[0.        , 0.        , 0.        , ..., 0.        , 0.        ,
           0.        ],
          [0.        , 0.        , 0.        , ..., 0.        , 0.        ,
           0.        ],
          [0.        , 0.        , 0.        , ..., 0.        , 0.        ,
           0.        ],
          ...,
          [0.        , 0.        , 0.        , ..., 0.        , 1.71697354,
           1.38650537],
          [0.        , 0.        , 0.        , ..., 1.71697354, 0.        ,
           1.27207112],
          [0.        , 0.        , 0.        , ..., 1.38650537, 1.27207112,
           0.        ]])


Statistical analysis for water bridges
-------------------------------------------------------------------------------

To visualize the results in a more accessible way, we can use 
:func:`.calcWaterBridgeMatrix` function which will show how often each residue 
were contributing to the water bridges in the trajectory.


.. ipython:: python
   :verbatim:

   calcBridgingResiduesHistogram(waterBridges)

.. figure:: images/traj_res_hist.png
   :scale: 60 %

.. parsed-literal::

   [('LEU96P', 1),
    ('MET63P', 1),
    ('PHE152P', 1),
    ('LEU29P', 1),
    ('PRO130P', 1),
    ('PHE85P', 1),
    ('PRO54P', 1),
    ('ILE16P', 1),
    ('CYS148P', 1),
    ('VAL25P', 1),
    ('ILE77P', 1),
    ('PRO20P', 2),
    ('ILE127P', 2),
    ('ILE68P', 2),
    ('GLY14P', 2),
    ('GLY67P', 2),
    ('ALA111P', 3),
    ('VAL73P', 3),
    ('ALA24P', 3),
    ('LEU115P', 3),
    ('PRO55P', 4),
    ('ALA74P', 4),
    ('PRO121P', 4),
    ('ASN15P', 4),
    ('LEU13P', 4),
    ('ILE51P', 5),
    ('THR31P', 5),
    ('TYR119P', 5),
    ('VAL106P', 5),
    ('SER103P', 5),
    ('SER43P', 5),
    ('CYS17P', 5),
    ('CYS62P', 5),
    ('THR78P', 5),
    ('ALA151P', 5),
    ('ASP56P', 5),
    ('GLU139P', 5),
    ('TYR142P', 6),
    ('GLU114P', 6),
    ('TYR87P', 6),
    ('PRO69P', 6),
    ('LEU153P', 6),
    ('ASP81P', 7),
    ('CYS90P', 7),
    ('SER7P', 7),
    ('SER118P', 7),
    ('TYR57P', 7),
    ('LYS112P', 7),
    ('HSE66P', 7),
    ('GLN33P', 7),
    ('THR140P', 8),
    ('GLN144P', 8),
    ('ASP98P', 8),
    ('LEU116P', 8),
    ('LYS64P', 8),
    ('GLY133P', 8),
    ('MET70P', 8),
    ('GLY52P', 8),
    ('ASP32P', 9),
    ('ILE113P', 9),
    ('LYS110P', 9),
    ('CYS109P', 9),
    ('LYS155P', 9),
    ('GLU80P', 9),
    ('ASN104P', 9),
    ('GLU23P', 10),
    ('MET91P', 10),
    ('ASN100P', 10),
    ('VAL41P', 10),
    ('ASP135P', 10),
    ('GLU93P', 10),
    ('THR84P', 10),
    ('GLN60P', 10),
    ('TRP39P', 10),
    ('GLU128P', 10),
    ('LYS28P', 10),
    ('SER61P', 10),
    ('GLU154P', 10),
    ('ASP86P', 10),
    ('ASP120P', 10),
    ('LYS107P', 11),
    ('SER136P', 11),
    ('THR108P', 11),
    ('ALA156P', 11),
    ('LYS6P', 11),
    ('GLN122P', 11),
    ('SER47P', 12),
    ('GLN143P', 12),
    ('ARG97P', 12),
    ('ASN38P', 12),
    ('HSE157P', 12),
    ('ASP129P', 12),
    ('ASN53P', 12),
    ('THR5P', 12),
    ('GLY48P', 12),
    ('ASN34P', 13),
    ('HSE72P', 13),
    ('ASP42P', 13),
    ('GLU50P', 13),
    ('LYS123P', 13),
    ('ASN134P', 13),
    ('ARG101P', 13),
    ('SER94P', 14),
    ('ILE126P', 14),
    ('GLU37P', 14),
    ('GLN76P', 14),
    ('SER71P', 14),
    ('LYS79P', 14),
    ('TYR131P', 14),
    ('TYR132P', 14),
    ('GLN124P', 14),
    ('ASP137P', 15),
    ('GLN105P', 15),
    ('THR46P', 15),
    ('GLY117P', 15),
    ('ASN95P', 15),
    ('LEU125P', 15),
    ('ARG75P', 15),
    ('ARG18P', 15),
    ('ARG65P', 15),
    ('ARG40P', 16),
    ('ARG147P', 16),
    ('ARG58P', 16),
    ('ARG27P', 16),
    ('ASP92P', 16),
    ('TYR49P', 16),
    ('LYS102P', 16),
    ('ARG150P', 16),
    ('SER36P', 16)]

*clip* option can be used to include different number of results on the histogram.


.. ipython:: python
   :verbatim:    

    calcBridgingResiduesHistogram(waterBridges, clip=25)

.. figure:: images/traj_res_hist2.png
   :scale: 60 %

If we are interested in one particular residue, we can also use
:func:`.calcWaterBridgesDistribution` to find their partners in water bridges. 
Below we can see results for arginine 147 or aspartic acid 92 from chain P.


.. ipython:: python
   :verbatim:

    calcWaterBridgesDistribution(waterBridges, 'ARG147P')

.. parsed-literal::

   [('GLN122P', 8),
    ('ARG150P', 7),
    ('GLN143P', 6),
    ('LYS123P', 6),
    ('GLN124P', 5),
    ('ASP120P', 5),
    ('GLN144P', 3),
    ('THR140P', 2)]

.. ipython:: python
   :verbatim:

    calcWaterBridgesDistribution(waterBridges, 'ASP92P')

.. parsed-literal::

   [('ARG18P', 11),
    ('ASN95P', 10),
    ('SER94P', 5),
    ('MET91P', 5),
    ('ASP129P', 4),
    ('LEU13P', 3),
    ('CYS90P', 1)]

Once we select a pair of residues which are supported by interactions with water 
molecules we can use :func:`.calcWaterBridgesDistribution` to obtain histograms 
with results such as distances between them *(metric='distance')*, the number of 
water molecules which were involved *(metric='waters')*, and information about 
residue part which was involved in water bridges, i.e. backbone or side chain 
*(metric='location')*. 

.. ipython:: python
   :verbatim:

   calcWaterBridgesDistribution(waterBridges,  'ASP92P', 'ARG18P', trajectory=trajectory, metric='distance')

.. figure:: images/traj_distribution.png
   :scale: 50 %

.. parsed-literal::

   [5.3736005,
    5.3736005,
    5.167575,
    2.681302,
    5.371548,
    2.6318514,
    3.0394073,
    4.0884595,
    5.4406505,
    3.4112484,
    2.805657,
    5.4176636,
    3.5104342,
    5.991175,
    5.470093,
    3.4345005,
    3.6427624]

.. ipython:: python
   :verbatim:

   calcWaterBridgesDistribution(waterBridges, 'ARG147P', 'GLN122P', metric='waters') 

.. figure:: images/traj_distribution2.png
   :scale: 60 %

.. parsed-literal::

   [2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2]

.. ipython:: python
   :verbatim:

   calcWaterBridgesDistribution(waterBridges, 'ARG147P', 'GLN122P', trajectory=trajectory, metric='location')

.. parsed-literal::

   {'ARG147P': {'backbone': 7, 'side': 86},
   'GLN122P': {'backbone': 21, 'side': 25}}


Save results as PDB file
-------------------------------------------------------------------------------

The :meth:`.interactionsTrajectory.calcProteinInteractionsTrajectory` method 
computes interactions using default parameters for interactions. However, it 
can be changed according to our needs. To do that, we need to recalculate the
selected types of interactions. 

The results can be storage as PDB file using :func:`.savePDBWaterBridges` 
(single PDB file, single frame) or using :func:`.savePDBWaterBridgesTrajectory`
to save all the results (large number of frames saved each independently).

5kqm_all_sci_multi_0.pdb  5kqm_all_sci_multi_4.pdb  
5kqm_all_sci_multi_1.pdb  5kqm_all_sci_multi_5.pdb  
5kqm_all_sci_multi_2.pdb  5kqm_all_sci_multi_6.pdb  
5kqm_all_sci_multi_3.pdb  5kqm_all_sci_multi_7.pdb  

5kqm_all_sci_multi_8.pdb   5kqm_all_sci_multi_12.pdb
5kqm_all_sci_multi_9.pdb   5kqm_all_sci_multi_13.pdb
5kqm_all_sci_multi_10.pdb  5kqm_all_sci_multi_14.pdb
5kqm_all_sci_multi_11.pdb  5kqm_all_sci_multi_15.pdb


Those results can be displayed in any program for visualization. The results 
for protein structure will be storage in beta column (average values of 
contributions of each residue in water bridging) and occupancy column 
(results for particular frame). Water molecules will be included in each frame.


.. ipython:: python
   :verbatim:

   savePDBWaterBridges(waterBridges[0], coords_traj, PDBtraj_file[:-4]+'_frame0.pdb')

   savePDBWaterBridgesTrajectory(waterBridges, coords_traj, filename=PDBtraj_file[:-4]+'_multi.pdb', trajectory=trajectory)


Results saved in PDB file can be displayed as follows:


.. figure:: images/Fig2.tga
   :scale: 50 %


