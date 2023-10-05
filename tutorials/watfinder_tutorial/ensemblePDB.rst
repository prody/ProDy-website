.. _watfinder_tutorial:

Water bridges detection in Ensemble PDB
===============================================================================

This time we will use multi-model PDB which contain 50 frames from MD 
simulations form PE-binding protein 1 (PDB code: *1beh*). Simulation were 
performed using NAMD and saved as multi-model PDB using VMD. Remember to align 
the protein structure before analyzing it. Otherwise when all structures will 
be uploaded to the visualization program they will be spread out in space.


Parse structure
-------------------------------------------------------------------------------

.. ipython:: python
   :verbatim:

   ens = 'pebp1_50frames.pdb'
   coords_ens = parsePDB(ens)
   bridgeFrames_ens = calcWaterBridgesTrajectory(coords_ens, coords_ens)

.. parsed-literal::

   @> 20195 atoms and 51 coordinate set(s) were parsed in 1.88s.
   @> Frame: 0
   @> 161 water bridges detected.
   @> Frame: 1
   @> 127 water bridges detected.
   @> Frame: 2
   @> 168 water bridges detected.
   @> Frame: 3
   @> 132 water bridges detected.
   @> Frame: 4
   @> 142 water bridges detected.
   @> Frame: 5
   @> 166 water bridges detected.
   @> Frame: 6
   @> 150 water bridges detected.
   @> Frame: 7
   @> 159 water bridges detected.
   @> Frame: 8
   @> 147 water bridges detected.
   @> Frame: 9
   @> 136 water bridges detected.
   @> Frame: 10
   @> 127 water bridges detected.
   @> Frame: 11
   @> 135 water bridges detected.
   @> Frame: 12
   @> 162 water bridges detected.
   @> Frame: 13
   @> 135 water bridges detected.
   @> Frame: 14
   @> 178 water bridges detected.
   @> Frame: 15
   @> 128 water bridges detected.
   @> Frame: 16
   @> 145 water bridges detected.
   @> Frame: 17
   @> 157 water bridges detected.
   @> Frame: 18
   @> 111 water bridges detected.
   @> Frame: 19
   @> 122 water bridges detected.
   @> Frame: 20
   @> 171 water bridges detected.
   @> Frame: 21
   @> 144 water bridges detected.
   @> Frame: 22
   @> 152 water bridges detected.
   @> Frame: 23
   @> 155 water bridges detected.
   @> Frame: 24
   @> 164 water bridges detected.
   @> Frame: 25
   @> 138 water bridges detected.
   @> Frame: 26
   @> 156 water bridges detected.
   @> Frame: 27
   @> 118 water bridges detected.
   @> Frame: 28
   @> 165 water bridges detected.
   @> Frame: 29
   @> 108 water bridges detected.
   @> Frame: 30
   @> 136 water bridges detected.
   @> Frame: 31
   @> 117 water bridges detected.
   @> Frame: 32
   @> 146 water bridges detected.
   @> Frame: 33
   @> 139 water bridges detected.
   @> Frame: 34
   @> 142 water bridges detected.
   @> Frame: 35
   @> 125 water bridges detected.
   @> Frame: 36
   @> 123 water bridges detected.
   @> Frame: 37
   @> 147 water bridges detected.
   @> Frame: 38
   @> 146 water bridges detected.
   @> Frame: 39
   @> 143 water bridges detected.
   @> Frame: 40
   @> 162 water bridges detected.
   @> Frame: 41
   @> 151 water bridges detected.
   @> Frame: 42
   @> 146 water bridges detected.
   @> Frame: 43
   @> 131 water bridges detected.
   @> Frame: 44
   @> 121 water bridges detected.
   @> Frame: 45
   @> 126 water bridges detected.
   @> Frame: 46
   @> 145 water bridges detected.
   @> Frame: 47
   @> 112 water bridges detected.
   @> Frame: 48
   @> 151 water bridges detected.
   @> Frame: 49
   @> 131 water bridges detected.
   @> Frame: 50
   @> 175 water bridges detected.

Analysis of the results is similar to the one presented in trajectory analysis. 
Below examples showing which residues are the most frequently involved in water 
bridges formation (:func:`.calcBridgingResiduesHistogram`), details of that 
interactions (:func:`.calcWaterBridgesStatistics`), and results saved as PDB 
structure for further visualization (:func:`.savePDBWaterBridgesTrajectory`). 
Other functions can be seen in the analysis of trajectory.

.. ipython:: python
   :verbatim:

   calcBridgingResiduesHistogram(bridgeFrames_ens)


.. figure:: images/ensemble_hist.png
   :scale: 60 %

.. parsed-literal::

   [('VAL34P', 1),        
    ('VAL177P', 1),
    ('PRO43P', 1),
    ('LEU41P', 2),
    ('MET92P', 2),
    ('VAL164P', 3),
    ('LEU14P', 3),
    ('TYR169P', 3),
    ('PHE154P', 4),
    ('THR167P', 4),
    ('PRO178P', 4),
    ('LEU25P', 5),
    ('SER104P', 6),
    ('PRO136P', 6),
    ('ILE53P', 7),
    ('GLN170P', 7),
    ('TYR29P', 8),
    ('TYR106P', 8),
    ('PRO50P', 8),
    ('GLY38P', 8),
    ('PRO71P', 9),
    ('LEU138P', 9),
    ('PRO24P', 10),
    ('ASP18P', 10),
    ('LEU103P', 10),
    ('GLN127P', 10),
    ('TYR181P', 11),
    ('ILE137P', 11),
    ('VAL107P', 12),
    ('LEU131P', 12),
    ('TYR158P', 13),
    ('PRO112P', 13),
    ('PRO11P', 13),
    ('PRO74P', 15),
    ('PRO163P', 16),
    ('LEU160P', 16),
    ('HSE85P', 16),
    ('LEU58P', 17),
    ('LEU184P', 17),
    ('PRO130P', 17),
    ('PRO21P', 18),
    ('VAL17P', 19),
    ('LEU63P', 19),
    ('ASP35P', 19),
    ('GLU36P', 20),
    ('TYR81P', 20),
    ('PRO79P', 20),
    ('VAL40P', 20),
    ('HSE86P', 21),
    ('TRP84P', 21),
    ('ASP175P', 23),
    ('GLU159P', 23),
    ('GLY116P', 23),
    ('LEU5P', 23),
    ('ASP174P', 25),
    ('GLU172P', 25),
    ('GLU182P', 26),
    ('THR51P', 26),
    ('GLU19P', 26),
    ('CYS133P', 27),
    ('GLY57P', 27),
    ('SER153P', 27),
    ('TRP55P', 27),
    ('LYS156P', 28),
    ('ILE97P', 28),
    ('ASP134P', 28),
    ('LYS148P', 29),
    ('GLY100P', 29),
    ('ALA73P', 29),
    ('VAL102P', 29),
    ('LYS179P', 30),
    ('ASP105P', 30),
    ('SER109P', 30),
    ('PHE149P', 31),
    ('LYS132P', 31),
    ('GLU135P', 31),
    ('GLY108P', 31),
    ('SER185P', 32),
    ('LEU117P', 32),
    ('HSE23P', 32),
    ('ALA162P', 32),
    ('GLY143P', 32),
    ('LYS7P', 32),
    ('GLY186P', 34),
    ('LEU12P', 34),
    ('ASP70P', 34),
    ('LYS113P', 35),
    ('LYS157P', 35),
    ('GLN45P', 35),
    ('THR42P', 35),
    ('SER99P', 35),
    ('ASN91P', 36),
    ('ASP4P', 37),
    ('GLY10P', 38),
    ('HSE26P', 39),
    ('TYR64P', 39),
    ('ASP72P', 39),
    ('ASP78P', 39),
    ('LYS77P', 40),
    ('GLU16P', 40),
    ('ALA152P', 40),
    ('ASP144P', 40),
    ('SER52P', 40),
    ('GLY61P', 40),
    ('SER60P', 40),
    ('GLN15P', 40),
    ('GLU83P', 41),
    ('VAL3P', 41),
    ('GLN22P', 41),
    ('ALA32P', 41),
    ('SER9P', 41),
    ('ASP59P', 42),
    ('LEU37P', 42),
    ('THR101P', 42),
    ('SER75P', 43),
    ('ASN95P', 43),
    ('THR69P', 43),
    ('THR28P', 43),
    ('THR44P', 43),
    ('GLN20P', 43),
    ('TRP8P', 43),
    ('TYR176P', 43),
    ('ASP56P', 44),
    ('GLY110P', 44),
    ('SER54P', 44),
    ('SER142P', 44),
    ('TRP173P', 44),
    ('ASP96P', 44),
    ('GLY114P', 44),
    ('SER6P', 44),
    ('HSE118P', 45),
    ('SER13P', 45),
    ('ASN140P', 45),
    ('SER139P', 45),
    ('GLY31P', 45),
    ('GLN183P', 45),
    ('ARG155P', 45),
    ('LYS47P', 45),
    ('ALA33P', 45),
    ('SER98P', 46),
    ('ALA30P', 47),
    ('HSE145P', 47),
    ('LYS39P', 47),
    ('LYS62P', 47),
    ('THR115P', 48),
    ('VAL151P', 48),
    ('ASP128P', 48),
    ('LYS93P', 48),
    ('TYR120P', 49),
    ('GLY147P', 49),
    ('GLU126P', 49),
    ('GLY94P', 49),
    ('ARG49P', 50),
    ('ASN48P', 50),
    ('ARG141P', 50),
    ('LYS150P', 50),
    ('ARG119P', 51),
    ('LYS80P', 51),
    ('ARG76P', 51),
    ('ARG161P', 51),
    ('ARG129P', 51),
    ('ARG82P', 51),
    ('ARG146P', 51)]

.. ipython:: python
   :verbatim:

   analysisAtomic_ens = calcWaterBridgesStatistics(bridgeFrames_ens, coords_ens)

   for item in analysisAtomic_ens.items():
      print(item)

.. parsed-literal::

   @> RES1           RES2           PERC      DIST_AVG  DIST_STD
   @> VAL3P          HSE26P         19.608    5.581     0.696
   @> ASP4P          SER6P          13.725    3.817     0.560
   @> SER6P          LYS7P          43.137    4.394     1.114
   @> LYS7P          GLU36P         1.961     6.088     0.000
   @> LYS7P          LEU37P         7.843     6.353     0.433
   @> GLY10P         SER13P         43.137    4.759     0.612
   @> GLY10P         ARG76P         11.765    5.309     0.586
   @> LEU12P         SER13P         45.098    2.767     0.080
   @> SER13P         GLU16P         25.490    4.449     1.133
   @> GLN15P         ASP18P         7.843     3.732     0.174
   @> GLU16P         ARG82P         45.098    4.550     1.086
   @> GLU16P         VAL17P         17.647    3.438     0.952
   @> GLU16P         LYS150P        21.569    5.056     0.929
   @> GLU16P         GLU83P         9.804     5.476     1.138
   @> GLU16P         ALA152P        7.843     7.307     0.450
   @> VAL17P         GLU83P         1.961     7.262     0.000
   @> VAL17P         LYS150P        13.725    6.303     0.572
   @> GLN22P         GLU126P        33.333    6.458     1.216
   @> GLN22P         HSE23P         37.255    4.738     0.669
   @> HSE23P         GLU126P        7.843     7.911     0.239
   @> PRO24P         ASP56P         17.647    5.592     0.910
   @> THR28P         SER52P         43.137    3.970     0.677
   @> THR28P         ILE53P         5.882     5.849     0.027
   @> TYR29P         THR51P         7.843     3.583     0.286
   @> ALA30P         ARG49P         17.647    5.206     0.304
   @> GLY31P         ALA32P         78.431    2.875     0.278
   @> ALA33P         LYS39P         82.353    4.668     0.469
   @> ASP35P         LYS39P         23.529    3.191     0.719
   @> THR42P         THR44P         52.941    4.123     0.501
   @> GLN45P         LYS47P         3.922     3.255     0.182
   @> GLN45P         ASN48P         7.843     4.808     0.469
   @> LYS47P         ASN48P         64.706    3.898     1.056
   @> LYS47P         GLN183P        37.255    4.918     1.111
   @> ASN48P         ARG49P         58.824    5.056     1.055
   @> ASN48P         ASP105P        9.804     8.198     0.487
   @> ASN48P         TYR106P        9.804     7.210     0.618
   @> ASN48P         ASN140P        25.490    6.728     0.817
   @> ARG49P         ASP105P        45.098    4.940     1.916
   @> ARG49P         TYR106P        7.843     6.351     0.321
   @> ARG49P         ASN140P        41.176    6.343     0.865
   @> SER52P         ILE53P         7.843     4.284     0.315
   @> SER54P         ASP56P         5.882     6.349     0.389
   @> ASP56P         GLY57P         33.333    3.059     0.594
   @> LEU58P         ASP59P         13.725    3.697     0.345
   ..
   ..

.. ipython:: python
   :verbatim:

   savePDBWaterBridgesTrajectory(bridgeFrames_ens, coords_ens, ens[:-4]+'_ens.pdb')

.. parsed-literal::

   @> All 51 coordinate sets are copied to pebp1_50frames Selection 'protein' + pebp1_50frames Selection 'same residue as...6074 4190 14360'.
   @> All 51 coordinate sets are copied to pebp1_50frames Selection 'protein' + pebp1_50frames Selection 'same residue as...9718 17936 7184'.
   @> All 51 coordinate sets are copied to pebp1_50frames Selection 'protein' + pebp1_50frames Selection 'same residue as...947 10043 11756'.
   @> All 51 coordinate sets are copied to pebp1_50frames Selection 'protein' + pebp1_50frames Selection 'same residue as...0099 12848 4175'.
   @> All 51 coordinate sets are copied to pebp1_50frames Selection 'protein' + pebp1_50frames Selection 'same residue as...6031 8645 18008'.
   ..
   ..


Detecting water centers
-------------------------------------------------------------------------------

Previous function generated multiple PDB in which we would found protein and 
water molecules for each frame that form water bridges with protein structure. 
Now we can use another function :func:`.findClusterCenters` which will extract 
water centers (they refer to the oxygens from water molecules that are forming 
clusters). We need to provide a file pattern as show below. Now all the files 
with prefix *'pebp1_50frames_ens_'* will be analyzed.


.. ipython:: python
   :verbatim:

   findClusterCenters('pebp1_50frames_ens_*.pdb')

.. parsed-literal::

   @> 3269 atoms and 1 coordinate set(s) were parsed in 0.11s.
   @> 3161 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 3173 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3173 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3218 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3251 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3215 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3230 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3230 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3224 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3158 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3176 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3218 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3284 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3227 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3251 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3233 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3254 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3197 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3197 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3209 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3236 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3194 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3227 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3179 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3215 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3284 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3188 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3176 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3203 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3227 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3269 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3191 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3245 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3225 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3261 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3221 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3233 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3209 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3167 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3221 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3275 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3167 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3218 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3191 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3221 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3251 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3209 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> 3278 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3239 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 3228 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> Results are saved in clusters_pebp1_50frames_ens_.pdb.

Function generated one PDB file with water centers. We used default values, 
such as distC (distance to other molecule) and numC (min number of molecules 
in a cluster), but those values could be changed if the molecules are more 
widely distributed or we would like to have more numerous clusters.
Moreover, this function can be applied on different type of molecules by using 
*selection* paramater. We can provide the whole molecule and by default 
the center of mass will be used as a reference.


Saved PDB files using :func:`.savePDBWaterBridgesTrajectory` in the previous
step can be upload to VMD or other program for visualization:

.. figure:: images/Fig3.tga
   :scale: 50 %


After uploading a new PDB file with water centers we can see the results as
follows:

.. figure:: images/Fig4.tga
   :scale: 50 %
