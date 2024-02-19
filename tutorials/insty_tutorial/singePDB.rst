.. _insty_tutorial:

Interactions/Stability Evaluation (InSty)
===============================================================================

This example shows how to perform Interactions/Stability Evaluation
(**InSty**) analysis for a small protein (<200 residues) called tyrosine
phosphatase LMW-PTP (**5KQM**) and visualize the results using Matplotlib_
library and VMD_ program. 
In the tutorial, we will use already preapared structure for
simulation (with hydrogens added). The same structure will be later
analyzed with the trajectory file to show how the analysis of interactions 
in the course of simulation can change. 

The tutorial will also include an example of a PDB structure directly
downloaded from Protein Data Bank (PDB) which requires adding the missing hydrogen
atoms to the protein and ligand structure. The analysis will be performed for
protein-ligand interactions.


Analysis of interactions for a single PDB structure
-------------------------------------------------------------------------------

We start by parsing PDB file with LMW-PTP **5kqm_all_sci.pdb** which is avalable
as the tutorial files. PDB file contains protein structures with water and 
counter ions prepared using VMD_ program.

Before that import everything from the ProDy packages.

.. ipython:: python

   from prody import *
   from pylab import *
   import matplotlib
   ion()   # turn interactive mode on


.. ipython:: python
   :verbatim:

   PDBfile = '5kqm_all_sci.pdb'
   coords = parsePDB(PDBfile)
   coords

.. parsed-literal::

   @> 19321 atoms and 1 coordinate set(s) were parsed in 0.23s.

For the analysis we will use only protein coordinates (*atoms*):

.. ipython:: python
   :verbatim:

   atoms = coords.select('protein')
   atoms

.. parsed-literal::

   @> 19321 atoms and 1 coordinate set(s) were parsed in 0.21s.


Compute all types of interactions
-------------------------------------------------------------------------------

In the next step, we instantiate an :class:`.Interactions` instance:

.. ipython:: python
   :verbatim:

   interactions = Interactions()

Now we can compute all available types of interactions (seven types: hydrogen
bonds, salt bridges, repulsive ionic bonding, Pi-cation, Pi-stacking,
hydrphobic interactions, and disulfide bonds) for protein structure by passing
selected atoms (*atoms*) to :meth:`.Interactions.calcProteinInteractions` method:

.. ipython:: python
   :verbatim:

   all_interactions = interactions.calcProteinInteractions(atoms)

.. parsed-literal::

   @> Calculating interations.
   @> Calculating hydrogen bonds.
   @>      DONOR (res chid atom)   <--->       ACCEPTOR (res chid atom)    Distance  Angle
   @>     ARG101    P      NH1_1516  <--->      ASP98    P      OD1_1463     2.0    33.1
   @>      HSE72    P      NE2_1042  <--->      ASN15    P       OD1_165     2.6    34.8
   @>     GLN143    P      NE2_2192  <--->     GLU139    P      OE2_2126     2.7     9.2
   @>      HSE66    P       NE2_957  <--->     GLU139    P      OE1_2125     2.7     6.4
   @>      ARG40    P         N_561  <--->       LYS6    P          O_37     2.7    17.1
   @>      ARG58    P         N_813  <--->      ASP56    P       OD1_788     2.7    30.0
   @>      ALA45    P         N_634  <--->      ARG75    P        O_1097     2.8    35.1
   @>      ASN53    P       ND2_747  <--->      GLU50    P       OE1_708     2.8    18.2
   @>      ALA74    P        N_1064  <--->      ASN53    P         O_751     2.8    21.3
   @>      ASP56    P         N_780  <--->      ILE16    P         O_189     2.8    27.0
   @>     LYS110    P       NZ_1667  <--->      THR84    P        O_1240     2.8    38.2
   @>     LEU116    P        N_1758  <--->      CYS90    P        O_1342     2.8    15.0
   @>     SER103    P        N_1546  <--->      LEU99    P        O_1485     2.8    29.1
   @>     ASN134    P        N_2045  <--->     ASP137    P      OD2_2091     2.8    22.6
   @>     PHE152    P        N_2321  <--->     CYS148    P        O_2275     2.8     8.3
   @>      ASN95    P        N_1398  <--->      ASP92    P      OD1_1368     2.8    12.6
   @>       LYS6    P          N_16  <--->      ASN38    P         O_536     2.8    25.0
   @>      ILE77    P        N_1115  <--->      ALA45    P         O_643     2.8    12.2
   @>      ARG58    P       NH2_832  <--->      ASP56    P       OD2_789     2.8    27.7
   @>      LEU99    P        N_1467  <--->      ASN95    P        O_1411     2.8    15.5
   @>     CYS149    P        N_2276  <--->     CYS145    P        O_2224     2.8     9.6
   @>      GLY52    P         N_731  <--->      ALA74    P        O_1073     2.8     6.6
   @>      ASP32    P         N_435  <--->      LYS28    P         O_385     2.8     8.8
   @>      ILE88    P        N_1294  <--->     LYS112    P        O_1704     2.8    17.7
   @>     GLN143    P        N_2180  <--->     GLU139    P        O_2128     2.8    21.7
   @>      ARG27    P         N_340  <--->      GLU23    P         O_293     2.8    15.4
   @>     TYR142    P        N_2159  <--->     PHE138    P        O_2113     2.9    14.2
   @>     GLY133    P        N_2038  <--->     PRO130    P        O_1995     2.9    25.4
   @>      PHE26    P         N_320  <--->      ALA22    P         O_278     2.9     4.9
   @>      ASN15    P       ND2_166  <--->      SER19    P        OG_232     2.9    32.1
   @>      ARG75    P      NH1_1090  <--->      ASP81    P      OD2_1194     2.9    19.7
   @>      ARG75    P      NH2_1093  <--->      ASP42    P       OD2_610     2.9    23.5
   @>      ARG97    P        N_1431  <--->      GLU93    P        O_1386     2.9    22.2
   @>      ARG65    P       NH2_941  <--->     GLU139    P      OE1_2125     2.9    32.3
   @>      VAL25    P         N_304  <--->      ILE21    P         O_268     2.9     8.2
   @>     LEU153    P        N_2341  <--->     CYS149    P        O_2286     2.9    12.5
   @>       SER7    P          N_38  <--->      ASP86    P      OD2_1270     2.9    39.9
   @>      ASP86    P        N_1261  <--->       SER7    P         OG_45     2.9    34.7
   @>      ARG58    P       NH2_832  <--->     TYR131    P        O_2016     2.9    33.1
   @>      THR46    P         N_644  <--->      CYS12    P         O_130     2.9    36.1
   @>     GLN144    P        N_2197  <--->     THR140    P        O_2142     2.9    23.3
   @>      THR78    P        N_1134  <--->      ASP81    P      OD2_1194     2.9    12.4
   @>      LEU89    P        N_1313  <--->       LEU9    P          O_83     2.9    29.5
   @>      THR31    P         N_421  <--->      ARG27    P         O_363     2.9    24.1
   @>      CYS90    P        N_1332  <--->     GLU114    P        O_1738     2.9    24.6
   @>     CYS148    P        N_2265  <--->     GLN144    P        O_2213     2.9     9.3
   @>      GLU23    P         N_279  <--->      SER19    P         O_235     2.9    15.4
   @>      ILE68    P         N_970  <--->      MET63    P         O_899     2.9    13.0
   @>      PHE10    P          N_84  <--->      ASP42    P         O_612     2.9    22.8
   @>     LYS112    P        N_1683  <--->      ASP86    P        O_1272     2.9    10.1
   @>      SER61    P         N_861  <--->      TYR57    P         O_812     2.9    35.1
   @>     CYS145    P        N_2214  <--->     VAL141    P        O_2158     2.9    15.9
   @>      ARG27    P       NH2_359  <--->      GLU23    P       OE2_291     2.9    31.5
   @>      LYS64    P         N_900  <--->      GLN60    P         O_860     2.9    22.9
   @>       LEU9    P          N_65  <--->      TYR87    P        O_1293     2.9    16.4
   @>      ASN38    P         N_523  <--->      ILE35    P         O_496     2.9    29.1
   @>      VAL11    P         N_104  <--->      LEU89    P        O_1331     2.9    29.7
   @>     ASN100    P        N_1486  <--->      LEU96    P        O_1430     2.9    10.3
   @>     GLN124    P        N_1881  <--->     ASP120    P      OD2_1825     2.9    27.5
   @>     LYS102    P        N_1524  <--->      ASP98    P        O_1466     2.9     9.3
   @>      GLN76    P      NE2_1110  <--->      THR46    P         O_657     2.9    31.4
   @>      ARG40    P       NH1_577  <--->      THR84    P      OG1_1233     2.9     8.4
   @>      ALA44    P         N_624  <--->      PHE10    P         O_103     2.9    33.2
   @>     GLU154    P        N_2360  <--->     ARG150    P        O_2310     3.0    22.6
   @>       VAL8    P          N_49  <--->      ARG40    P         O_584     3.0    25.0
   @>      MET63    P         N_883  <--->      GLY59    P         O_843     3.0    18.3
   @>      GLN60    P         N_844  <--->      ASP56    P         O_791     3.0    35.5
   @>      ILE35    P         N_478  <--->      VAL30    P         O_420     3.0    23.5
   @>     VAL146    P        N_2225  <--->     TYR142    P        O_2179     3.0    31.5
   @>      ARG58    P       NH1_829  <--->     TYR131    P        O_2016     3.0    38.1
   @>      ASN53    P         N_738  <--->      GLU50    P         O_711     3.0    28.6
   @>     ARG101    P        N_1500  <--->      ARG97    P        O_1454     3.0    32.3
   @>      ARG18    P       NH1_217  <--->     ILE127    P        O_1954     3.0    26.0
   @>      ARG75    P        N_1074  <--->      ASN15    P       OD1_165     3.0    25.1
   @>     GLN144    P      NE2_2209  <--->     ILE126    P        O_1935     3.0    18.3
   @>      ASN34    P         N_464  <--->      THR31    P         O_434     3.0    18.2
   @>      ASN15    P       ND2_166  <--->      SER43    P        OG_620     3.0    25.7
   @>      ARG58    P        NE_826  <--->      ASP56    P       OD1_788     3.0    22.2
   @>      ARG27    P       NH1_356  <--->      GLU23    P       OE2_291     3.0    36.9
   @>     ILE127    P        N_1936  <--->      MET91    P        O_1359     3.0    17.6
   @>     TYR119    P       OH_1808  <--->     HSE157    P        N_2407     3.0    28.1
   @>     HSE157    P        N_2407  <--->     TYR119    P       OH_1808     3.0    19.2
   @>     GLU139    P        N_2114  <--->     ASP135    P        O_2070     3.0    27.9
   @>      LEU29    P         N_386  <--->      VAL25    P         O_319     3.0    19.1
   @>      SER47    P         N_658  <--->      LEU13    P         O_149     3.0    28.8
   @>      VAL30    P         N_405  <--->      PHE26    P         O_339     3.0    17.7
   @>     GLN105    P        N_1571  <--->     LYS102    P        O_1545     3.0    19.7
   @>     SER118    P        N_1784  <--->     LEU115    P        O_1757     3.1    21.4
   @>     LYS155    P        N_2375  <--->     ALA151    P        O_2320     3.1    21.3
   @>     GLU114    P        N_1724  <--->      ILE88    P        O_1312     3.1    24.2
   @>     ASP120    P        N_1816  <--->     GLY117    P        O_1783     3.1    12.7
   @>      CYS62    P         N_872  <--->      ARG58    P         O_836     3.1    20.4
   @>      ARG18    P       NH1_217  <--->      ASP92    P      OD2_1369     3.1     4.3
   @>      ALA24    P         N_294  <--->      PRO20    P         O_249     3.1    29.9
   @>     ARG150    P        N_2287  <--->     VAL146    P        O_2240     3.1    12.7
   @>      LYS28    P         N_364  <--->      ALA24    P         O_303     3.1    20.0
   @>     VAL141    P        N_2143  <--->     ASP137    P        O_2093     3.1    18.5
   @>      ASP98    P        N_1455  <--->      SER94    P        O_1397     3.1    19.6
   @>      LEU96    P        N_1412  <--->      ASP92    P        O_1371     3.1    36.3
   @>      ALA22    P         N_269  <--->      ARG18    P         O_224     3.1    21.9
   @>     ALA151    P        N_2311  <--->     ARG147    P        O_2264     3.1    15.6
   @>      GLY67    P         N_963  <--->      LYS64    P         O_921     3.1    22.8
   @>      ASP42    P         N_601  <--->       VAL8    P          O_64     3.1    35.6
   @>      ARG65    P         N_922  <--->      SER61    P         O_871     3.1    23.4
   @>      TRP39    P         N_537  <--->      SER36    P         O_507     3.1    15.2
   @>     LYS123    P        N_1859  <--->     ASP120    P        O_1827     3.1    18.7
   @>      MET91    P        N_1343  <--->      ASN95    P      OD1_1406     3.2    39.0
   @>     THR140    P        N_2129  <--->     SER136    P        O_2081     3.2    30.3
   @>      PHE85    P        N_1241  <--->      ASP81    P        O_1196     3.2    20.2
   @>      ASN15    P         N_157  <--->      CYS12    P        SG_127     3.2    37.5
   @>     ALA111    P        N_1673  <--->      PHE82    P        O_1216     3.2    20.6
   @>     ARG147    P        N_2241  <--->     GLN143    P        O_2196     3.2    12.1
   @>      ARG75    P      NH2_1093  <--->      ASP81    P      OD1_1193     3.2    29.3
   @>     LYS112    P       NZ_1699  <--->     HSE157    P      OT2_2424     3.3    28.7
   @>     ARG147    P      NH1_2257  <--->     GLN124    P      OE1_1892     3.3    29.9
   @>     PHE138    P        N_2094  <--->     ASN134    P        O_2058     3.3    31.0
   @>       SER7    P         OG_45  <--->      THR84    P        O_1240     3.3    35.5
   @>      CYS12    P         N_120  <--->      ALA44    P         O_633     3.3    36.1
   @>      SER19    P         N_225  <--->      CYS12    P        SG_127     3.3     8.0
   @>      PHE82    P        N_1197  <--->      LYS79    P        O_1169     3.4    37.7
   @>      ASP81    P        N_1185  <--->      THR78    P      OG1_1140     3.5    39.5
   @>     LYS102    P       NZ_1540  <--->      ASP98    P      OD2_1464     3.5    26.1
   @>     ARG147    P      NH2_2260  <--->     GLN124    P      OE1_1892     3.5    33.9
   @>     VAL106    P        N_1588  <--->     SER103    P        O_1556     3.5    34.2
   @> Number of detected hydrogen bonds: 124.
   @> Calculating salt bridges.
   @>      HSE66    P         NE2_957  <--->     GLU139    P   OE1_2125_2126     2.8
   @>      ASP81    P   OD1_1193_1194  <--->      ARG75    P   NH1_1090_1093     2.9
   @>      ASP32    P     OD1_443_444  <--->      LYS28    P          NZ_380     3.0
   @>     ARG101    P   NH1_1516_1519  <--->      ASP98    P   OD1_1463_1464     3.1
   @>      ARG27    P     NH1_356_359  <--->      GLU23    P     OE1_290_291     3.7
   @>     GLU139    P   OE1_2125_2126  <--->      ARG65    P     NH1_938_941     3.8
   @>     LYS102    P         NZ_1540  <--->      ASP98    P   OD1_1463_1464     3.9
   @>      ARG58    P     NH1_829_832  <--->      ASP56    P     OD1_788_789     3.9
   @>      ARG18    P     NH1_217_220  <--->      ASP92    P   OD1_1368_1369     4.1
   @>     GLU114    P   OE1_1735_1736  <--->     LYS112    P         NZ_1699     4.1
   @>     ASP120    P   OD1_1824_1825  <--->     ARG147    P   NH1_2257_2260     4.2
   @>     LYS110    P         NZ_1667  <--->      ASP86    P   OD1_1269_1270     4.2
   @>     GLU114    P   OE1_1735_1736  <--->     HSE157    P        NE2_2418     4.4
   @>      ARG18    P     NH1_217_220  <--->     ASP129    P   OD1_1978_1979     4.6
   @>      ARG75    P   NH1_1090_1093  <--->      ASP42    P     OD1_609_610     4.6
   @>      GLU23    P     OE1_290_291  <--->      HSE72    P        NE2_1042     5.0
   @> Number of detected salt bridges: 16.
   @> Calculating repulsive ionic bonding.
   @>     ARG101    P   NH1_1516_1519  <--->     LYS102    P         NZ_1540     4.3
   @> Number of detected Repulsive Ionic Bonding interactions: 1.
   @> Calculating Pi stacking interactions.
   @>      HSE66       P             953_954_955_957_959  <--->     TYR142       P   2166_2167_2169_2171_2174_2176     3.9   162.1
   @>     HSE157       P2414_2415_2416_2418_2420_2423_2424  <--->     TYR119       P   1802_1803_1805_1807_1810_1812     4.4     3.0
   @>      TRP39       P         549_550_551_553_555_557  <--->      PHE26       P         327_328_330_332_334_336     4.8    75.5
   @>     TYR131       P   2003_2004_2006_2008_2011_2013  <--->     TYR132       P   2024_2025_2027_2029_2032_2034     4.9    91.4
   @> Number of detected Pi stacking interactions: 4.
   @> Calculating cation-Pi interactions.
   @>      PHE85   P   1248_1249_1251_1253_1255_1257  <--->      ARG40   P                     NH1_577_580     3.7
   @>      HSE66   P             953_954_955_957_959  <--->      ARG65   P                     NH1_938_941     4.5
   @>     HSE157   P2414_2415_2416_2418_2420_2423_2424  <--->     LYS112   P                         NZ_1699     4.8
   @> Number of detected cation-pi interactions: 3.
   @> Hydrophobic Overlaping Areas are computed.
   @> Calculating hydrophobic interactions.
   @>      TYR87    P       OH_128614s  <--->     ALA156    P       CB_2401     3.0    22.0
   @>      MET63    P        CE_89414s  <--->      ALA24    P        CB_298     3.3     5.2
   @>      ILE68    P       CG2_97614s  <--->      MET63    P        CE_894     3.3    52.4
   @>     TYR142    P       CZ_217114s  <--->     VAL146    P      CG2_2235     3.5    49.7
   @>      PHE10    P        CD1_9214s  <--->      ALA22    P        CB_273     3.5    31.2
   @>       LYS6    P         CD_2614s  <--->      TRP39    P       CZ2_555     3.5    68.7
   @>      VAL30    P       CG1_41114s  <--->      PHE26    P       CE2_336     3.6    21.1
   @>     ALA111    P       CB_167714s  <--->      ILE88    P       CD_1307     3.6    21.2
   @>      VAL11    P       CG2_11414s  <--->      ILE88    P      CG2_1300     3.6     9.3
   @>      VAL41    P       CG2_59514s  <--->      PHE26    P       CD2_334     3.6    16.6
   @>     PHE152    P      CE1_233114s  <--->     ALA156    P       CB_2401     3.7    17.5
   @>     VAL106    P      CG2_159814s  <--->      LYS79    P       CG_1155     3.7    25.1
   @>      ILE77    P       CD_112814s  <--->      LEU99    P      CD2_1480     3.7    12.0
   @>      PHE82    P      CD1_120514s  <--->      ILE88    P       CD_1307     3.7    17.6
   @>     LEU116    P      CD2_177114s  <--->     ILE127    P       CD_1949     3.7    17.4
   @>       VAL8    P        CG1_5514s  <--->      PHE26    P       CE2_336     3.7    12.1
   @>      LEU96    P      CD1_142114s  <--->     ILE113    P      CG2_1711     3.7    17.0
   @>       LEU9    P        CD2_7814s  <--->      ILE77    P       CD_1128     3.7    15.4
   @>      LEU89    P      CD1_132214s  <--->       VAL8    P        CG2_59     3.8    15.9
   @>     ILE126    P       CD_193014s  <--->     LEU125    P      CD1_1907     3.8    54.2
   @>     VAL141    P      CG1_214914s  <--->     ILE127    P      CG2_1942     3.9    11.5
   @>      MET91    P       SD_135314s  <--->     ILE127    P       CD_1949     3.9    35.9
   @>      ALA44    P        CB_62814s  <--->       LEU9    P        CD1_74     3.9    15.1
   @>      VAL25    P       CG2_31414s  <--->     TYR142    P      CE1_2169     3.9    12.0
   @>      ILE21    P       CG2_25614s  <--->      MET63    P        SD_893     4.0    20.8
   @>     LEU153    P      CD1_235014s  <--->      TRP39    P       NE1_547     4.0     9.4
   @>      PHE85    P       CZ_125314s  <--->       LEU9    P        CD1_74     4.0    32.1
   @>      ILE35    P        CD_49114s  <--->      TRP39    P       NE1_547     4.0    26.0
   @>      LEU29    P       CD1_39514s  <--->      VAL25    P       CG1_310     4.1    19.7
   @>      ALA74    P       CB_106814s  <--->      ILE16    P       CG2_177     4.1     6.7
   @>      ARG75    P       CG_108114s  <--->      ALA44    P        CB_628     4.1    36.2
   @>      ARG18    P        CG_20814s  <--->     VAL141    P      CG1_2149     4.1    20.3
   @>     LYS102    P       CD_153414s  <--->      ILE77    P      CG2_1121     4.1    17.5
   @>     TYR119    P      CE1_180514s  <--->      LEU89    P      CD2_1326     4.1    11.6
   @>      ARG40    P        CG_56814s  <--->      PHE85    P      CE2_1257     4.3    60.9
   @>      LYS28    P        CG_37114s  <--->      ILE68    P        CD_983     4.3    21.8
   @>     PHE138    P      CD2_210814s  <--->      ILE21    P        CD_263     4.3     6.6
   @>     TYR131    P      CE1_200614s  <--->      ILE16    P        CD_184     4.3     8.9
   @>      ARG58    P        CG_82014s  <--->     PHE138    P      CE1_2104     4.5    59.4
   @> Number of detected hydrophobic interactions: 39.
   @> Calculating disulfide bonds.
   @> Number of detected disulfide bonds: 0.

All types of interactions will be displayed on the screen with all types of
information such as distance or angle (if applied).

Moreover, we will have access to the details of each interaction type
using the following methods: 

:meth:`.Interactions.getHydrogenBonds` - hydrogen bonds:

.. ipython:: python
   :verbatim:   

   interactions.getHydrogenBonds()

.. parsed-literal::

   [['ARG101', 'NH1_1516', 'P', 'ASP98', 'OD1_1463', 'P', 1.998, 33.1238],
    ['HSE72', 'NE2_1042', 'P', 'ASN15', 'OD1_165', 'P', 2.5997, 34.752],
    ['GLN143', 'NE2_2192', 'P', 'GLU139', 'OE2_2126', 'P', 2.7287, 9.1823],
    ['HSE66', 'NE2_957', 'P', 'GLU139', 'OE1_2125', 'P', 2.7314, 6.3592],
    ['ARG40', 'N_561', 'P', 'LYS6', 'O_37', 'P', 2.7479, 17.1499],
    ['ARG58', 'N_813', 'P', 'ASP56', 'OD1_788', 'P', 2.7499, 29.9737],
    ['ALA45', 'N_634', 'P', 'ARG75', 'O_1097', 'P', 2.7609, 35.0983],
    ['ASN53', 'ND2_747', 'P', 'GLU50', 'OE1_708', 'P', 2.7702, 18.2336],
    ['ALA74', 'N_1064', 'P', 'ASN53', 'O_751', 'P', 2.7782, 21.3375],
    ['ASP56', 'N_780', 'P', 'ILE16', 'O_189', 'P', 2.7793, 27.0481],
    ['LYS110', 'NZ_1667', 'P', 'THR84', 'O_1240', 'P', 2.7977, 38.2213],
    ['LEU116', 'N_1758', 'P', 'CYS90', 'O_1342', 'P', 2.8072, 15.0239],
    ['SER103', 'N_1546', 'P', 'LEU99', 'O_1485', 'P', 2.8075, 29.107],
    ['ASN134', 'N_2045', 'P', 'ASP137', 'OD2_2091', 'P', 2.8132, 22.562],
    ['PHE152', 'N_2321', 'P', 'CYS148', 'O_2275', 'P', 2.8141, 8.2562],
    ['ASN95', 'N_1398', 'P', 'ASP92', 'OD1_1368', 'P', 2.8148, 12.5701],
    ['LYS6', 'N_16', 'P', 'ASN38', 'O_536', 'P', 2.8178, 25.0305],
    ['ILE77', 'N_1115', 'P', 'ALA45', 'O_643', 'P', 2.8179, 12.1855],
    ['ARG58', 'NH2_832', 'P', 'ASP56', 'OD2_789', 'P', 2.8204, 27.6617],
    ['LEU99', 'N_1467', 'P', 'ASN95', 'O_1411', 'P', 2.8205, 15.4867],
    ['CYS149', 'N_2276', 'P', 'CYS145', 'O_2224', 'P', 2.8247, 9.5914],
    ['GLY52', 'N_731', 'P', 'ALA74', 'O_1073', 'P', 2.832, 6.6442],
    ['ASP32', 'N_435', 'P', 'LYS28', 'O_385', 'P', 2.8357, 8.8318],
    ['ILE88', 'N_1294', 'P', 'LYS112', 'O_1704', 'P', 2.8429, 17.7147],
    ['GLN143', 'N_2180', 'P', 'GLU139', 'O_2128', 'P', 2.8445, 21.6714],
    ['ARG27', 'N_340', 'P', 'GLU23', 'O_293', 'P', 2.8446, 15.4167],
    ['TYR142', 'N_2159', 'P', 'PHE138', 'O_2113', 'P', 2.8515, 14.2061],
    ['GLY133', 'N_2038', 'P', 'PRO130', 'O_1995', 'P', 2.854, 25.4301],
    ['PHE26', 'N_320', 'P', 'ALA22', 'O_278', 'P', 2.8541, 4.8732],
    ['ASN15', 'ND2_166', 'P', 'SER19', 'OG_232', 'P', 2.8592, 32.1244],
    ['ARG75', 'NH1_1090', 'P', 'ASP81', 'OD2_1194', 'P', 2.8632, 19.6664],
    ['ARG75', 'NH2_1093', 'P', 'ASP42', 'OD2_610', 'P', 2.8649, 23.5083],
    ['ARG97', 'N_1431', 'P', 'GLU93', 'O_1386', 'P', 2.8654, 22.24],
    ['ARG65', 'NH2_941', 'P', 'GLU139', 'OE1_2125', 'P', 2.8655, 32.3239],
    ['VAL25', 'N_304', 'P', 'ILE21', 'O_268', 'P', 2.8666, 8.2255],
    ['LEU153', 'N_2341', 'P', 'CYS149', 'O_2286', 'P', 2.8707, 12.4931],
    ['SER7', 'N_38', 'P', 'ASP86', 'OD2_1270', 'P', 2.8732, 39.8838],
    ['ASP86', 'N_1261', 'P', 'SER7', 'OG_45', 'P', 2.8753, 34.7426],
    ['ARG58', 'NH2_832', 'P', 'TYR131', 'O_2016', 'P', 2.8815, 33.1098],
    ['THR46', 'N_644', 'P', 'CYS12', 'O_130', 'P', 2.883, 36.1279],
    ['GLN144', 'N_2197', 'P', 'THR140', 'O_2142', 'P', 2.8836, 23.2545],
    ['THR78', 'N_1134', 'P', 'ASP81', 'OD2_1194', 'P', 2.8869, 12.4465],
    ['LEU89', 'N_1313', 'P', 'LEU9', 'O_83', 'P', 2.8946, 29.5105],
    ['THR31', 'N_421', 'P', 'ARG27', 'O_363', 'P', 2.896, 24.1287],
    ['CYS90', 'N_1332', 'P', 'GLU114', 'O_1738', 'P', 2.8975, 24.576],
    ['CYS148', 'N_2265', 'P', 'GLN144', 'O_2213', 'P', 2.8976, 9.3165],
    ['GLU23', 'N_279', 'P', 'SER19', 'O_235', 'P', 2.8979, 15.4146],
    ['ILE68', 'N_970', 'P', 'MET63', 'O_899', 'P', 2.8986, 12.9904],
    ['PHE10', 'N_84', 'P', 'ASP42', 'O_612', 'P', 2.9026, 22.751],
    ['LYS112', 'N_1683', 'P', 'ASP86', 'O_1272', 'P', 2.912, 10.1158],
    ['SER61', 'N_861', 'P', 'TYR57', 'O_812', 'P', 2.9132, 35.1196],
    ['CYS145', 'N_2214', 'P', 'VAL141', 'O_2158', 'P', 2.9144, 15.8507],
    ['ARG27', 'NH2_359', 'P', 'GLU23', 'OE2_291', 'P', 2.9199, 31.5487],
    ['LYS64', 'N_900', 'P', 'GLN60', 'O_860', 'P', 2.9211, 22.8782],
    ['LEU9', 'N_65', 'P', 'TYR87', 'O_1293', 'P', 2.9229, 16.439],
    ['ASN38', 'N_523', 'P', 'ILE35', 'O_496', 'P', 2.9255, 29.091],
    ['VAL11', 'N_104', 'P', 'LEU89', 'O_1331', 'P', 2.9316, 29.7192],
    ['ASN100', 'N_1486', 'P', 'LEU96', 'O_1430', 'P', 2.933, 10.3321],
    ['GLN124', 'N_1881', 'P', 'ASP120', 'OD2_1825', 'P', 2.9333, 27.4547],
    ['LYS102', 'N_1524', 'P', 'ASP98', 'O_1466', 'P', 2.9361, 9.2855],
    ['GLN76', 'NE2_1110', 'P', 'THR46', 'O_657', 'P', 2.9381, 31.3836],
    ['ARG40', 'NH1_577', 'P', 'THR84', 'OG1_1233', 'P', 2.9482, 8.3748],
    ['ALA44', 'N_624', 'P', 'PHE10', 'O_103', 'P', 2.9499, 33.1772],
    ['GLU154', 'N_2360', 'P', 'ARG150', 'O_2310', 'P', 2.956, 22.5898],
    ['VAL8', 'N_49', 'P', 'ARG40', 'O_584', 'P', 2.9631, 25.0079],
    ['MET63', 'N_883', 'P', 'GLY59', 'O_843', 'P', 2.9733, 18.2731],
    ['GLN60', 'N_844', 'P', 'ASP56', 'O_791', 'P', 2.9795, 35.5229],
    ['ILE35', 'N_478', 'P', 'VAL30', 'O_420', 'P', 2.9811, 23.5092],
    ['VAL146', 'N_2225', 'P', 'TYR142', 'O_2179', 'P', 2.9914, 31.4798],
    ['ARG58', 'NH1_829', 'P', 'TYR131', 'O_2016', 'P', 2.9942, 38.0937],
    ['ASN53', 'N_738', 'P', 'GLU50', 'O_711', 'P', 2.995, 28.587],
    ['ARG101', 'N_1500', 'P', 'ARG97', 'O_1454', 'P', 2.9952, 32.2712],
    ['ARG18', 'NH1_217', 'P', 'ILE127', 'O_1954', 'P', 2.9957, 25.9507],
    ['ARG75', 'N_1074', 'P', 'ASN15', 'OD1_165', 'P', 3.0026, 25.0853],
    ['GLN144', 'NE2_2209', 'P', 'ILE126', 'O_1935', 'P', 3.0038, 18.2744],
    ['ASN34', 'N_464', 'P', 'THR31', 'O_434', 'P', 3.0041, 18.2465],
    ['ASN15', 'ND2_166', 'P', 'SER43', 'OG_620', 'P', 3.0129, 25.6996],
    ['ARG58', 'NE_826', 'P', 'ASP56', 'OD1_788', 'P', 3.017, 22.2284],
    ['ARG27', 'NH1_356', 'P', 'GLU23', 'OE2_291', 'P', 3.0175, 36.9343],
    ['ILE127', 'N_1936', 'P', 'MET91', 'O_1359', 'P', 3.018, 17.5601],
    ['TYR119', 'OH_1808', 'P', 'HSE157', 'N_2407', 'P', 3.0224, 28.0923],
    ['HSE157', 'N_2407', 'P', 'TYR119', 'OH_1808', 'P', 3.0224, 19.1804],
    ['GLU139', 'N_2114', 'P', 'ASP135', 'O_2070', 'P', 3.0245, 27.9246],
    ['LEU29', 'N_386', 'P', 'VAL25', 'O_319', 'P', 3.0299, 19.109],
    ['SER47', 'N_658', 'P', 'LEU13', 'O_149', 'P', 3.0386, 28.8029],
    ['VAL30', 'N_405', 'P', 'PHE26', 'O_339', 'P', 3.0394, 17.6883],
    ['GLN105', 'N_1571', 'P', 'LYS102', 'O_1545', 'P', 3.0464, 19.6806],
    ['SER118', 'N_1784', 'P', 'LEU115', 'O_1757', 'P', 3.051, 21.4045],
    ['LYS155', 'N_2375', 'P', 'ALA151', 'O_2320', 'P', 3.0555, 21.3245],
    ['GLU114', 'N_1724', 'P', 'ILE88', 'O_1312', 'P', 3.059, 24.1605],
    ['ASP120', 'N_1816', 'P', 'GLY117', 'O_1783', 'P', 3.0623, 12.6661],
    ['CYS62', 'N_872', 'P', 'ARG58', 'O_836', 'P', 3.0651, 20.443],
    ['ARG18', 'NH1_217', 'P', 'ASP92', 'OD2_1369', 'P', 3.0679, 4.2778],
    ['ALA24', 'N_294', 'P', 'PRO20', 'O_249', 'P', 3.0751, 29.9487],
    ['ARG150', 'N_2287', 'P', 'VAL146', 'O_2240', 'P', 3.078, 12.7022],
    ['LYS28', 'N_364', 'P', 'ALA24', 'O_303', 'P', 3.0783, 19.9504],
    ['VAL141', 'N_2143', 'P', 'ASP137', 'O_2093', 'P', 3.081, 18.4812],
    ['ASP98', 'N_1455', 'P', 'SER94', 'O_1397', 'P', 3.0844, 19.56],
    ['LEU96', 'N_1412', 'P', 'ASP92', 'O_1371', 'P', 3.085, 36.3254],
    ['ALA22', 'N_269', 'P', 'ARG18', 'O_224', 'P', 3.088, 21.873],
    ['ALA151', 'N_2311', 'P', 'ARG147', 'O_2264', 'P', 3.0991, 15.5713],
    ['GLY67', 'N_963', 'P', 'LYS64', 'O_921', 'P', 3.122, 22.7833],
    ['ASP42', 'N_601', 'P', 'VAL8', 'O_64', 'P', 3.1331, 35.5671],
    ['ARG65', 'N_922', 'P', 'SER61', 'O_871', 'P', 3.1339, 23.3682],
    ['TRP39', 'N_537', 'P', 'SER36', 'O_507', 'P', 3.1343, 15.1776],
    ['LYS123', 'N_1859', 'P', 'ASP120', 'O_1827', 'P', 3.1375, 18.6589],
    ['MET91', 'N_1343', 'P', 'ASN95', 'OD1_1406', 'P', 3.1581, 39.0427],
    ['THR140', 'N_2129', 'P', 'SER136', 'O_2081', 'P', 3.1742, 30.2937],
    ['PHE85', 'N_1241', 'P', 'ASP81', 'O_1196', 'P', 3.1845, 20.2243],
    ['ASN15', 'N_157', 'P', 'CYS12', 'SG_127', 'P', 3.2043, 37.4576],
    ['ALA111', 'N_1673', 'P', 'PHE82', 'O_1216', 'P', 3.2054, 20.58],
    ['ARG147', 'N_2241', 'P', 'GLN143', 'O_2196', 'P', 3.2416, 12.0678],
    ['ARG75', 'NH2_1093', 'P', 'ASP81', 'OD1_1193', 'P', 3.2447, 29.3403],
    ['LYS112', 'NZ_1699', 'P', 'HSE157', 'OT2_2424', 'P', 3.2687, 28.6743],
    ['ARG147', 'NH1_2257', 'P', 'GLN124', 'OE1_1892', 'P', 3.3008, 29.8529],
    ['PHE138', 'N_2094', 'P', 'ASN134', 'O_2058', 'P', 3.3062, 31.0247],
    ['SER7', 'OG_45', 'P', 'THR84', 'O_1240', 'P', 3.3227, 35.5232],
    ['CYS12', 'N_120', 'P', 'ALA44', 'O_633', 'P', 3.3349, 36.1006],
    ['SER19', 'N_225', 'P', 'CYS12', 'SG_127', 'P', 3.339, 8.0034],
    ['PHE82', 'N_1197', 'P', 'LYS79', 'O_1169', 'P', 3.3527, 37.7265],
    ['ASP81', 'N_1185', 'P', 'THR78', 'OG1_1140', 'P', 3.4526, 39.5114],
    ['LYS102', 'NZ_1540', 'P', 'ASP98', 'OD2_1464', 'P', 3.4548, 26.1223],
    ['ARG147', 'NH2_2260', 'P', 'GLN124', 'OE1_1892', 'P', 3.4691, 33.8944],
    ['VAL106', 'N_1588', 'P', 'SER103', 'O_1556', 'P', 3.4974, 34.2367]]

:meth:`.Interactions.getSaltBridges` - salt bridges (residues with oposite
charges):

.. ipython:: python
   :verbatim:   

   interactions.getSaltBridges()

.. parsed-literal::

   [['HSE66', 'NE2_957', 'P', 'GLU139', 'OE1_2125_2126', 'P', 2.8359],
    ['ASP81', 'OD1_1193_1194', 'P', 'ARG75', 'NH1_1090_1093', 'P', 2.9163],
    ['ASP32', 'OD1_443_444', 'P', 'LYS28', 'NZ_380', 'P', 3.037],
    ['ARG101', 'NH1_1516_1519', 'P', 'ASP98', 'OD1_1463_1464', 'P', 3.0699],
    ['ARG27', 'NH1_356_359', 'P', 'GLU23', 'OE1_290_291', 'P', 3.7148],
    ['GLU139', 'OE1_2125_2126', 'P', 'ARG65', 'NH1_938_941', 'P', 3.7799],
    ['LYS102', 'NZ_1540', 'P', 'ASP98', 'OD1_1463_1464', 'P', 3.9359],
    ['ARG58', 'NH1_829_832', 'P', 'ASP56', 'OD1_788_789', 'P', 3.9486],
    ['ARG18', 'NH1_217_220', 'P', 'ASP92', 'OD1_1368_1369', 'P', 4.0693],
    ['GLU114', 'OE1_1735_1736', 'P', 'LYS112', 'NZ_1699', 'P', 4.0787],
    ['ASP120', 'OD1_1824_1825', 'P', 'ARG147', 'NH1_2257_2260', 'P', 4.1543],
    ['LYS110', 'NZ_1667', 'P', 'ASP86', 'OD1_1269_1270', 'P', 4.1879],
    ['GLU114', 'OE1_1735_1736', 'P', 'HSE157', 'NE2_2418', 'P', 4.3835],
    ['ARG18', 'NH1_217_220', 'P', 'ASP129', 'OD1_1978_1979', 'P', 4.5608],
    ['ARG75', 'NH1_1090_1093', 'P', 'ASP42', 'OD1_609_610', 'P', 4.5612],
    ['GLU23', 'OE1_290_291', 'P', 'HSE72', 'NE2_1042', 'P', 4.99]]

:meth:`.Interactions.getRepulsiveIonicBonding` - repulsive ionic bonding
(between residues with the same charges):

.. ipython:: python
   :verbatim:

   interactions.getRepulsiveIonicBonding()

.. parsed-literal::

   [['ARG101', 'NH1_1516_1519', 'P', 'LYS102', 'NZ_1540', 'P', 4.2655]]

:meth:`.Interactions.getPiStacking` - Pi-stacking interactions:

.. ipython:: python
   :verbatim:

   interactions.getPiStacking()

.. parsed-literal::

   [['HSE66',
     '953_954_955_957_959',
     'P',
     'TYR142',
     '2166_2167_2169_2171_2174_2176',
     'P',
     3.8882,
     162.1245],
    ['HSE157',
     '2414_2415_2416_2418_2420_2423_2424',
     'P',
     'TYR119',
     '1802_1803_1805_1807_1810_1812',
     'P',
     4.3605,
     3.0062],
    ['TRP39',
     '549_550_551_553_555_557',
     'P',
     'PHE26',
     '327_328_330_332_334_336',
     'P',
     4.8394,
     75.4588],
    ['TYR131',
     '2003_2004_2006_2008_2011_2013',
     'P',
     'TYR132',
     '2024_2025_2027_2029_2032_2034',
     'P',
     4.8732,
     91.4358]]

:meth:`.Interactions.getPiCation` - Pi-cation:

.. ipython:: python
   :verbatim:

   interactions.getPiCation()

.. parsed-literal::

   [['PHE85',
     '1248_1249_1251_1253_1255_1257',
     'P',
     'ARG40',
     'NH1_577_580',
     'P',
     3.6523],
    ['HSE66', '953_954_955_957_959', 'P', 'ARG65', 'NH1_938_941', 'P', 4.5323],
    ['HSE157',
     '2414_2415_2416_2418_2420_2423_2424',
     'P',
     'LYS112',
     'NZ_1699',
     'P',
     4.828]]

:meth:`.Interactions.getHydrophohic` - hydrophobic interactions:

.. ipython:: python
   :verbatim:

   interactions.getHydrophohic()

.. parsed-literal::

   [['TYR87', 'OH_1286', 'P', 'ALA156', 'CB_2401', 'P', 3.0459],
    ['MET63', 'CE_894', 'P', 'ALA24', 'CB_298', 'P', 3.3105],
    ['ILE68', 'CG2_976', 'P', 'MET63', 'CE_894', 'P', 3.3306],
    ['TYR142', 'CZ_2171', 'P', 'VAL146', 'CG2_2235', 'P', 3.4815],
    ['PHE10', 'CD1_92', 'P', 'ALA22', 'CB_273', 'P', 3.5334],
    ['LYS6', 'CD_26', 'P', 'TRP39', 'CZ2_555', 'P', 3.5427],
    ['VAL30', 'CG1_411', 'P', 'PHE26', 'CE2_336', 'P', 3.5603],
    ['ALA111', 'CB_1677', 'P', 'ILE88', 'CD_1307', 'P', 3.5627],
    ['VAL11', 'CG2_114', 'P', 'ILE88', 'CG2_1300', 'P', 3.6386],
    ['VAL41', 'CG2_595', 'P', 'PHE26', 'CD2_334', 'P', 3.6448],
    ['PHE152', 'CE1_2331', 'P', 'ALA156', 'CB_2401', 'P', 3.6594],
    ['VAL106', 'CG2_1598', 'P', 'LYS79', 'CG_1155', 'P', 3.6828],
    ['ILE77', 'CD_1128', 'P', 'LEU99', 'CD2_1480', 'P', 3.6917],
    ['PHE82', 'CD1_1205', 'P', 'ILE88', 'CD_1307', 'P', 3.692],
    ['LEU116', 'CD2_1771', 'P', 'ILE127', 'CD_1949', 'P', 3.7057],
    ['VAL8', 'CG1_55', 'P', 'PHE26', 'CE2_336', 'P', 3.7106],
    ['LEU125', 'CD1_1907', 'P', 'LEU115', 'CD1_1748', 'P', 3.7115],
    ['MET70', 'CE_1014', 'P', 'MET63', 'CG_890', 'P', 3.7262],
    ['LEU96', 'CD1_1421', 'P', 'ILE113', 'CG2_1711', 'P', 3.7263],
    ['LEU9', 'CD2_78', 'P', 'ILE77', 'CD_1128', 'P', 3.745],
    ['LEU89', 'CD1_1322', 'P', 'VAL8', 'CG2_59', 'P', 3.7672],
    ['ILE126', 'CD_1930', 'P', 'LEU125', 'CD1_1907', 'P', 3.7885],
    ['VAL141', 'CG1_2149', 'P', 'ILE127', 'CG2_1942', 'P', 3.8659],
    ['MET91', 'SD_1353', 'P', 'ILE127', 'CD_1949', 'P', 3.8864],
    ['ALA44', 'CB_628', 'P', 'LEU9', 'CD1_74', 'P', 3.8992],
    ['VAL25', 'CG2_314', 'P', 'TYR142', 'CE1_2169', 'P', 3.92],
    ['ILE21', 'CG2_256', 'P', 'MET63', 'SD_893', 'P', 3.9614],
    ['LEU153', 'CD1_2350', 'P', 'TRP39', 'NE1_547', 'P', 3.967],
    ['PHE85', 'CZ_1253', 'P', 'LEU9', 'CD1_74', 'P', 4.0119],
    ['ILE35', 'CD_491', 'P', 'TRP39', 'NE1_547', 'P', 4.0172],
    ['LEU29', 'CD1_395', 'P', 'VAL25', 'CG1_310', 'P', 4.0642],
    ['ALA74', 'CB_1068', 'P', 'ILE16', 'CG2_177', 'P', 4.0772],
    ['ARG75', 'CG_1081', 'P', 'ALA44', 'CB_628', 'P', 4.0853],
    ['ARG18', 'CG_208', 'P', 'VAL141', 'CG1_2149', 'P', 4.104],
    ['LYS102', 'CD_1534', 'P', 'ILE77', 'CG2_1121', 'P', 4.1048],
    ['TYR119', 'CE1_1805', 'P', 'LEU89', 'CD2_1326', 'P', 4.1435],
    ['ARG40', 'CG_568', 'P', 'PHE85', 'CE2_1257', 'P', 4.2669],
    ['LYS28', 'CG_371', 'P', 'ILE68', 'CD_983', 'P', 4.2707],
    ['PHE138', 'CD2_2108', 'P', 'ILE21', 'CD_263', 'P', 4.3082],
    ['LYS112', 'CG_1690', 'P', 'TYR87', 'CE1_1283', 'P', 4.3083],
    ['TYR131', 'CE1_2006', 'P', 'ILE16', 'CD_184', 'P', 4.3352],
    ['ARG58', 'CG_820', 'P', 'PHE138', 'CE1_2104', 'P', 4.4781]]

:meth:`.Interactions.getDisulfideBonds` - disulfide bonds:

.. ipython:: python
   :verbatim:

   interactions.getDisulfideBonds()

.. parsed-literal::

   []

To display residues with the biggest number of potential interactions and their
types, we can use :meth:`.Interactions.getFrequentInteractions` method:

.. ipython:: python
   :verbatim:

   frequent_interactions = interactions.getFrequentInteractions(contacts_min=3)
   frequent_interactions

.. parsed-literal::

   @> The most frequent interactions between:
   @> LEU9P  <--->  hp:ALA44P  hp:PHE85P  hb:LEU89P
   @> CYS12P  <--->  hb:ASN15P  hb:SER19P  hb:THR46P
   @> ILE16P  <--->  hb:ASP56P  hp:ALA74P  hp:TYR131P
   @> PHE26P  <--->  hp:VAL8P  hp:VAL30P  ps:TRP39P  hp:VAL41P
   @> TRP39P  <--->  hp:LYS6P  hp:ILE35P  hp:LEU153P
   @> MET63P  <--->  hp:ILE21P  hp:ILE68P  hp:MET70P
   @> ASP81P  <--->  hb:ARG75P  hb:THR78P  hb:PHE85P
   @> THR84P  <--->  hb:SER7P  hb:ARG40P  hb:LYS110P
   @> ASP86P  <--->  hb:SER7P  sb:LYS110P  hb:LYS112P
   @> ILE88P  <--->  hp:VAL11P  hp:PHE82P  hp:ALA111P  hb:GLU114P
   @> ASP92P  <--->  sb:ARG18P  hb:ASN95P  hb:LEU96P
   @> LYS112P  <--->  hb:ILE88P  sb:GLU114P  pc:HSE157P
   @> ILE127P  <--->  hb:ARG18P  hp:MET91P  hp:LEU116P  hp:VAL141P
   @> Legend: hb-hydrogen bond, sb-salt bridge, rb-repulsive ionic bond, ps-Pi stacking interaction,pc-Cation-Pi interaction, hp-hydrophobic interaction, dibs-disulfide bonds
   @> The biggest number of interactions: 4

The value of *contacts_min* can be modified to display residues with smaller
number of interactions. 


Visualize interactions in VMD
-------------------------------------------------------------------------------

We can generate tcl files for visualizing each type of interaction with VMD_ 
using the :func:`.showProteinInteractions_VMD` function in the following way:

.. ipython:: python
   :verbatim:

   showProteinInteractions_VMD(atoms, interactions.getHydrogenBonds(), color='blue', filename='HBs.tcl')
   showProteinInteractions_VMD(atoms, interactions.getSaltBridges(), color='yellow',filename='SBs.tcl')
   showProteinInteractions_VMD(atoms, interactions.getRepulsiveIonicBonding(), color='red',filename='RIB.tcl')
   showProteinInteractions_VMD(atoms, interactions.getPiStacking(), color='green',filename='PiStacking.tcl') 
   showProteinInteractions_VMD(atoms, interactions.getPiCation(), color='orange',filename='PiCation.tcl') 
   showProteinInteractions_VMD(atoms, interactions.getHydrophobic(), color='silver',filename='HPh.tcl')
   showProteinInteractions_VMD(atoms, interactions.getDisulfideBonds(), color='black',filename='DiBs.tcl') 

.. parsed-literal::

   @> TCL file saved
   @> TCL file saved
   @> TCL file saved
   @> TCL file saved
   @> TCL file saved
   @> TCL file saved
   @> Lack of results
   @> TCL file saved

A TCL file will be saved and can be used in VMD_ after uploading the PDB file
with protein structure **5kqm_all_sci.pdb** and by running the following command 
line instruction in the VMD_ *TKConsole* (*VMD Main*) for Linux, Windows and Mac users: 

::  play HBs.tcl

The tcl file contains a method for drawing lines between selected pairs of 
residues. Those residues are also displayed.

.. figure:: images/HBs.png
   :scale: 60 %


::  play SBs.tcl

.. figure:: images/SBs.png
   :scale: 60 %


::  play RIB.tcl

.. figure:: images/RIB.png
   :scale: 60 %


::  play PiStacking.tcl

.. figure:: images/PiStacking.png
   :scale: 60 %


::  play PiCation.tcl

.. figure:: images/PiCation.png
   :scale: 60 %


::  play HPh.tcl

.. figure:: images/Hydrophobic.png
   :scale: 60 %


Additional selections
-------------------------------------------------------------------------------

From the predicted interactions we can select only interactions assigned to
a certain regions, chains or between different chains.

We can compute them by adding additional parameters to the selected
function. See examples below:

.. ipython:: python
   :verbatim:

   interactions.getSaltBridges(selection='chain P')

.. parsed-literal::

   [['HSE66', 'NE2_957', 'P', 'GLU139', 'OE1_2125_2126', 'P', 2.8359],
    ['ASP81', 'OD1_1193_1194', 'P', 'ARG75', 'NH1_1090_1093', 'P', 2.9163],
    ['ASP32', 'OD1_443_444', 'P', 'LYS28', 'NZ_380', 'P', 3.037],
    ['ARG101', 'NH1_1516_1519', 'P', 'ASP98', 'OD1_1463_1464', 'P', 3.0699],
    ['ARG27', 'NH1_356_359', 'P', 'GLU23', 'OE1_290_291', 'P', 3.7148],
    ['GLU139', 'OE1_2125_2126', 'P', 'ARG65', 'NH1_938_941', 'P', 3.7799],
    ['LYS102', 'NZ_1540', 'P', 'ASP98', 'OD1_1463_1464', 'P', 3.9359],
    ['ARG58', 'NH1_829_832', 'P', 'ASP56', 'OD1_788_789', 'P', 3.9486],
    ['ARG18', 'NH1_217_220', 'P', 'ASP92', 'OD1_1368_1369', 'P', 4.0693],
    ['GLU114', 'OE1_1735_1736', 'P', 'LYS112', 'NZ_1699', 'P', 4.0787],
    ['ASP120', 'OD1_1824_1825', 'P', 'ARG147', 'NH1_2257_2260', 'P', 4.1543],
    ['LYS110', 'NZ_1667', 'P', 'ASP86', 'OD1_1269_1270', 'P', 4.1879],
    ['GLU114', 'OE1_1735_1736', 'P', 'HSE157', 'NE2_2418', 'P', 4.3835],
    ['ARG18', 'NH1_217_220', 'P', 'ASP129', 'OD1_1978_1979', 'P', 4.5608],
    ['ARG75', 'NH1_1090_1093', 'P', 'ASP42', 'OD1_609_610', 'P', 4.5612],
    ['GLU23', 'OE1_290_291', 'P', 'HSE72', 'NE2_1042', 'P', 4.99]]

.. ipython:: python
   :verbatim:

   interactions.getRepulsiveIonicBonding(selection='resid 102')

.. parsed-literal::

   [['ARG101', 'NH1_1516_1519', 'P', 'LYS102', 'NZ_1540', 'P', 4.2655]]

.. ipython:: python
   :verbatim:

   interactions.getPiStacking(selection='chain P and resid 26')

.. parsed-literal::

   [['TRP39',
  '549_550_551_553_555_557',
  'P',
  'PHE26',
  '327_328_330_332_334_336',
  'P',
  4.8394,
  75.4588]]

It can be done for all kinds of interactions as well. The function will
return a list of interactions with following order:

    (1) Hydrogen bonds
    (2) Salt Bridges
    (3) RepulsiveIonicBonding 
    (4) Pi stacking interactions
    (5) Pi-cation interactions
    (6) Hydrophobic interactions
    (7) Disulfide bonds

.. ipython:: python
   :verbatim:

   allRes_20to50 = interactions.getInteractions(selection='resid 20 to 50')
   allRes_20to50

.. parsed-literal::

   [[['ARG40', 'N_561', 'P', 'LYS6', 'O_37', 'P', 2.7479, 17.1499],
     ['ALA45', 'N_634', 'P', 'ARG75', 'O_1097', 'P', 2.7609, 35.0983],
     ['ASN53', 'ND2_747', 'P', 'GLU50', 'OE1_708', 'P', 2.7702, 18.2336],
     ['LYS6', 'N_16', 'P', 'ASN38', 'O_536', 'P', 2.8178, 25.0305],
     ['ILE77', 'N_1115', 'P', 'ALA45', 'O_643', 'P', 2.8179, 12.1855],
     ['ASP32', 'N_435', 'P', 'LYS28', 'O_385', 'P', 2.8357, 8.8318],
     ['ARG27', 'N_340', 'P', 'GLU23', 'O_293', 'P', 2.8446, 15.4167],
     ['PHE26', 'N_320', 'P', 'ALA22', 'O_278', 'P', 2.8541, 4.8732],
     ['ARG75', 'NH2_1093', 'P', 'ASP42', 'OD2_610', 'P', 2.8649, 23.5083],
     ['VAL25', 'N_304', 'P', 'ILE21', 'O_268', 'P', 2.8666, 8.2255],
     ['THR46', 'N_644', 'P', 'CYS12', 'O_130', 'P', 2.883, 36.1279],
     ['THR31', 'N_421', 'P', 'ARG27', 'O_363', 'P', 2.896, 24.1287],
     ['GLU23', 'N_279', 'P', 'SER19', 'O_235', 'P', 2.8979, 15.4146],
     ['PHE10', 'N_84', 'P', 'ASP42', 'O_612', 'P', 2.9026, 22.751],
     ['ARG27', 'NH2_359', 'P', 'GLU23', 'OE2_291', 'P', 2.9199, 31.5487],
     ['ASN38', 'N_523', 'P', 'ILE35', 'O_496', 'P', 2.9255, 29.091],
     ['GLN76', 'NE2_1110', 'P', 'THR46', 'O_657', 'P', 2.9381, 31.3836],
     ['ARG40', 'NH1_577', 'P', 'THR84', 'OG1_1233', 'P', 2.9482, 8.3748],
     ['ALA44', 'N_624', 'P', 'PHE10', 'O_103', 'P', 2.9499, 33.1772],
     ['VAL8', 'N_49', 'P', 'ARG40', 'O_584', 'P', 2.9631, 25.0079],
     ['ILE35', 'N_478', 'P', 'VAL30', 'O_420', 'P', 2.9811, 23.5092],
     ['ASN53', 'N_738', 'P', 'GLU50', 'O_711', 'P', 2.995, 28.587],
     ['ASN34', 'N_464', 'P', 'THR31', 'O_434', 'P', 3.0041, 18.2465],
     ['ASN15', 'ND2_166', 'P', 'SER43', 'OG_620', 'P', 3.0129, 25.6996],
     ['ARG27', 'NH1_356', 'P', 'GLU23', 'OE2_291', 'P', 3.0175, 36.9343],
     ['LEU29', 'N_386', 'P', 'VAL25', 'O_319', 'P', 3.0299, 19.109],
     ['SER47', 'N_658', 'P', 'LEU13', 'O_149', 'P', 3.0386, 28.8029],
     ['VAL30', 'N_405', 'P', 'PHE26', 'O_339', 'P', 3.0394, 17.6883],
     ['ALA24', 'N_294', 'P', 'PRO20', 'O_249', 'P', 3.0751, 29.9487],
     ['LYS28', 'N_364', 'P', 'ALA24', 'O_303', 'P', 3.0783, 19.9504],
     ['ALA22', 'N_269', 'P', 'ARG18', 'O_224', 'P', 3.088, 21.873],
     ['ASP42', 'N_601', 'P', 'VAL8', 'O_64', 'P', 3.1331, 35.5671],
     ['TRP39', 'N_537', 'P', 'SER36', 'O_507', 'P', 3.1343, 15.1776],
     ['CYS12', 'N_120', 'P', 'ALA44', 'O_633', 'P', 3.3349, 36.1006]],
    [['ASP32', 'OD1_443_444', 'P', 'LYS28', 'NZ_380', 'P', 3.037],
     ['ARG27', 'NH1_356_359', 'P', 'GLU23', 'OE1_290_291', 'P', 3.7148],
     ['ARG75', 'NH1_1090_1093', 'P', 'ASP42', 'OD1_609_610', 'P', 4.5612],
     ['GLU23', 'OE1_290_291', 'P', 'HSE72', 'NE2_1042', 'P', 4.99]],
    [],
    [['TRP39',
      '549_550_551_553_555_557',
      'P',
      'PHE26',
      '327_328_330_332_334_336',
      'P',
      4.8394,
      75.4588]],
    [['PHE85',
      '1248_1249_1251_1253_1255_1257',
      'P',
      'ARG40',
      'NH1_577_580',
      'P',
      3.6523]],
    [['MET63', 'CE_894', 'P', 'ALA24', 'CB_298', 'P', 3.3105],
     ['PHE10', 'CD1_92', 'P', 'ALA22', 'CB_273', 'P', 3.5334],
     ['LYS6', 'CD_26', 'P', 'TRP39', 'CZ2_555', 'P', 3.5427],
     ['VAL30', 'CG1_411', 'P', 'PHE26', 'CE2_336', 'P', 3.5603],
     ['VAL41', 'CG2_595', 'P', 'PHE26', 'CD2_334', 'P', 3.6448],
     ['VAL8', 'CG1_55', 'P', 'PHE26', 'CE2_336', 'P', 3.7106],
     ['ALA44', 'CB_628', 'P', 'LEU9', 'CD1_74', 'P', 3.8992],
     ['VAL25', 'CG2_314', 'P', 'TYR142', 'CE1_2169', 'P', 3.92],
     ['ILE21', 'CG2_256', 'P', 'MET63', 'SD_893', 'P', 3.9614],
     ['LEU153', 'CD1_2350', 'P', 'TRP39', 'NE1_547', 'P', 3.967],
     ['ILE35', 'CD_491', 'P', 'TRP39', 'NE1_547', 'P', 4.0172],
     ['LEU29', 'CD1_395', 'P', 'VAL25', 'CG1_310', 'P', 4.0642],
     ['ARG75', 'CG_1081', 'P', 'ALA44', 'CB_628', 'P', 4.0853],
     ['ARG40', 'CG_568', 'P', 'PHE85', 'CE2_1257', 'P', 4.2669],
     ['LYS28', 'CG_371', 'P', 'ILE68', 'CD_983', 'P', 4.2707],
     ['PHE138', 'CD2_2108', 'P', 'ILE21', 'CD_263', 'P', 4.3082]],
    []]

The list of hydrogen bonds, salt bridges and other types of interactions can
be displayed as follows:

.. ipython:: python
   :verbatim:

   allRes_20to50[0]

.. parsed-literal::

   [['ARG40', 'N_561', 'P', 'LYS6', 'O_37', 'P', 2.7479, 17.1499],
    ['ALA45', 'N_634', 'P', 'ARG75', 'O_1097', 'P', 2.7609, 35.0983],
    ['ASN53', 'ND2_747', 'P', 'GLU50', 'OE1_708', 'P', 2.7702, 18.2336],
    ['LYS6', 'N_16', 'P', 'ASN38', 'O_536', 'P', 2.8178, 25.0305],
    ['ILE77', 'N_1115', 'P', 'ALA45', 'O_643', 'P', 2.8179, 12.1855],
    ['ASP32', 'N_435', 'P', 'LYS28', 'O_385', 'P', 2.8357, 8.8318],
    ['ARG27', 'N_340', 'P', 'GLU23', 'O_293', 'P', 2.8446, 15.4167],
    ['PHE26', 'N_320', 'P', 'ALA22', 'O_278', 'P', 2.8541, 4.8732],
    ['ARG75', 'NH2_1093', 'P', 'ASP42', 'OD2_610', 'P', 2.8649, 23.5083],
    ['VAL25', 'N_304', 'P', 'ILE21', 'O_268', 'P', 2.8666, 8.2255],
    ['THR46', 'N_644', 'P', 'CYS12', 'O_130', 'P', 2.883, 36.1279],
    ['THR31', 'N_421', 'P', 'ARG27', 'O_363', 'P', 2.896, 24.1287],
    ['GLU23', 'N_279', 'P', 'SER19', 'O_235', 'P', 2.8979, 15.4146],
    ['PHE10', 'N_84', 'P', 'ASP42', 'O_612', 'P', 2.9026, 22.751],
    ['ARG27', 'NH2_359', 'P', 'GLU23', 'OE2_291', 'P', 2.9199, 31.5487],
    ['ASN38', 'N_523', 'P', 'ILE35', 'O_496', 'P', 2.9255, 29.091],
    ['GLN76', 'NE2_1110', 'P', 'THR46', 'O_657', 'P', 2.9381, 31.3836],
    ['ARG40', 'NH1_577', 'P', 'THR84', 'OG1_1233', 'P', 2.9482, 8.3748],
    ['ALA44', 'N_624', 'P', 'PHE10', 'O_103', 'P', 2.9499, 33.1772],
    ['VAL8', 'N_49', 'P', 'ARG40', 'O_584', 'P', 2.9631, 25.0079],
    ['ILE35', 'N_478', 'P', 'VAL30', 'O_420', 'P', 2.9811, 23.5092],
    ['ASN53', 'N_738', 'P', 'GLU50', 'O_711', 'P', 2.995, 28.587],
    ['ASN34', 'N_464', 'P', 'THR31', 'O_434', 'P', 3.0041, 18.2465],
    ['ASN15', 'ND2_166', 'P', 'SER43', 'OG_620', 'P', 3.0129, 25.6996],
    ['ARG27', 'NH1_356', 'P', 'GLU23', 'OE2_291', 'P', 3.0175, 36.9343],
    ['LEU29', 'N_386', 'P', 'VAL25', 'O_319', 'P', 3.0299, 19.109],
    ['SER47', 'N_658', 'P', 'LEU13', 'O_149', 'P', 3.0386, 28.8029],
    ['VAL30', 'N_405', 'P', 'PHE26', 'O_339', 'P', 3.0394, 17.6883],
    ['ALA24', 'N_294', 'P', 'PRO20', 'O_249', 'P', 3.0751, 29.9487],
    ['LYS28', 'N_364', 'P', 'ALA24', 'O_303', 'P', 3.0783, 19.9504],
    ['ALA22', 'N_269', 'P', 'ARG18', 'O_224', 'P', 3.088, 21.873],
    ['ASP42', 'N_601', 'P', 'VAL8', 'O_64', 'P', 3.1331, 35.5671],
    ['TRP39', 'N_537', 'P', 'SER36', 'O_507', 'P', 3.1343, 15.1776],
    ['CYS12', 'N_120', 'P', 'ALA44', 'O_633', 'P', 3.3349, 36.1006]]

Salt Bridges:

.. ipython:: python
   :verbatim:

   allRes_20to50[1]

.. parsed-literal::

   [['ASP32', 'OD1_443_444', 'P', 'LYS28', 'NZ_380', 'P', 3.037],
    ['ARG27', 'NH1_356_359', 'P', 'GLU23', 'OE1_290_291', 'P', 3.7148],
    ['ARG75', 'NH1_1090_1093', 'P', 'ASP42', 'OD1_609_610', 'P', 4.5612],
    ['GLU23', 'OE1_290_291', 'P', 'HSE72', 'NE2_1042', 'P', 4.99]]

We can also select one particular residue of our interest:

.. ipython:: python
   :verbatim:

   interactions.getPiCation(selection='resid 85')

.. parsed-literal::

   [['PHE85',
     '1248_1249_1251_1253_1255_1257',
     'P',
     'ARG40',
     'NH1_577_580',
     'P',
     3.6523]]

.. ipython:: python
   :verbatim:

   interactions.getHydrophobic(selection='resid 26 to 100')

.. parsed-literal::

   [['TYR87', 'OH_1286', 'P', 'ALA156', 'CB_2401', 'P', 3.0459],
    ['MET63', 'CE_894', 'P', 'ALA24', 'CB_298', 'P', 3.3105],
    ['ILE68', 'CG2_976', 'P', 'MET63', 'CE_894', 'P', 3.3306],
    ['LYS6', 'CD_26', 'P', 'TRP39', 'CZ2_555', 'P', 3.5427],
    ['VAL30', 'CG1_411', 'P', 'PHE26', 'CE2_336', 'P', 3.5603],
    ['ALA111', 'CB_1677', 'P', 'ILE88', 'CD_1307', 'P', 3.5627],
    ['VAL11', 'CG2_114', 'P', 'ILE88', 'CG2_1300', 'P', 3.6386],
    ['VAL41', 'CG2_595', 'P', 'PHE26', 'CD2_334', 'P', 3.6448],
    ['VAL106', 'CG2_1598', 'P', 'LYS79', 'CG_1155', 'P', 3.6828],
    ['ILE77', 'CD_1128', 'P', 'LEU99', 'CD2_1480', 'P', 3.6917],
    ['PHE82', 'CD1_1205', 'P', 'ILE88', 'CD_1307', 'P', 3.692],
    ['VAL8', 'CG1_55', 'P', 'PHE26', 'CE2_336', 'P', 3.7106],
    ['MET70', 'CE_1014', 'P', 'MET63', 'CG_890', 'P', 3.7262],
    ['LEU96', 'CD1_1421', 'P', 'ILE113', 'CG2_1711', 'P', 3.7263],
    ['LEU9', 'CD2_78', 'P', 'ILE77', 'CD_1128', 'P', 3.745],
    ['LEU89', 'CD1_1322', 'P', 'VAL8', 'CG2_59', 'P', 3.7672],
    ['MET91', 'SD_1353', 'P', 'ILE127', 'CD_1949', 'P', 3.8864],
    ['ALA44', 'CB_628', 'P', 'LEU9', 'CD1_74', 'P', 3.8992],
    ['ILE21', 'CG2_256', 'P', 'MET63', 'SD_893', 'P', 3.9614],
    ['LEU153', 'CD1_2350', 'P', 'TRP39', 'NE1_547', 'P', 3.967],
    ['PHE85', 'CZ_1253', 'P', 'LEU9', 'CD1_74', 'P', 4.0119],
    ['ILE35', 'CD_491', 'P', 'TRP39', 'NE1_547', 'P', 4.0172],
    ['LEU29', 'CD1_395', 'P', 'VAL25', 'CG1_310', 'P', 4.0642],
    ['ALA74', 'CB_1068', 'P', 'ILE16', 'CG2_177', 'P', 4.0772],
    ['ARG75', 'CG_1081', 'P', 'ALA44', 'CB_628', 'P', 4.0853],
    ['LYS102', 'CD_1534', 'P', 'ILE77', 'CG2_1121', 'P', 4.1048],
    ['TYR119', 'CE1_1805', 'P', 'LEU89', 'CD2_1326', 'P', 4.1435],
    ['ARG40', 'CG_568', 'P', 'PHE85', 'CE2_1257', 'P', 4.2669],
    ['LYS28', 'CG_371', 'P', 'ILE68', 'CD_983', 'P', 4.2707],
    ['LYS112', 'CG_1690', 'P', 'TYR87', 'CE1_1283', 'P', 4.3083],
    ['ARG58', 'CG_820', 'P', 'PHE138', 'CE1_2104', 'P', 4.4781]]


Change selection criteria for interaction type
-------------------------------------------------------------------------------

The :meth:`.Interactions.buildInteractionMatrix` method computes interactions 
using default parameters for interactions. However, it can be changed
according to our needs. To do that, we need to recalculate the selected type
of interactions. 

We can do it using the following functions: :func:`.calcHydrogenBonds`,
:func:`.calcHydrogenBonds`, :func:`.calcSaltBridges`,
:func:`.calcRepulsiveIonicBonding`, :func:`.calcPiStacking`,
:func:`.calcPiCation`, :func:`.calcHydrophohic`,
:func:`.calcDisulfideBonds`, and use
:meth:`.Interactions.setNewHydrogenBonds`,
:meth:`.Interactions.setNewSaltBridges`,
:meth:`.Interactions.setNewRepulsiveIonicBonding`,
:meth:`.Interactions.setNewPiStacking`,
:meth:`.Interactions.setNewPiCation`,
:meth:`.Interactions.setNewHydrophohic`,
:meth:`.Interactions.setNewDisulfideBonds` method to replace it in the main
Instance. 

For example:

.. ipython:: python
   :verbatim:

   newHydrogenBonds2 = calcHydrogenBonds(atoms, distA=2.8, angle=30, cutoff_dist=15)
   interactions.setNewHydrogenBonds(newHydrogenBonds2)

.. parsed-literal::

   @> Calculating hydrogen bonds.
   @>      DONOR (res chid atom)   <--->       ACCEPTOR (res chid atom)    Distance  Angle
   @>     GLN143    P      NE2_2192  <--->     GLU139    P      OE2_2126     2.7     9.2
   @>      HSE66    P       NE2_957  <--->     GLU139    P      OE1_2125     2.7     6.4
   @>      ARG40    P         N_561  <--->       LYS6    P          O_37     2.7    17.1
   @>      ARG58    P         N_813  <--->      ASP56    P       OD1_788     2.7    30.0
   @>      ASN53    P       ND2_747  <--->      GLU50    P       OE1_708     2.8    18.2
   @>      ALA74    P        N_1064  <--->      ASN53    P         O_751     2.8    21.3
   @>      ASP56    P         N_780  <--->      ILE16    P         O_189     2.8    27.0
   @> Number of detected hydrogen bonds: 7.
   @> Hydrogen Bonds are replaced
   
.. ipython:: python
   :verbatim:

   interactions.getHydrogenBonds()

.. parsed-literal::

   [['GLN143', 'NE2_2192', 'P', 'GLU139', 'OE2_2126', 'P', 2.7287, 9.1823],
    ['HSE66', 'NE2_957', 'P', 'GLU139', 'OE1_2125', 'P', 2.7314, 6.3592],
    ['ARG40', 'N_561', 'P', 'LYS6', 'O_37', 'P', 2.7479, 17.1499],
    ['ARG58', 'N_813', 'P', 'ASP56', 'OD1_788', 'P', 2.7499, 29.9737],
    ['ASN53', 'ND2_747', 'P', 'GLU50', 'OE1_708', 'P', 2.7702, 18.2336],
    ['ALA74', 'N_1064', 'P', 'ASN53', 'O_751', 'P', 2.7782, 21.3375],
    ['ASP56', 'N_780', 'P', 'ILE16', 'O_189', 'P', 2.7793, 27.0481]]

.. ipython:: python
   :verbatim:

   sb2 = calcSaltBridges(atoms, distA=6)
   interactions.setNewSaltBridges(sb2)

   rib2 = calcRepulsiveIonicBonding(atoms, distA=9)
   interactions.setNewRepulsiveIonicBonding(rib2)

   picat2 = calcPiCation(atoms, distA=7)
   interactions.setNewPiCation(picat2)

.. parsed-literal::

   @> Calculating salt bridges.
   @>      HSE66    P         NE2_957  <--->     GLU139    P   OE1_2125_2126     2.8
   @>      ASP81    P   OD1_1193_1194  <--->      ARG75    P   NH1_1090_1093     2.9
   @>      ASP32    P     OD1_443_444  <--->      LYS28    P          NZ_380     3.0
   @>     ARG101    P   NH1_1516_1519  <--->      ASP98    P   OD1_1463_1464     3.1
   @>      ARG27    P     NH1_356_359  <--->      GLU23    P     OE1_290_291     3.7
   @>     GLU139    P   OE1_2125_2126  <--->      ARG65    P     NH1_938_941     3.8
   @>     LYS102    P         NZ_1540  <--->      ASP98    P   OD1_1463_1464     3.9
   @>      ARG58    P     NH1_829_832  <--->      ASP56    P     OD1_788_789     3.9
   @>      ARG18    P     NH1_217_220  <--->      ASP92    P   OD1_1368_1369     4.1
   @>     GLU114    P   OE1_1735_1736  <--->     LYS112    P         NZ_1699     4.1
   @>     ASP120    P   OD1_1824_1825  <--->     ARG147    P   NH1_2257_2260     4.2
   @>     LYS110    P         NZ_1667  <--->      ASP86    P   OD1_1269_1270     4.2
   @>     GLU114    P   OE1_1735_1736  <--->     HSE157    P        NE2_2418     4.4
   @>      ARG18    P     NH1_217_220  <--->     ASP129    P   OD1_1978_1979     4.6
   @>      ARG75    P   NH1_1090_1093  <--->      ASP42    P     OD1_609_610     4.6
   @>      GLU23    P     OE1_290_291  <--->      HSE72    P        NE2_1042     5.0
   @>      ASP42    P     OD1_609_610  <--->      HSE72    P        NE2_1042     5.4
   @>      ASP81    P   OD1_1193_1194  <--->      ARG40    P     NH1_577_580     5.8
   @> Number of detected salt bridges: 18.
   @> Salt Bridges are replaced
   @> Calculating repulsive ionic bonding.
   @>      ASP42    P     OD1_609_610  <--->      ASP81    P   OD1_1193_1194     6.7
   @>      GLU80    P   OE1_1181_1182  <--->      ASP81    P   OD1_1193_1194     7.0
   @>      ASP92    P   OD1_1368_1369  <--->     ASP129    P   OD1_1978_1979     7.6
   @>     LYS110    P         NZ_1667  <--->      ARG40    P     NH1_577_580     7.8
   @>      ASP92    P   OD1_1368_1369  <--->      GLU93    P   OE1_1383_1384     8.6
   @>     GLU128    P   OE1_1966_1967  <--->     ASP137    P   OD1_2090_2091     8.9
   @> Number of detected Repulsive Ionic Bonding interactions: 6.
   @> Repulsive Ionic Bonding are replaced
   @> Calculating cation-Pi interactions.
   @>      PHE85   P   1248_1249_1251_1253_1255_1257  <--->      ARG40   P                     NH1_577_580     3.7
   @>      HSE66   P             953_954_955_957_959  <--->      ARG65   P                     NH1_938_941     4.5
   @>     HSE157   P2414_2415_2416_2418_2420_2423_2424  <--->     LYS112   P                         NZ_1699     4.8
   @>     PHE138   P   2101_2102_2104_2106_2108_2110  <--->      ARG58   P                     NH1_829_832     5.0
   @>     TYR131   P   2003_2004_2006_2008_2011_2013  <--->      ARG58   P                     NH1_829_832     5.1
   @>      PHE85   P   1248_1249_1251_1253_1255_1257  <--->      ARG75   P                   NH1_1090_1093     6.3
   @>      TRP39   P         549_550_551_553_555_557  <--->       LYS6   P                           NZ_32     6.6
   @>      TYR87   P   1280_1281_1283_1285_1288_1290  <--->     LYS112   P                         NZ_1699     6.8
   @> Number of detected cation-pi interactions: 8.
   @> Pi-Cation interactions are replaced


Assess the functional significance of a residue
-------------------------------------------------------------------------------

For assessing the functional significance of each residue in protein
structure, we counted the number of possible contacts based on:

    (1) Hydrogen bonds (HBs)
    (2) Salt Bridges (SBs)
    (3) Repulsive Ionic Bonding (RIB)  
    (4) Pi stacking interactions (PiStack)
    (5) Pi-cation interactions (PiCat) 
    (6) Hydrophobic interactions (HPh) 
    (7) Disulfide Bonds (DiBs)


To compute the weighted interactions use the 
:meth:`.Interactions.buildInteractionMatrix` method:

.. ipython:: python
   :verbatim:

   matrix = interactions.buildInteractionMatrix()

.. parsed-literal::

   @> Calculating interactions

The results can be displayed in the following way:

.. ipython:: python
   :verbatim:

    import matplotlib.pylab as plt
    showAtomicMatrix(matrix, atoms=atoms.ca, cmap='seismic', markersize=8)
    plt.xlabel('Residue')
    plt.ylabel('Residue')
    plt.clim([-3,3])

.. figure:: images/single_imshow.png
   :scale: 60 %

The total number of interaction for each residue can be displayed on the plot using
:func:`.showCumulativeInteractionTypes()` function.

.. ipython:: python
   :verbatim:

   interactions.showCumulativeInteractionTypes()

.. parsed-literal::

   @> Calculating interactions
   @> Calculating interactions
   @> Calculating interactions
   @> Calculating interactions
   @> Calculating interactions
   @> Calculating interactions
   @> Calculating interactions

The results with the higest number of possible contacts can be saved in PDB
file. They will be restored in Occupancy column and display in VMD_.

.. figure:: images/single_bar_plot.png
   :scale: 60 %

.. ipython:: python
   :verbatim:

   interactions.saveInteractionsPDB(filename='5kqm_meanMatrix.pdb')

.. parsed-literal::


Visualize number of interactions onto 3D structure
-------------------------------------------------------------------------------

The number of the interaction can be saved to a PDB file in the
*Occupancy* column by using :meth:`.Interactions.saveInteractionsPDB`
method. Then the score would be displayed in color in any available graphical
program, for example, in VMD_.
 
.. ipython:: python
   :verbatim:

   interactions.saveInteractionsPDB(filename='5kqm_meanMatrix.pdb')

.. parsed-literal::

   @> PDB file saved.

A file *5kqm_meanMatrix.pdb* will be saved and can be used in VMD_ by 
uploading PDB structure and displaying it with *Coloring Method*
*Occupancy*. By default blue colors correspond to the highest values but we
can change it in *VMD Main* -> *Graphics* -> *Color Controls* -> *Color
Scale* -> *Method* to *BWR*.

.. figure:: images/fig1.png
   :scale: 60 %


Exclude some interaction types from calculations
-------------------------------------------------------------------------------

For analysis we can exclude some of the interaction types by assigning zero
to the type of interactions (HBs - hydrogen bonds, SBs - salt bridges, RIB -
repulsive ionic bonding, PiCat - Pi-Cation, PiStack - Pi-Stacking, HPh -
hydrophobic interactions and finally DiBs - disulfide bonds). 

.. ipython:: python
   :verbatim:

   matrix = interactions.buildInteractionMatrix(RIB=0, HBs=0, HPh=0, DiBs=0)

.. parsed-literal::

   @> Calculating interactions

The results can be displayed in a similar way:

.. ipython:: python
   :verbatim:

    showAtomicMatrix(matrix, atoms=atoms.ca, cmap='seismic', markersize=8)
    plt.xlabel('Residue')
    plt.ylabel('Residue')
    plt.clim([-3,3])


.. figure:: images/single_imshow2.png
   :scale: 60 %
