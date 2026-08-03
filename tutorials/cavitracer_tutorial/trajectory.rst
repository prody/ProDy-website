.. _cavitracer_single:

Detection of channels in molecular dynamics (MD) trajectory
===============================================================================

Analysis of the trajectory will be performed on a short MD trajectory
containing over 200 frames of simulation performed for the vesicular monoamine 
transporter VMAT2. This protein contains 460 residues. During the analysis, 
only protein structure will be taken into consideration.

Before analyzing the trajectory, its need to be parsed (see more details
in `Trajectory tutorial`_).

.. ipython:: python
   :verbatim:

   PDBfile = 'caseStudy2.pdb'
   DCDfile = 'caseStudy2.dcd'
   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)

.. parsed-literal::

   @> 56457 atoms and 1 coordinate set(s) were parsed in 0.69s.

To analyze only protein structure for the analysis, we will select protein
structure:

.. ipython:: python
   :verbatim:

   protein = atoms.select("protein")
   dcd.setAtoms(protein)

Next, to detect channels in MD trajectories, we need to use
:func:`.scalcChannelsMultipleFrames` function. To speed up the calculations,
``max_proc`` parameter can be set up to a higher number. In the example
below, we will use 4 processors to perform the calulations. The results will
saved with the prefix ``"chls_dcd"``, and all the channels will be storage
separately (``separate=True``). 

.. ipython:: python
   :verbatim:

   channels4, surfaces4=calcChannelsMultipleFrames(atoms, dcd, 
		output_path='chls_dcd', separate=True, max_proc=4)

.. parsed-literal::

   @> Frame/model: 0
   @> Frame/model: 14
   @> Frame/model: 28
   @> Frame/model: 42
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.29s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.29s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.30s.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.30s.
   @> Delaunay tessellation of 77434 points constructed in 4.00s.
   @> Delaunay tessellation of 77434 points constructed in 4.04s.
   @> Delaunay tessellation of 77434 points constructed in 4.10s.
   @> Delaunay tessellation of 77434 points constructed in 4.16s.
   @> Surface and inner simplices filtered in 3.77s.
   @> Surface and inner simplices filtered in 3.86s.
   @> Surface and inner simplices filtered in 3.97s.
   @> Surface and inner simplices filtered in 3.84s.
   @> 6 surface cavities detected and filtered in 0.55s.
   @> 7 surface cavities detected and filtered in 0.57s.
   @> 11 surface cavities detected and filtered in 0.55s.
   @> 12 surface cavities detected and filtered in 0.55s.
   @> Channel pathfinding (graph Dijkstra) over 7 cavities completed in 0.99s.
   @> Detected 25 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 9.81s.
   @> Frame/model: 29
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Channel pathfinding (graph Dijkstra) over 11 cavities completed in 1.04s.
   @> Detected 20 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 9.91s.
   @> Frame/model: 15
   ..
   ..
   @> Frame/model: 195
   @> The atoms supplied to calcChannels contain protein atoms only.
   @> Substituted 5986 atoms with 77434 homogeneous balls of radius 1.20 A in 0.31s.
   @> Surface and inner simplices filtered in 2.70s.
   @> 9 surface cavities detected and filtered in 0.55s.
   @> Delaunay tessellation of 77434 points constructed in 3.00s.
   @> Channel pathfinding (graph Dijkstra) over 9 cavities completed in 1.28s.
   @> Detected 30 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 8.12s.
   @> Delaunay tessellation of 77434 points constructed in 3.24s.
   @> Surface and inner simplices filtered in 2.93s.
   @> Surface and inner simplices filtered in 2.69s.
   @> 9 surface cavities detected and filtered in 0.53s.
   @> 6 surface cavities detected and filtered in 0.50s.
   @> Channel pathfinding (graph Dijkstra) over 9 cavities completed in 1.43s.
   @> Detected 35 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 8.31s.
   @> Channel pathfinding (graph Dijkstra) over 6 cavities completed in 2.36s.
   @> Detected 36 channels.
   @> Saving multiple results to directory ..
   @> Channel calculation completed in 9.20s.


All the details about the predicted channels can be displayed using
:func:`.getChannelParametersMultipleFrames`.

.. ipython:: python
   :verbatim:

   getChannelParametersMultipleFrames(channels4, param_file_name='DATA_chls_dcd')

.. parsed-literal::

   @> Frame/model: 0
   @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> channel 0: 	42.62 		5.23 		0.91
   @> channel 1: 	46.26 		6.49 		0.91
   @> channel 2: 	41.15 		6.61 		0.91
   @> channel 3: 	39.01 		6.52 		0.95
   @> channel 4: 	39.8 		6.83 		0.91
   @> channel 5: 	28.15 		6.23 		0.85
   @> channel 6: 	31.99 		7.13 		0.85
   @> channel 7: 	43.34 		9.08 		0.8
   @> channel 8: 	42.4 		10.06 		0.99
   @> channel 9: 	57.35 		10.05 		0.94
   @> channel 10: 	28.1 		7.24 		0.56
   @> channel 11: 	54.05 		11.69 		0.96
   @> channel 12: 	64.98 		12.59 		0.82
   @> channel 13: 	76.95 		14.62 		0.94
   @> channel 14: 	72.35 		14.99 		0.94
   @> channel 15: 	48.2 		12.18 		0.56
   @> channel 16: 	60.7 		15.73 		0.81
   @> channel 17: 	61.95 		17.84 		0.85
   @> channel 18: 	487.25 		46.63 		0.93
   @> channel 19: 	141.6 		27.14 		0.85
   @> channel 20: 	133.72 		26.46 		0.85
   @> channel 21: 	97.65 		20.64 		0.53
   @> channel 22: 	64.82 		18.52 		0.53
   @> channel 23: 	434.36 		46.08 		0.93
   @> channel 24: 	72.47 		19.67 		0.53
   @> channel 25: 	420.49 		52.08 		0.9
   @> channel 26: 	377.6 		51.75 		0.9
   @> channel 27: 	547.44 		62.65 		0.9
   @> channel 28: 	520.97 		61.04 		0.9
   @> channel 29: 	563.87 		65.18 		0.9
   @> channel 30: 	579.84 		66.15 		0.9
   @> channel 31: 	587.81 		66.92 		0.9
   @> Frame/model: 1
   @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> channel 0: 	25.35 		5.48 		0.88
   @> channel 1: 	33.95 		6.76 		0.9
   @> channel 2: 	39.95 		7.56 		0.89
   @> channel 3: 	23.72 		6.63 		0.87
   @> channel 4: 	36.81 		10.0 		0.9
   @> channel 5: 	80.37 		13.18 		0.79
   @> channel 6: 	52.55 		11.86 		0.9
   @> channel 7: 	49.69 		13.07 		0.79
   @> channel 8: 	75.65 		14.36 		0.64
   @> channel 9: 	108.27 		22.83 		0.88
   @> channel 10: 	96.26 		21.57 		0.78
   @> channel 11: 	69.35 		18.87 		0.58
   @> channel 12: 	42.33 		14.67 		0.5
   @> channel 13: 	136.31 		26.84 		0.84
   @> channel 14: 	127.65 		26.21 		0.84
   @> channel 15: 	517.36 		47.06 		0.81
   @> channel 16: 	136.2 		28.73 		0.81
   @> channel 17: 	455.47 		48.29 		0.81
   @> channel 18: 	110.0 		27.4 		0.59
   @> channel 19: 	355.87 		48.16 		0.81
   @> channel 20: 	308.57 		48.47 		0.81
   @> channel 21: 	464.03 		55.84 		0.81
   @> channel 22: 	498.0 		58.37 		0.81
   @> channel 23: 	489.61 		58.64 		0.81
   @> channel 24: 	506.86 		60.38 		0.81
   ..
   ..
   [([5.227329141803516,
      6.487400298449215,
      6.605144066108894,
      6.516275048865574,
      6.82887082510323,
      6.225022727043587,
      7.130896203430573,
      9.0820395745718,
      10.055076137621302,
      10.046965975880378,
      7.2448800420085835,
      11.685716270931396,
      12.588204876727449,
      14.615599846655114,
      14.987739219808653,
      12.177197068381641,
      15.728096263328084,
      17.83669237347921,
      46.63095154889754,
      27.14201931744085,
      26.46087867366709,
      20.638332276673168,
      18.520993854173597,
      46.07627263551035,
      19.666847649036676,
      52.08168239809312,
      51.754432428231695,
      62.65440497610238,
      61.038934300273354,
      65.17632012481089,
      66.14572267592729,
      66.91973888604468],
     [0.9062093848777095,
      0.9062093848777095,
      0.9091246523541974,
      0.9466973308923545,
      0.9062093848777095,
      0.8450995780137116,
      0.8450995780137116,
      0.7964742507765301,
      0.9877666026933867,
      0.9440099272394715,
      0.5627200930417462,
      0.9580799362892518,
      0.8168412673605527,
      0.9368914605676404,
      0.9368914605676404,
      0.5634107045739882,
      0.8126297939221512,
      0.8458444201300528,
      0.9288219861959515,
      0.8471905678106821,
      0.8471905678106821,
      0.5336380547900532,
      0.5336380547900532,
      0.9288219861959515,
      0.5336380547900532,
      0.8980676489273869,
      0.8980676489273869,
      0.8980374597036274,
      0.8980374597036274,
      0.8980374597036274,
      0.8980374597036274,
      0.8980374597036274],
      ..
      ..
     [33.93878201451358,
      53.57138814058363,
      51.65102405007965,
      60.40988921509134,
      93.47767217117577,
      153.06628278359594,
      65.78081695746604,
      188.95516489959445,
      160.83223704342512,
      174.62962363999816,
      62.95855748022218,
      71.36168073076173,
      140.79568166555896,
      68.43219227709203,
      114.26271670969929,
      746.4890229908065,
      784.8926565168398,
      699.5382690630242,
      770.3751166765595,
      805.863155925841,
      695.3510801096598,
      767.8987104643048,
      621.7119196863381,
      797.251345235155,
      694.1672952596066,
      820.4920143280556,
      749.9826774750843,
      758.222986561308,
      788.8971708945872,
      539.6884081368177,
      794.5395102358551,
      537.7938757590493,
      473.89250290133845,
      371.30162647018585,
      383.5343643574686,
      381.43842455263564,
      324.36179374323757])]


Assigning the function's output to the ``results`` variable grants access 
to these parameters, as shown below.

.. ipython:: python
   :verbatim:

   results = getChannelParametersMultipleFrames(channels4)

.. parsed-literal::

   @> Frame/model: 0
   @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> channel 0: 	42.62 		5.23 		0.91
   @> channel 1: 	46.26 		6.49 		0.91
   @> channel 2: 	41.15 		6.61 		0.91
   @> channel 3: 	39.01 		6.52 		0.95
   @> channel 4: 	39.8 		6.83 		0.91
   @> channel 5: 	28.15 		6.23 		0.85
   @> channel 6: 	31.99 		7.13 		0.85
   @> channel 7: 	43.34 		9.08 		0.8
   @> channel 8: 	42.4 		10.06 		0.99
   @> channel 9: 	57.35 		10.05 		0.94
   @> channel 10: 	28.1 		7.24 		0.56
   @> channel 11: 	54.05 		11.69 		0.96
   @> channel 12: 	64.98 		12.59 		0.82
   @> channel 13: 	76.95 		14.62 		0.94
   @> channel 14: 	72.35 		14.99 		0.94
   @> channel 15: 	48.2 		12.18 		0.56
   @> channel 16: 	60.7 		15.73 		0.81
   @> channel 17: 	61.95 		17.84 		0.85
   @> channel 18: 	487.25 		46.63 		0.93
   @> channel 19: 	141.6 		27.14 		0.85
   @> channel 20: 	133.72 		26.46 		0.85
   @> channel 21: 	97.65 		20.64 		0.53
   @> channel 22: 	64.82 		18.52 		0.53
   @> channel 23: 	434.36 		46.08 		0.93
   @> channel 24: 	72.47 		19.67 		0.53
   @> channel 25: 	420.49 		52.08 		0.9
   @> channel 26: 	377.6 		51.75 		0.9
   @> channel 27: 	547.44 		62.65 		0.9
   @> channel 28: 	520.97 		61.04 		0.9
   @> channel 29: 	563.87 		65.18 		0.9
   @> channel 30: 	579.84 		66.15 		0.9
   @> channel 31: 	587.81 		66.92 		0.9
   @> Frame/model: 1
   @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> channel 0: 	25.35 		5.48 		0.88
   @> channel 1: 	33.95 		6.76 		0.9
   @> channel 2: 	39.95 		7.56 		0.89
   @> channel 3: 	23.72 		6.63 		0.87
   @> channel 4: 	36.81 		10.0 		0.9
   @> channel 5: 	80.37 		13.18 		0.79
   @> channel 6: 	52.55 		11.86 		0.9
   @> channel 7: 	49.69 		13.07 		0.79
   @> channel 8: 	75.65 		14.36 		0.64
   @> channel 9: 	108.27 		22.83 		0.88
   @> channel 10: 	96.26 		21.57 		0.78
   @> channel 11: 	69.35 		18.87 		0.58
   @> channel 12: 	42.33 		14.67 		0.5
   @> channel 13: 	136.31 		26.84 		0.84
   @> channel 14: 	127.65 		26.21 		0.84
   @> channel 15: 	517.36 		47.06 		0.81
   @> channel 16: 	136.2 		28.73 		0.81
   @> channel 17: 	455.47 		48.29 		0.81
   @> channel 18: 	110.0 		27.4 		0.59
   @> channel 19: 	355.87 		48.16 		0.81
   @> channel 20: 	308.57 		48.47 		0.81
   @> channel 21: 	464.03 		55.84 		0.81
   @> channel 22: 	498.0 		58.37 		0.81
   @> channel 23: 	489.61 		58.64 		0.81
   @> channel 24: 	506.86 		60.38 		0.81
   ..
   ..
   @> Frame/model: 210
   @> Channel ID: 	Volume [Å³] 	Length [Å] 	Bottleneck [Å]
   @> channel 0: 	33.94 		5.84 		0.91
   @> channel 1: 	53.57 		6.06 		0.9
   @> channel 2: 	51.65 		8.36 		0.92
   @> channel 3: 	60.41 		10.15 		0.9
   @> channel 4: 	93.48 		14.96 		0.9
   @> channel 5: 	153.07 		20.89 		0.85
   @> channel 6: 	65.78 		15.06 		0.78
   @> channel 7: 	188.96 		22.78 		0.85
   @> channel 8: 	160.83 		22.3 		0.85
   @> channel 9: 	174.63 		23.52 		0.85
   @> channel 10: 	62.96 		15.47 		0.8
   @> channel 11: 	71.36 		16.67 		0.78
   @> channel 12: 	140.8 		23.1 		0.85
   @> channel 13: 	68.43 		16.98 		0.8
   @> channel 14: 	114.26 		26.12 		0.75
   @> channel 15: 	746.49 		56.91 		0.82
   @> channel 16: 	784.89 		59.45 		0.82
   @> channel 17: 	699.54 		57.48 		0.82
   @> channel 18: 	770.38 		59.6 		0.82
   @> channel 19: 	805.86 		61.61 		0.82
   @> channel 20: 	695.35 		58.25 		0.82
   @> channel 21: 	767.9 		61.46 		0.82
   @> channel 22: 	621.71 		56.17 		0.82
   @> channel 23: 	797.25 		63.75 		0.82
   @> channel 24: 	694.17 		58.74 		0.82
   @> channel 25: 	820.49 		65.8 		0.82
   @> channel 26: 	749.98 		61.69 		0.82
   @> channel 27: 	758.22 		61.31 		0.81
   @> channel 28: 	788.9 		64.56 		0.82
   @> channel 29: 	539.69 		56.59 		0.82
   @> channel 30: 	794.54 		68.51 		0.82
   @> channel 31: 	537.79 		57.32 		0.82
   @> channel 32: 	473.89 		53.43 		0.73
   @> channel 33: 	371.3 		55.49 		0.82
   @> channel 34: 	383.53 		57.85 		0.82
   @> channel 35: 	381.44 		57.9 		0.82
   @> channel 36: 	324.36 		53.64 		0.64


Results for the first frame in the trajectory:

.. ipython:: python
   :verbatim:

   results[0]

.. parsed-literal::

   ([5.227329141803516,
     6.487400298449215,
     6.605144066108894,
     6.516275048865574,
     6.82887082510323,
     6.225022727043587,
     7.130896203430573,
     9.0820395745718,
     10.055076137621302,
     10.046965975880378,
     7.2448800420085835,
     11.685716270931396,
     12.588204876727449,
     14.615599846655114,
     14.987739219808653,
     12.177197068381641,
     15.728096263328084,
     17.83669237347921,
     46.63095154889754,
     27.14201931744085,
     26.46087867366709,
     20.638332276673168,
     18.520993854173597,
     46.07627263551035,
     19.666847649036676,
     52.08168239809312,
     51.754432428231695,
     62.65440497610238,
     61.038934300273354,
     65.17632012481089,
     66.14572267592729,
     66.91973888604468],
    [0.9062093848777095,
     0.9062093848777095,
     0.9091246523541974,
     0.9466973308923545,
     0.9062093848777095,
     0.8450995780137116,
     0.8450995780137116,
     0.7964742507765301,
     0.9877666026933867,
     0.9440099272394715,
     0.5627200930417462,
     0.9580799362892518,
     0.8168412673605527,
     0.9368914605676404,
     0.9368914605676404,
     0.5634107045739882,
     0.8126297939221512,
     0.8458444201300528,
     0.9288219861959515,
     0.8471905678106821,
     0.8471905678106821,
     0.5336380547900532,
     0.5336380547900532,
     0.9288219861959515,
     0.5336380547900532,
     0.8980676489273869,
     0.8980676489273869,
     0.8980374597036274,
     0.8980374597036274,
     0.8980374597036274,
     0.8980374597036274,
     0.8980374597036274],
    [42.618844639015094,
     46.26452069417938,
     41.145690082853065,
     39.00955277636218,
     39.7986371406858,
     28.151929526035254,
     31.988857209237942,
     43.33960313799116,
     42.39734128980631,
     57.34659622602214,
     28.095445326847372,
     54.05182624768031,
     64.97821132673613,
     76.94831436706494,
     72.34929204813164,
     48.195036012352446,
     60.699269451579816,
     61.954394708434094,
     487.24788755195357,
     141.60215966117727,
     133.71951470120246,
     97.64763069231812,
     64.81924213531595,
     434.36430272846616,
     72.47178046077198,
     420.49095799060683,
     377.5965996540025,
     547.4403640282828,
     520.965184770951,
     563.8703645100817,
     579.8383939046562,
     587.8109741174737])


To obtain information about the lengths of the channels detected in the 
first frame in the trajectory (#0):

.. ipython:: python
   :verbatim:

   results[0][0]

.. parsed-literal::

   [5.227329141803516,
    6.487400298449215,
    6.605144066108894,
    6.516275048865574,
    6.82887082510323,
    6.225022727043587,
    7.130896203430573,
    9.0820395745718,
    10.055076137621302,
    10.046965975880378,
    7.2448800420085835,
    11.685716270931396,
    12.588204876727449,
    14.615599846655114,
    14.987739219808653,
    12.177197068381641,
    15.728096263328084,
    17.83669237347921,
    46.63095154889754,
    27.14201931744085,
    26.46087867366709,
    20.638332276673168,
    18.520993854173597,
    46.07627263551035,
    19.666847649036676,
    52.08168239809312,
    51.754432428231695,
    62.65440497610238,
    61.038934300273354,
    65.17632012481089,
    66.14572267592729,
    66.91973888604468]

Bottlenecks of the channels in the first frame in the trajectory (#0):

.. ipython:: python
   :verbatim:

   results[0][1]

.. parsed-literal::

   [0.9062093848777095,
    0.9062093848777095,
    0.9091246523541974,
    0.9466973308923545,
    0.9062093848777095,
    0.8450995780137116,
    0.8450995780137116,
    0.7964742507765301,
    0.9877666026933867,
    0.9440099272394715,
    0.5627200930417462,
    0.9580799362892518,
    0.8168412673605527,
    0.9368914605676404,
    0.9368914605676404,
    0.5634107045739882,
    0.8126297939221512,
    0.8458444201300528,
    0.9288219861959515,
    0.8471905678106821,
    0.8471905678106821,
    0.5336380547900532,
    0.5336380547900532,
    0.9288219861959515,
    0.5336380547900532,
    0.8980676489273869,
    0.8980676489273869,
    0.8980374597036274,
    0.8980374597036274,
    0.8980374597036274,
    0.8980374597036274,
    0.8980374597036274]

Volume of the channels detected in the first frame in the trajectory:

.. ipython:: python
   :verbatim:

   results[0][2]

.. parsed-literal::

   [42.618844639015094,
    46.26452069417938,
    41.145690082853065,
    39.00955277636218,
    39.7986371406858,
    28.151929526035254,
    31.988857209237942,
    43.33960313799116,
    42.39734128980631,
    57.34659622602214,
    28.095445326847372,
    54.05182624768031,
    64.97821132673613,
    76.94831436706494,
    72.34929204813164,
    48.195036012352446,
    60.699269451579816,
    61.954394708434094,
    487.24788755195357,
    141.60215966117727,
    133.71951470120246,
    97.64763069231812,
    64.81924213531595,
    434.36430272846616,
    72.47178046077198,
    420.49095799060683,
    377.5965996540025,
    547.4403640282828,
    520.965184770951,
    563.8703645100817,
    579.8383939046562,
    587.8109741174737]

Once we have access to ``results``, we can display bottleneck data for all
the channels in the following way:

.. ipython:: python
   :verbatim:

   import matplotlib.pylab as plt
   all_Bottleneck = []
   for nr_i,i in enumerate(results):
      all_Bottleneck.extend(results[nr_i][1])

   plt.hist(all_Bottleneck)
   plt.show()

.. parsed-literal::

   (array([  15.,   20.,   68.,  105.,  186.,  475., 1295., 3483.,  953.,
             43.]),
    array([0.10125895, 0.20444863, 0.30763831, 0.410828  , 0.51401768,
           0.61720736, 0.72039704, 0.82358673, 0.92677641, 1.02996609,
           1.13315577]),
    <BarContainer object of 10 artists>)


.. figure:: images/cavitracer_figure16.jpg
   :scale: 50 %


To have access to a particular frame, we need to use :meth:`.getFrame`.
Below, we will select third frame from the simulation (counting from 0):

.. ipython:: python
   :verbatim:

   frame3 = dcd.getFrame(2)
   frame3

.. parsed-literal::

   <Frame: 2 from caseStudy2 (selected 5986 of 56457 atoms)>


To obtain infromation about residues for a particular frame, the following
operations should be performed:

.. ipython:: python
   :verbatim:

   protein_frame3 = protein.copy()
   protein_frame3.setCoords(frame3.getCoords())

.. ipython:: python
   :verbatim:

   getChannelResidueNames(protein_frame3, channels4[2], residues_file_name='DCD_fr3_res')

.. parsed-literal::

   @> 6076 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6101 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6056 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6076 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6081 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6046 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6056 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6071 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6121 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6181 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6171 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6211 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6366 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6366 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6346 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6406 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6396 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6406 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6431 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6436 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6466 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6461 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6466 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> 6466 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6521 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6596 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6536 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 6536 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> Channel residues were saved to: DCD_fr3_res_Residues_All_channels.txt

   ['channel0: MET310, ALA444, PHE449, PRO450, LEU452, MET453, ILE456',
    'channel1: THR212, ASP214, ARG217, MET221, ARG357, MET403, PRO404, GLY407, TYR408, ASP411, TYR422',
    'channel2: ARG19, VAL210, TYR211, THR212, ASP213, GLU216, ARG217',
    'channel3: TYR293, ARG357, ILE405, TYR408, LEU409, PHE469, LEU470, PRO473',
    'channel4: LEU285, LEU288, LEU289, ILE294, SER420, VAL421, ALA423, ILE424, VAL427',
    'channel5: TYR293, ILE405, TYR408, LEU409, PHE469, LEU470, ARG471, SER472, PRO473',
    'channel6: MET310, ALA444, PHE449, PRO450, LEU452, MET453, ILE456',
    'channel7: LEU288, LEU289, LYS290, ASP291, ILE294, LEU295, ALA298, ILE424',
    'channel8: LEU151, ASN154, LYS281, GLY282, THR283, PRO284, LEU288, GLY407, VAL410, VAL417, TYR418, SER420, VAL421',
    'channel9: SER300, PHE303, ALA304, GLY364, MET365, VAL368, ALA394, ILE395, MET397, VAL398, ILE459, ASP460, PHE463',
    'channel10: SER300, PHE303, ALA304, GLY364, MET365, VAL368, ALA394, ILE395, MET397, VAL398, ILE459, ASP460, LEU462, PHE463',
    'channel11: MET204, ARG217, GLY218, MET221, GLY222, ALA224, LEU225, GLY349, ALA352, HSP353, LYS354, GLY356, ARG357, LEU359, CYS360, LEU363, ASP399, SER400, SER401, MET403, PRO404',
    'channel12: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, GLU312, LEU315, ILE317, GLN329, VAL332, ALA333, PHE334, PRO336, ALA337, TYR341, ILE381, ILE385, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel13: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, GLU312, LEU315, ILE317, GLN329, ALA333, PHE334, PRO336, ALA337, TYR341, ILE381, TYR382, ILE385, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel14: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, GLU312, LEU315, ILE317, VAL332, ALA333, PHE334, PRO336, ALA337, TYR341, ILE381, ILE385, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel15: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, GLU312, LEU315, ILE317, ARG326, GLN329, ALA333, PHE334, PRO336, ALA337, TYR341, ILE381, ILE385, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel16: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, GLU312, LEU315, ILE317, ARG326, GLN329, ALA333, PHE334, PRO336, ALA337, TYR341, ILE381, TYR382, ILE385, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel17: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, GLU312, LEU315, PRO316, ILE317, THR322, ARG326, GLN329, ALA333, PHE334, PRO336, ALA337, TYR341, ILE381, ILE385, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel18: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, GLU312, LEU315, PRO316, ILE317, MET319, GLU321, THR322, ARG326, GLN329, ALA333, PHE334, PRO336, ALA337, TYR341, ILE381, ILE385, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel19: ASP33, ASN34, LEU37, THR38, VAL40, VAL41, PRO42, ILE44, PHE135, LYS138, GLN142, ARG189, SER196, SER200, LEU225, LEU228, PRO236, SER240, ALA297, ALA298, ILE301, GLU312, ILE317, TRP318, LYS327, TRP328, LEU330, GLY331, ALA333, PHE334, TYR341, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel20: ASP33, ASN34, LEU37, THR38, VAL40, VAL41, PRO42, ILE44, PHE135, LYS138, GLN142, ARG189, SER196, SER200, LEU225, LEU228, PRO236, SER240, TYR243, GLU244, ALA297, ALA298, ILE301, GLU312, ILE317, TRP318, LYS327, TRP328, LEU330, GLY331, ALA333, PHE334, TYR341, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel21: ASP33, ASN34, LEU37, THR38, VAL40, VAL41, PRO42, ILE44, PHE135, LYS138, GLN142, ARG189, SER196, SER200, LEU225, LEU228, PRO236, SER240, TYR243, GLU244, ALA297, ALA298, ILE301, GLU312, ILE317, TRP318, LYS327, LEU330, GLY331, ALA333, PHE334, TYR341, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel22: ASP33, ASN34, LEU37, THR38, VAL40, VAL41, PRO42, ILE44, PRO45, PHE135, LYS138, GLN142, ARG189, SER196, SER200, LEU225, LEU228, PRO236, SER240, TYR243, ALA297, ALA298, ILE301, GLU312, ILE317, TRP318, LYS327, LEU330, GLY331, ALA333, PHE334, TYR341, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel23: ASP33, ASN34, LEU37, THR38, VAL40, VAL41, PRO42, ILE44, PRO45, PHE135, LYS138, GLN142, ARG189, SER196, SER200, LEU225, LEU228, PRO236, SER240, TYR243, ALA297, ALA298, ILE301, GLU312, ILE317, TRP318, MET323, LYS327, LEU330, GLY331, ALA333, PHE334, TYR341, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433',
    'channel24: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, LEU311, GLU312, LEU315, ALA333, PHE334, PRO336, ALA337, TYR341, SER371, ILE372, LEU373, ILE375, PRO376, LEU384, ILE385, PRO387, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433, PRO450, TRP451, MET453, THR454',
    'channel25: ASP33, ASN34, LEU37, THR38, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, PHE303, ALA304, MET306, GLY307, ILE308, LEU311, GLU312, LEU315, ALA333, PHE334, PRO336, ALA337, TYR341, VAL368, ILE372, LEU384, ILE385, PRO387, ASN388, VAL391, ILE395, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433, ILE456, ILE459, ASP460',
    'channel26: ASP33, ASN34, LEU37, THR38, LYS122, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, LEU311, GLU312, LEU315, ALA333, PHE334, PRO336, ALA337, TYR341, SER371, ILE372, LEU373, ILE375, PRO376, LYS379, LEU384, ILE385, PRO387, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433, PRO450, MET453, THR454',
    'channel27: ASP33, ASN34, LEU37, THR38, LYS122, LYS138, GLN142, SER196, SER200, LEU225, LEU228, VAL232, ALA297, ALA298, ILE301, LEU311, GLU312, LEU315, ALA333, PHE334, PRO336, ALA337, TYR341, SER371, ILE372, LEU373, ILE375, PRO376, LEU384, ILE385, PRO387, ASN388, ASP399, MET402, MET403, MET406, ALA425, ASP426, PHE429, TYR433, PRO450, MET453, THR454']


To display the results for a particular frame, we need to create a model
using :func:`.getVmdModel` function and provide a path to VMD_.

.. ipython:: python
   :verbatim:

   vmd_path = '/usr/local/bin/vmd'
   model3_traj = getVmdModel(vmd_path, frame3)

.. parsed-literal::

   @> Model created successfully.

Once the model is created, the channels can be displayed together with the
protein structure using :func:`.showChannels`.

.. ipython:: python
   :verbatim:

   showChannels(channels4[2], model=model3_traj)

.. figure:: images/cavitracer_figure15.jpg
   :scale: 50 %

To select channels with particular localization, we can further use
:func:`.selectChannelBySelection` function and specify which region will be
analyzed, for example, using residues, as shown below. In the example, we
are selecting channels that are generated from MD trajectory using
``"chls_dcd*_chl*.pqr"`` pattern. List with ``pqr`` files
called ``pqr_files_channels`` is further use by
:func:`.selectChannelBySelection`.

.. ipython:: python
   :verbatim:
   
   from pathlib import Path
   pqr_files_channels = [i.name for i in Path(".").glob("chls_dcd*_chl*.pqr") if i.is_file()]
   pqr_files_channels

.. parsed-literal::

   ['chls_dcd118_chl19.pqr',
    'chls_dcd185_chl31.pqr',
    'chls_dcd13_chl10.pqr',
    'chls_dcd186_chl16.pqr',
    'chls_dcd47_chl16.pqr',
    'chls_dcd3_chl27.pqr',
    'chls_dcd184_chl18.pqr',
    'chls_dcd52_chl21.pqr',
    'chls_dcd174_chl5.pqr',
    'chls_dcd16_chl17.pqr',
    'chls_dcd163_chl24.pqr',
    'chls_dcd173_chl11.pqr',
    'chls_dcd204_chl6.pqr',
    'chls_dcd63_chl39.pqr',
    'chls_dcd183_chl4.pqr',
    ..
    'chls_dcd90_chl11.pqr',
    'chls_dcd124_chl18.pqr',
    'chls_dcd145_chl19.pqr',
    'chls_dcd22_chl7.pqr',
    'chls_dcd125_chl14.pqr',
    'chls_dcd46_chl12.pqr',
    'chls_dcd29_chl35.pqr',
    'chls_dcd12_chl9.pqr',
    ...]


.. ipython:: python
   :verbatim:

   selectChannelBySelection(atoms, pqr_files=pqr_files_channels, 
			residue_sele='resid 135 37 and backbone', 
                        folder_name="Selected_channel1", 
			distA=4.0)

.. parsed-literal::

   @> 390 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 455 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 105 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 325 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 115 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 455 atoms and 1 coordinate sets were parsed in 0.01s.
   @> Filtered files are now in: Selected_channel1
   @> 410 atoms and 1 coordinate sets were parsed in 0.01s.
   @> Filtered files are now in: Selected_channel1
   @> 565 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 70 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 380 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 215 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 180 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 85 atoms and 1 coordinate sets were parsed in 0.00s.
   ..
   ..
   @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 210 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 85 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 475 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 140 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 65 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 400 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 440 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: Selected_channel1
   @> 155 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 490 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 75 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 450 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: Selected_channel1
   @> 390 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 100 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 60 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 430 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 335 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: Selected_channel1
   @> 130 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 170 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 165 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 130 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 435 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Filtered files are now in: Selected_channel1
   @> 545 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 385 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 505 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 570 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 70 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Selected files: 
   @> chls_dcd3_chl27.pqr chls_dcd184_chl18.pqr chls_dcd75_chl25.pqr chls_dcd85_chl16.pqr chls_dcd140_chl12.pqr chls_dcd174_chl19.pqr chls_dcd204_chl35.pqr chls_dcd66_chl17.pqr chls_dcd97_chl19.pqr chls_dcd3_chl24.pqr chls_dcd44_chl27.pqr chls_dcd205_chl19.pqr chls_dcd154_chl18.pqr chls_dcd134_chl13.pqr chls_dcd135_chl26.pqr chls_dcd116_chl29.pqr chls_dcd181_chl20.pqr chls_dcd14_chl19.pqr chls_dcd69_chl21.pqr chls_dcd87_chl13.pqr chls_dcd120_chl5.pqr chls_dcd209_chl14.pqr chls_dcd101_chl21.pqr chls_dcd95_chl10.pqr chls_dcd142_chl12.pqr chls_dcd143_chl14.pqr chls_dcd113_chl21.pqr chls_dcd10_chl37.pqr chls_dcd135_chl18.pqr chls_dcd12_chl17.pqr chls_dcd17_chl24.pqr chls_dcd3_chl20.pqr chls_dcd124_chl7.pqr chls_dcd72_chl16.pqr chls_dcd195_chl21.pqr chls_dcd54_chl28.pqr chls_dcd14_chl17.pqr chls_dcd136_chl29.pqr chls_dcd153_chl11.pqr chls_dcd31_chl21.pqr chls_dcd83_chl18.pqr chls_dcd77_chl23.pqr chls_dcd21_chl16.pqr chls_dcd74_chl29.pqr chls_dcd80_chl20.pqr chls_dcd93_chl20.pqr chls_dcd209_chl23.pqr chls_dcd74_chl31.pqr chls_dcd5_chl19.pqr chls_dcd118_chl13.pqr chls_dcd157_chl20.pqr chls_dcd65_chl33.pqr chls_dcd8_chl29.pqr chls_dcd3_chl19.pqr chls_dcd102_chl19.pqr chls_dcd21_chl15.pqr chls_dcd173_chl9.pqr chls_dcd9_chl23.pqr chls_dcd8_chl30.pqr chls_dcd82_chl17.pqr chls_dcd131_chl25.pqr chls_dcd101_chl29.pqr chls_dcd29_chl14.pqr chls_dcd103_chl25.pqr chls_dcd127_chl21.pqr chls_dcd19_chl27.pqr chls_dcd116_chl22.pqr chls_dcd122_chl17.pqr chls_dcd90_chl19.pqr chls_dcd145_chl19.pqr chls_dcd125_chl14.pqr chls_dcd151_chl31.pqr chls_dcd111_chl18.pqr chls_dcd191_chl19.pqr chls_dcd88_chl26.pqr chls_dcd153_chl9.pqr chls_dcd97_chl20.pqr chls_dcd196_chl18.pqr chls_dcd85_chl15.pqr chls_dcd106_chl19.pqr chls_dcd133_chl15.pqr chls_dcd126_chl24.pqr chls_dcd17_chl25.pqr chls_dcd16_chl20.pqr chls_dcd173_chl17.pqr chls_dcd129_chl17.pqr chls_dcd206_chl20.pqr chls_dcd137_chl23.pqr chls_dcd58_chl23.pqr chls_dcd135_chl24.pqr chls_dcd104_chl23.pqr chls_dcd24_chl9.pqr chls_dcd0_chl26.pqr chls_dcd92_chl21.pqr chls_dcd173_chl13.pqr chls_dcd100_chl20.pqr chls_dcd22_chl19.pqr chls_dcd4_chl26.pqr chls_dcd80_chl16.pqr chls_dcd19_chl23.pqr chls_dcd22_chl20.pqr chls_dcd0_chl25.pqr chls_dcd16_chl25.pqr chls_dcd19_chl25.pqr chls_dcd11_chl22.pqr chls_dcd165_chl17.pqr chls_dcd1_chl20.pqr chls_dcd146_chl13.pqr chls_dcd95_chl14.pqr chls_dcd83_chl17.pqr chls_dcd101_chl20.pqr chls_dcd14_chl18.pqr chls_dcd74_chl23.pqr chls_dcd150_chl26.pqr chls_dcd22_chl16.pqr chls_dcd150_chl27.pqr chls_dcd155_chl17.pqr chls_dcd114_chl17.pqr chls_dcd109_chl23.pqr chls_dcd59_chl22.pqr chls_dcd146_chl14.pqr chls_dcd209_chl12.pqr chls_dcd148_chl14.pqr chls_dcd94_chl25.pqr chls_dcd48_chl29.pqr chls_dcd149_chl27.pqr chls_dcd184_chl15.pqr chls_dcd194_chl27.pqr chls_dcd19_chl21.pqr chls_dcd46_chl29.pqr chls_dcd45_chl12.pqr chls_dcd44_chl26.pqr chls_dcd105_chl23.pqr chls_dcd118_chl15.pqr chls_dcd76_chl20.pqr chls_dcd60_chl26.pqr chls_dcd186_chl20.pqr chls_dcd97_chl34.pqr chls_dcd66_chl12.pqr chls_dcd175_chl19.pqr chls_dcd107_chl36.pqr chls_dcd143_chl16.pqr chls_dcd113_chl9.pqr chls_dcd197_chl30.pqr chls_dcd6_chl30.pqr chls_dcd100_chl21.pqr chls_dcd172_chl16.pqr chls_dcd148_chl18.pqr chls_dcd210_chl20.pqr chls_dcd206_chl12.pqr chls_dcd3_chl29.pqr chls_dcd68_chl14.pqr chls_dcd210_chl17.pqr chls_dcd130_chl14.pqr chls_dcd45_chl15.pqr chls_dcd92_chl16.pqr chls_dcd171_chl13.pqr chls_dcd132_chl22.pqr chls_dcd173_chl27.pqr chls_dcd67_chl15.pqr chls_dcd60_chl27.pqr chls_dcd155_chl11.pqr chls_dcd209_chl11.pqr chls_dcd201_chl14.pqr chls_dcd165_chl19.pqr chls_dcd131_chl24.pqr chls_dcd59_chl24.pqr chls_dcd157_chl26.pqr chls_dcd93_chl19.pqr chls_dcd59_chl26.pqr chls_dcd12_chl16.pqr chls_dcd21_chl18.pqr chls_dcd125_chl9.pqr chls_dcd201_chl20.pqr chls_dcd68_chl15.pqr chls_dcd5_chl20.pqr chls_dcd33_chl41.pqr chls_dcd50_chl26.pqr chls_dcd203_chl19.pqr chls_dcd144_chl16.pqr chls_dcd121_chl15.pqr chls_dcd9_chl22.pqr chls_dcd132_chl24.pqr chls_dcd157_chl27.pqr chls_dcd148_chl23.pqr chls_dcd122_chl18.pqr chls_dcd167_chl13.pqr chls_dcd119_chl22.pqr chls_dcd20_chl23.pqr chls_dcd3_chl25.pqr chls_dcd59_chl23.pqr chls_dcd22_chl18.pqr chls_dcd65_chl29.pqr chls_dcd11_chl25.pqr chls_dcd68_chl21.pqr chls_dcd120_chl10.pqr chls_dcd87_chl17.pqr chls_dcd22_chl21.pqr chls_dcd151_chl26.pqr chls_dcd1_chl19.pqr chls_dcd149_chl28.pqr chls_dcd70_chl18.pqr chls_dcd170_chl25.pqr chls_dcd30_chl17.pqr chls_dcd206_chl24.pqr chls_dcd27_chl29.pqr chls_dcd180_chl23.pqr chls_dcd60_chl24.pqr chls_dcd19_chl26.pqr chls_dcd188_chl11.pqr chls_dcd182_chl15.pqr chls_dcd13_chl29.pqr chls_dcd176_chl17.pqr chls_dcd107_chl41.pqr chls_dcd72_chl18.pqr chls_dcd55_chl31.pqr chls_dcd204_chl32.pqr chls_dcd114_chl13.pqr chls_dcd97_chl22.pqr chls_dcd119_chl21.pqr chls_dcd182_chl17.pqr chls_dcd149_chl24.pqr chls_dcd20_chl22.pqr chls_dcd191_chl20.pqr chls_dcd152_chl17.pqr chls_dcd113_chl22.pqr chls_dcd181_chl16.pqr chls_dcd59_chl29.pqr chls_dcd81_chl24.pqr chls_dcd147_chl17.pqr chls_dcd115_chl23.pqr chls_dcd33_chl39.pqr chls_dcd180_chl27.pqr chls_dcd45_chl18.pqr chls_dcd19_chl24.pqr chls_dcd197_chl29.pqr chls_dcd158_chl17.pqr chls_dcd157_chl25.pqr chls_dcd44_chl30.pqr chls_dcd110_chl17.pqr chls_dcd175_chl22.pqr chls_dcd4_chl24.pqr chls_dcd146_chl16.pqr chls_dcd157_chl32.pqr chls_dcd179_chl31.pqr chls_dcd153_chl8.pqr chls_dcd112_chl13.pqr chls_dcd95_chl8.pqr chls_dcd180_chl25.pqr chls_dcd10_chl35.pqr chls_dcd121_chl21.pqr chls_dcd123_chl19.pqr chls_dcd177_chl22.pqr chls_dcd53_chl21.pqr chls_dcd57_chl25.pqr chls_dcd92_chl20.pqr chls_dcd188_chl15.pqr chls_dcd157_chl28.pqr chls_dcd120_chl6.pqr chls_dcd96_chl13.pqr chls_dcd71_chl19.pqr chls_dcd127_chl13.pqr chls_dcd198_chl5.pqr chls_dcd3_chl28.pqr chls_dcd206_chl11.pqr chls_dcd113_chl11.pqr chls_dcd150_chl23.pqr chls_dcd90_chl21.pqr chls_dcd87_chl18.pqr chls_dcd8_chl27.pqr chls_dcd196_chl25.pqr chls_dcd16_chl23.pqr chls_dcd10_chl32.pqr chls_dcd178_chl11.pqr chls_dcd189_chl4.pqr chls_dcd107_chl32.pqr chls_dcd145_chl27.pqr chls_dcd194_chl26.pqr chls_dcd207_chl6.pqr chls_dcd91_chl13.pqr chls_dcd8_chl28.pqr chls_dcd134_chl15.pqr chls_dcd187_chl10.pqr chls_dcd126_chl25.pqr chls_dcd191_chl21.pqr chls_dcd202_chl22.pqr chls_dcd16_chl22.pqr chls_dcd31_chl20.pqr chls_dcd194_chl29.pqr chls_dcd181_chl17.pqr chls_dcd202_chl20.pqr chls_dcd107_chl40.pqr chls_dcd60_chl20.pqr chls_dcd16_chl21.pqr chls_dcd164_chl19.pqr chls_dcd21_chl17.pqr chls_dcd57_chl36.pqr chls_dcd187_chl8.pqr chls_dcd123_chl14.pqr chls_dcd173_chl18.pqr chls_dcd164_chl28.pqr chls_dcd15_chl27.pqr chls_dcd75_chl15.pqr chls_dcd98_chl20.pqr chls_dcd154_chl21.pqr chls_dcd60_chl28.pqr chls_dcd18_chl22.pqr chls_dcd128_chl16.pqr chls_dcd99_chl15.pqr chls_dcd111_chl16.pqr chls_dcd141_chl19.pqr chls_dcd10_chl34.pqr chls_dcd205_chl18.pqr chls_dcd141_chl26.pqr chls_dcd62_chl22.pqr chls_dcd107_chl37.pqr chls_dcd126_chl28.pqr chls_dcd81_chl23.pqr chls_dcd127_chl17.pqr chls_dcd195_chl9.pqr chls_dcd173_chl20.pqr chls_dcd200_chl19.pqr chls_dcd207_chl9.pqr chls_dcd111_chl9.pqr chls_dcd210_chl24.pqr chls_dcd198_chl9.pqr chls_dcd165_chl18.pqr chls_dcd76_chl25.pqr chls_dcd198_chl17.pqr chls_dcd6_chl29.pqr chls_dcd17_chl21.pqr chls_dcd202_chl18.pqr chls_dcd14_chl16.pqr chls_dcd28_chl15.pqr chls_dcd80_chl21.pqr chls_dcd197_chl33.pqr chls_dcd110_chl16.pqr chls_dcd170_chl26.pqr chls_dcd108_chl19.pqr chls_dcd175_chl9.pqr chls_dcd168_chl18.pqr chls_dcd39_chl27.pqr chls_dcd136_chl24.pqr chls_dcd180_chl21.pqr chls_dcd10_chl36.pqr chls_dcd130_chl19.pqr chls_dcd3_chl18.pqr chls_dcd144_chl24.pqr chls_dcd152_chl18.pqr chls_dcd196_chl22.pqr chls_dcd130_chl17.pqr chls_dcd200_chl17.pqr chls_dcd167_chl16.pqr chls_dcd116_chl28.pqr chls_dcd181_chl21.pqr chls_dcd187_chl11.pqr chls_dcd193_chl21.pqr chls_dcd180_chl20.pqr chls_dcd124_chl16.pqr chls_dcd13_chl27.pqr chls_dcd3_chl21.pqr chls_dcd10_chl33.pqr chls_dcd156_chl23.pqr chls_dcd84_chl13.pqr chls_dcd117_chl15.pqr chls_dcd131_chl26.pqr chls_dcd148_chl19.pqr chls_dcd3_chl26.pqr chls_dcd71_chl14.pqr chls_dcd173_chl22.pqr chls_dcd157_chl18.pqr chls_dcd149_chl26.pqr chls_dcd55_chl32.pqr chls_dcd184_chl13.pqr chls_dcd104_chl26.pqr chls_dcd116_chl30.pqr chls_dcd33_chl40.pqr chls_dcd93_chl16.pqr chls_dcd143_chl13.pqr chls_dcd119_chl19.pqr chls_dcd9_chl18.pqr chls_dcd18_chl23.pqr chls_dcd88_chl22.pqr chls_dcd57_chl22.pqr chls_dcd125_chl15.pqr chls_dcd71_chl20.pqr chls_dcd156_chl17.pqr chls_dcd185_chl28.pqr chls_dcd146_chl12.pqr chls_dcd57_chl38.pqr chls_dcd82_chl20.pqr chls_dcd192_chl16.pqr chls_dcd3_chl31.pqr chls_dcd162_chl19.pqr chls_dcd29_chl13.pqr chls_dcd3_chl30.pqr chls_dcd102_chl24.pqr chls_dcd64_chl25.pqr chls_dcd20_chl21.pqr chls_dcd106_chl18.pqr chls_dcd127_chl27.pqr chls_dcd167_chl22.pqr chls_dcd9_chl21.pqr chls_dcd12_chl18.pqr chls_dcd173_chl24.pqr chls_dcd162_chl17.pqr chls_dcd71_chl12.pqr chls_dcd3_chl23.pqr chls_dcd15_chl26.pqr chls_dcd4_chl25.pqr chls_dcd101_chl24.pqr chls_dcd112_chl15.pqr chls_dcd124_chl21.pqr chls_dcd114_chl15.pqr chls_dcd11_chl20.pqr chls_dcd207_chl11.pqr chls_dcd132_chl17.pqr chls_dcd117_chl14.pqr chls_dcd91_chl11.pqr chls_dcd107_chl35.pqr chls_dcd158_chl12.pqr chls_dcd109_chl29.pqr chls_dcd19_chl22.pqr chls_dcd145_chl26.pqr chls_dcd178_chl16.pqr chls_dcd109_chl19.pqr chls_dcd15_chl25.pqr chls_dcd188_chl10.pqr chls_dcd5_chl18.pqr chls_dcd147_chl19.pqr chls_dcd71_chl15.pqr chls_dcd88_chl32.pqr chls_dcd70_chl21.pqr chls_dcd94_chl23.pqr chls_dcd105_chl22.pqr chls_dcd4_chl23.pqr chls_dcd129_chl18.pqr chls_dcd128_chl20.pqr chls_dcd89_chl15.pqr chls_dcd175_chl20.pqr chls_dcd127_chl20.pqr chls_dcd101_chl19.pqr chls_dcd86_chl10.pqr chls_dcd133_chl12.pqr chls_dcd83_chl13.pqr chls_dcd193_chl15.pqr chls_dcd202_chl21.pqr chls_dcd73_chl14.pqr chls_dcd109_chl28.pqr chls_dcd66_chl13.pqr chls_dcd171_chl19.pqr chls_dcd171_chl11.pqr chls_dcd147_chl18.pqr chls_dcd172_chl11.pqr chls_dcd55_chl28.pqr chls_dcd6_chl28.pqr chls_dcd15_chl10.pqr chls_dcd161_chl26.pqr chls_dcd112_chl18.pqr chls_dcd156_chl16.pqr chls_dcd8_chl31.pqr chls_dcd192_chl20.pqr chls_dcd197_chl35.pqr chls_dcd133_chl20.pqr chls_dcd139_chl21.pqr chls_dcd118_chl20.pqr chls_dcd96_chl15.pqr
   @> If newly created files are empty please check whether the parameter names are: PDB_id+_Parameters_All_channels.txt

Filtered pqr files will be stored in ``"Selected_channel1"``. Once the files
are filtered, the :func:`.calcChannelSurfaceOverlaps` function can be used to
display the overlapping surface that is shared by filtered channels, which
will be saved when using ``output_file_name`` option.

.. ipython:: python
   :verbatim:

   calcChannelSurfaceOverlaps(pqr_files="./Selected_channel1", 
			output_file_name='overlapping_surf_traj.pdb',
			max_proc=4)

.. parsed-literal::

   @> Number of PQR files: 467
   @> Resolution: 0.5
   @> max_proc: 4
   @> Calculating overlaps using 4 processes.
   @> 450 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 455 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 395 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 455 atoms and 1 coordinate sets were parsed in 0.01s.
   @> 400 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 285 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 410 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 385 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 480 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 475 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 295 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 450 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 495 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 415 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 455 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 445 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 500 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 185 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 390 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 425 atoms and 1 coordinate sets were parsed in 0.00s.
   ..
   ..
   @> 505 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 435 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 405 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 600 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 345 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 450 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 345 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 370 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 490 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 440 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 450 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 335 atoms and 1 coordinate sets were parsed in 0.00s.
   @> 435 atoms and 1 coordinate sets were parsed in 0.00s.
   @> Overlap written to: overlapping_surf_traj.pdb
   @> Number of occupied overlap voxels: 68347

   'overlapping_surf_traj.pdb'


.. figure:: images/cavitracer_figure17.jpg
   :scale: 50 %

.. figure:: images/cavitracer_figure18.jpg
   :scale: 50 %


.. _Trajectory tutorial: http://www.bahargroup.org/prody/tutorials/trajectory_analysis/

