.. _wbfinder_tutorial:
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

   PDBtraj_file = "5kqm_all_sci.pdb"
   coords_traj = parsePDB(PDBtraj_file)
   trajectory = parseDCD("NAMD_D2_sample.dcd")


The analysis od water bridges can be performed on selected frames by using 
*start_frame* or *stop_frame*. 


.. ipython:: python

   wb_traj = calcWaterBridgesTrajectory(coords_traj, trajectory, start_frame=5, stop_frame=15, output='info')


Because of the number of data there results will not be displayed. We can have 
an access to the raw data by using *output='info'*.


.. ipython:: python

   wb_traj


Save the results
-------------------------------------------------------------------------------

The results can be saved using :func:`.saveWaterBridges` in two formats. Txt 
file will contain all the results for analysis and can be visualized in text 
editor, and wb file will restore data for further analysis. It can be uploaded 
using :func:`.parseWaterBridges` as shown below.


.. ipython:: python

   waterBridges_save = calcWaterBridgesTrajectory(coords_traj, trajectory, stop_frame=15)
   saveWaterBridges(waterBridges_save,'wb_saved.txt')
   saveWaterBridges(waterBridges_save,'wb_saved.wb')

To upload wb file use :func:`.parseWaterBridges` and protein coordinates 
as follows:


.. ipython:: python

   uploaded_results = parseWaterBridges('wb_saved.wb',coords_traj)


Uploaded results are of type atomic (.wb file), therefore it can be used for 
analysis later. The atomic output can be transformed to 
detailed information using :func:`.getWaterBridgesInfoOutput`.


Analysis of the results
-------------------------------------------------------------------------------

To perform analysis we can not use *output='info'*, therefore we will run 
the calculations again. This time we will run first 15 frames of the simulation.


.. ipython:: python

   waterBridges = calcWaterBridgesTrajectory(coords_traj, trajectory, stop_frame=15)


Information about residues contributiong to water bridges
-------------------------------------------------------------------------------

Analysis of the data can be performed using :func:`.calcWaterBridgesStatistics`.
The analysis presented below gave information about pairs of residues involved 
in water bridges, how often they occur, and the average distance between them. 
Standard deviation provides information on how the distance was changing during 
the simulation.Additionally, the analysis can be saved by using *filename* option.


.. ipython:: python
   
   analysisAtomic = calcWaterBridgesStatistics(waterBridges, trajectory, filename='data.txt')

   for item in analysisAtomic.items():
      print(item)


To have an access to the data we can use :func:`.getWaterBridgeStatInfo`.


.. ipython:: python
   
   getWaterBridgeStatInfo(analysisAtomic, coords_traj)


To obtain maps of interactions for protein structure, we can use 
:func:`.showWaterBridgeMatrix` which is equipted in three paramaters: 
*'percentage'* (how often two residues were forming water bridges), 
*'distAvg'* (how close there were), and *'distStd'* (how stable that 
water bridge was).


.. ipython:: python
   
   showWaterBridgeMatrix(analysisAtomic, 'percentage')


.. ipython:: python
   
   showWaterBridgeMatrix(analysisAtomic, 'distAvg')


.. ipython:: python
   
   showWaterBridgeMatrix(analysisAtomic, 'distStd')


Raw data of the matrices can be obatined with :func:`.calcWaterBridgeMatrix`. 
The type of the matrix can be selected among: *'percentage'*, *'distAvg'*, *'distStd'*.


.. ipython:: python

    M1 = calcWaterBridgeMatrix(analysisAtomic, 'percentage')
    M2 = calcWaterBridgeMatrix(analysisAtomic, 'distAvg')
    M3 = calcWaterBridgeMatrix(analysisAtomic, 'distStd')


.. ipython:: python

   M1


.. ipython:: python

   M2


.. ipython:: python

   M3


Statistical analysis for water bridges
-------------------------------------------------------------------------------

To visualize the results in a more accessible way, we can use 
:func:`.calcWaterBridgeMatrix` function which will show how often each residue 
were contributing to the water bridges in the trajectory.


.. ipython:: python

    calcBridgingResiduesHistogram(waterBridges)

*clip* option can be used to include different number of results on the histogram.


.. ipython:: python
    
    calcBridgingResiduesHistogram(waterBridges, clip=25)


If we are interested in one particular residue, we can also use
:func:`.calcWaterBridgesDistribution` to find their partners in water bridges. 
Below we can see results for arginine 147 or aspartic acid 92 from chain P.


.. ipython:: python

    calcWaterBridgesDistribution(waterBridges, 'ARG147P')


.. ipython:: python

    calcWaterBridgesDistribution(waterBridges, 'ASP92P')


Once we select a pair of residues which are supported by interactions with water 
molecules we can use :func:`.calcWaterBridgesDistribution` to obtain histograms 
with results such as distances between them *(metric='distance')*, the number of 
water molecules which were involved *(metric='waters')*, and information about 
residue part which was involved in water bridges, i.e. backbone or side chain 
*(metric='location')*. 


.. ipython:: python

   calcWaterBridgesDistribution(waterBridges,  'ASP92P', 'ARG18P', trajectory=trajectory, metric='distance')


.. ipython:: python

   calcWaterBridgesDistribution(waterBridges, 'ARG147P', 'GLN122P', metric='waters') 


.. ipython:: python

   calcWaterBridgesDistribution(waterBridges, 'ARG147P', 'GLN122P', trajectory=trajectory, metric='location')



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

   savePDBWaterBridges(waterBridges[0], coords_traj, PDBtraj_file[:-4]+'_frame0.pdb')

   savePDBWaterBridgesTrajectory(waterBridges, coords_traj, filename=PDBtraj_file[:-4]+'_multi.pdb', trajectory=trajectory)


Results saved in PDB file can be displayed as follows:


.. figure:: images/Fig2.tga
   :scale: 50 %


