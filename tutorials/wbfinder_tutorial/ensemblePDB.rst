.. _wbfinder_tutorial:

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

   ens = 'pebp1_50frames.pdb'
   coords_ens = parsePDB(ens)
   bridgeFrames_ens = calcWaterBridgesTrajectory(coords_ens, coords_ens)


Analysis of the results is similar to the one presented in trajectory analysis. 
Below examples showing which residues are the most frequently involved in water 
bridges formation (:func:`.calcBridgingResiduesHistogram`), details of that 
interactions (:func:`.calcWaterBridgesStatistics`), and results saved as PDB 
structure for further visualization (:func:`.savePDBWaterBridgesTrajectory`). 
Other functions can be seen in the analysis of trajectory.


.. ipython:: python

   calcBridgingResiduesHistogram(bridgeFrames_ens)


.. ipython:: python

   analysisAtomic_ens = calcWaterBridgesStatistics(bridgeFrames_ens, coords_ens)

   for item in analysisAtomic_ens.items():
      print(item)


.. ipython:: python

   savePDBWaterBridgesTrajectory(bridgeFrames_ens, coords_ens, ens[:-4]+'_ens.pdb')



Detecting water centers
-------------------------------------------------------------------------------

Previous function generated multiple PDB in which we would found protein and 
water molecules for each frame that form water bridges with protein structure. 
Now we can use another function :func:`.findClusterCenters` which will extract 
water centers (they refer to the oxygens from water molecules that are forming 
clusters). We need to provide a file pattern as show below. Now all the files 
with prefix *'pebp1_50frames_ens_'* will be analyzed.


.. ipython:: python

   findClusterCenters('pebp1_50frames_ens_*.pdb')


Function generated one PDB file with water centers. We used default values, 
such as distC (distance to other molecule) and numC (min number of molecules 
in a cluster), but those values could be changed if the molecules are more 
widely distributed or we would like to have more numerous clusters.
Moreover, this function can be applied on different type of molecules by using 
*selection* paramater. We can provide the whole molecule and by default 
the center of mass will be used as a reference.


Saved PDB files using :func:`.savePDBWaterBridgesTrajectory` in the previous
step can be upload to VMD or other program for visualization:

.. figure:: images/Fig4.tga
   :scale: 50 %


After uploading a new PDB file with water centers we can see the results as
follows:

.. figure:: images/Fig4.tga
   :scale: 50 %
