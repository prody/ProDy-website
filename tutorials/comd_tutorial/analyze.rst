.. _analyze:

Collective Molecular Dynamics Analysis
======================================

After running coMD simulation, the results will be prepared in the same 
folder named in setup section. You need to download the following files
from that folder:

1. :file:`initial_ionized.pdb`
2. :file:`final_ionized.pdb`
3. :file:`initial_trajectory.pdb`
4. :file:`final_trajectory.pdb`

Preparing Trajectory Files
--------------------------

We recommend that you will follow this tutorial by typing commands in an
IPython session, e.g.::

  $ ipython

or with pylab environment::

  $ ipython --pylab


First, we will make necessary imports from ProDy and Matplotlib
packages.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

coMD simulations will create two different trajectory and we need to 
use those two trajectories to analyze our simulations. However, the
simulation boxes have different number of atoms due to solvation and
ionization procedure in the beginning of simulations. We need to load
trajectories and related structure files. 

.. ipython:: python

	dcd1 = Trajectory('initial_trajectory.dcd')
	dcd2 = Trajectory('final_trajectory.dcd')
	dcd1
	dcd2
	structure1 = parsePDB('initial_ionized.pdb')
	structure2 = parsePDB('final_ionized.pdb')
	structure1
	structure2

It is required to link the structure and atoms however, the size of
initial and final structure may be different in terms of total atoms. 
However, the total number of C-alpha atoms must be equal. The structure 
and trajectory is linked with following commands:

.. ipython:: python

	dcd1.setCoords(structure1)
	dcd1.setAtoms(structure1.calpha)
	dcd1
	dcd2.setCoords(structure2)
	dcd2.setAtoms(structure2.calpha)
	dcd2

It is not available to concatanate trajectories inside Python (It is available 
for ProDy command line application). Therefore, we prepare two trajectory files for only C-alpha atoms with following commands:

.. ipython:: python

	writeDCD('initial_filtered.dcd', dcd1)
	writeDCD('final_filtered.dcd', dcd2)

Concatenating Trajectory Files
------------------------------

After having two trajectory files, we can start analyzing multiple :files:`.dcd` files. First, the trajectory files will be concatanated. 

.. ipython:: python

	traj = Trajectory('initial_filtered.dcd')
	traj.addFile('final_filtered.dcd')

Principal Component Analysis
----------------------------

For concatanated trajectory we created a PCA object, created covariance matrix and calculated eigenvalues and eigenvectors. 

.. ipython:: python

	pca = PCA('Adelynate Kinase coMD')
	pca.buildCovariance(traj)
	pca.calcModes()

The first half of the trajectory is for initial structure and the second half of the trajectory is for final structure. Those two trajectories are seperated. 

.. ipython:: python

	forward = traj[0:40]
	backward = traj[40:]

Before visualizing the trajectories on principal components, it is required to translate those structures on top by using :func:`.superpose` function. 

.. ipython:: python

	forward.superpose()
	backward.superpose()

Visualization of Trajectories
-----------------------------

Finally, the trajectories can be plotted by using :func:`showProjection` function:

.. ipython:: python

	showProjection(forward, pca[:3], color='red', marker='.');
	showProjection(backward, pca[:3], color='blue', marker='.');
	showProjection(forward[0], pca[:3], color='red', marker='o');
	showProjection(backward[0], pca[:3], color='blue', marker='o');

The plots will be in the following form: 

.. figure:: _static/figures/comd_3d_out.png
	:scale: 80%

Now we calculated the modes and we can write them to a :file:`.nmd` file for viewing in normal mode wizard. 

.. ipython:: python
	
	writeNMD('ake_pca.nmd', pca[:3], structure1.select('calpha'))
