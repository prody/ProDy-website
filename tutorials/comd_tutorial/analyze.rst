.. _analyze:

Collective Molecular Dynamics Analysis
======================================

After running a coMD simulation, the results will be prepared in the same 
folder named in the setup section. You need to download the following files
from that folder:

1. :file:`initial_ionized.pdb`
2. :file:`final_ionized.pdb`
3. :file:`initial_trajectory.dcd`
4. :file:`final_trajectory.dcd`

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

coMD simulations will create two different trajectories and we need to 
use those two trajectories to analyze our simulations. However, the
simulation boxes have different number of atoms due to the solvation and
ionization procedure in the beginning of simulations. We therefore need to 
load trajectories with their related structure files. 

.. ipython:: python

	dcd1 = Trajectory('initial_trajectory.dcd')
	dcd2 = Trajectory('final_trajectory.dcd')
	dcd1
	dcd2
	structure1 = parsePDB('initial_ionized.pdb')
	structure2 = parsePDB('final_ionized.pdb')
	structure1
	structure2

In order to analyze the two trajectories together, we must have the 
same number of atoms in the sets that we analyze. To do this we link 
the trajectories to there corresponding initial structures and set 
the C-alpha atoms as the active ones with following commands:

.. ipython:: python

	dcd1.setCoords(structure1)
	dcd1.setAtoms(structure1.calpha)
	dcd1
	dcd2.setCoords(structure2)
	dcd2.setAtoms(structure2.calpha)
	dcd2

Concatenating Trajectory Files
------------------------------

It is not possible to concatenate Trajectory objects inside Python but 
must instead be done using files. We therefore write out new files 
filtered to contain only C-alpha atoms as set in the previous step:

.. ipython:: python

	writeDCD('initial_filtered.dcd', dcd1)
	writeDCD('final_filtered.dcd', dcd2)

One way to combined trajectories into the same object is the :meth:`.Trajectory.addFile` method 
of :class:`.Trajectory` objects.

.. ipython:: python

	traj = Trajectory('initial_filtered.dcd')
	traj.addFile('final_filtered.dcd')

Alternatively we can create an :class:`.Ensemble` using :func:`.parseDCD`,
which gives us the flexibility to do things like reversing the final 
trajectory to create something we can view in VMD_ rather than having
the trajectories both run towards the shared intermediate. 
We do this as follows:

.. ipython:: python

	combined_traj = parseDCD('initial_filtered.dcd')
	w2_traj = parseDCD('final_filtered.dcd')

    for i in reversed(range(len(w2_traj))):
        combined_traj.addCoordset(w2_traj.getConformation(i))

	combined_traj.superpose()

    writeDCD('combined_trajectory.dcd', combined_traj)

We also write out a :file:`pdb` file containing just the C-alpha atoms, which can 
be loaded into VMD_ together with this combined trajectory for visualization.

.. ipython:: python

	writePDB('initial_filtered.pdb', structure1.ca)


Principal Component Analysis
----------------------------

We next perform PCA on the concatenated trajectory as follows. 

.. ipython:: python

	pca = PCA('Adelynate Kinase coMD')
	pca.buildCovariance(combined_traj)
	pca.calcModes()

The first half of the trajectory is from the initial structure and the second half of the trajectory is from the final structure. 
We can identify these two trajectories as follows. 

.. ipython:: python

	forward = combined_traj[0:40]
	backward = combined_traj[40:]


Visualization of Trajectories
-----------------------------

Finally, the trajectories can be plotted by using the :func:`showProjection` function:

.. ipython:: python

	showProjection(forward, pca[:3], color='red', marker='.');
	showProjection(backward, pca[:3], color='blue', marker='.');
	showProjection(forward[0], pca[:3], color='red', marker='o');
	showProjection(backward[0], pca[:3], color='blue', marker='o');

The plots will be in the following form: 

.. figure:: images/comd_3d_out.png
	:scale: 80%

Having calculated the modes, we can write them to a :file:`.nmd` file for viewing in NMWiz_. 

.. ipython:: python
	
	writeNMD('ake_pca.nmd', pca, structure1.ca)

.. _NMWiz: http://prody.csb.pitt.edu/nmwiz/