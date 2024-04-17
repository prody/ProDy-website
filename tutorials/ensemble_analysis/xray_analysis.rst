.. _pca-xray-analysis:


Analysis
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example is continued from :ref:`pca-xray-calculations`.  The aim of this
part is to perform a quantitative comparison of experimental and theoretical
data and to print/save the numerical data that were presented in [AB09]_.


We start by importing everything from the ProDy package and loading data 
from the previous part. These steps can be skipped if you are continuing 
in the same iPython session.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

   pca = loadModel('p38_xray.pca.npz')
   anm = loadModel('5uoj.anm.npz')

Variance along PCs
-------------------------------------------------------------------------------

Of interest is the fraction of variance that is explained by principal
components, which are the dominant modes of variability in the dataset.
We can print this information to screen for top 3 PCs as follows:

.. ipython:: python

   for mode in pca[:3]:
       var = calcFractVariance(mode)*100
       print('{0:s}  % variance = {1:.2f}'.format(mode, var))

These data were included in Table 1 in [AB09]_.

Collectivity of modes
-------------------------------------------------------------------------------

Collectivity of a normal mode ([BR95]_) can be obtained using
:func:`.calcCollectivity`:

.. ipython:: python

   for mode in pca[:3]:    # Print PCA mode collectivity
       coll = calcCollectivity(mode)
       print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))

We can also calculate the collectivity of ANM modes:

.. ipython:: python

   for mode in pca[:3]:    # Print PCA mode collectivity
       coll = calcCollectivity(mode)
       print('{0:s}  collectivity = {1:.2f}'.format(mode, coll))

This shows that top PCA and ANM modes are highly collective.


Save numeric data
-------------------------------------------------------------------------------

:class:`.ANM` and :class:`.PCA` instances store calculated numeric data.
Their class documentation lists methods that return eigenvalue, eigenvector,
covariance matrix etc. data to the user. Such data can easily be written into
text files for analysis using external software. The function is to use is
:func:`.writeArray`:

.. ipython:: python

   writeArray('p38_PCA_eigvecs.txt', pca.getEigvecs() ) # PCA eigenvectors
   writeModes('p38_ANM_modes.txt', anm) # This function is based on writeArray

