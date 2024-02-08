Ensemble analysis: dynamics from PCA of multiple structures.
===============================================================================

This example shows how to run an ensemble analysis workflow in Scipion by 
importing several structures and building an ensemble, and analysing it via 
PCA. We will use the SARS-CoV-2 spike as in [KJ23]_.

See [KJ23]_ for more information about the Scipion-EM-ProDy plugin and more examples.

Getting started
-------------------------------------------------------------------------------
You can open the Scipion Projects window from the command line as follows::

   scipion3

From here, create a new project and call it prody_tutorial_ensemble. Then,
select the ProDy protocols tree from the dropdown on the left.

Importing structures and building an ensemble
-------------------------------------------------------------------------------

As with any ProDy pipeline, we start by parsing atomic structures. In this case, we use the second protocol from the Scipion core 
to importing multiple atomic structures at once. 

.. figure:: images/core/1_import_set_of_atomstructs.png
   :scale: 80
