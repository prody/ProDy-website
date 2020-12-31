Calculations
===============================================================================

Here are the required imports again. You do not need to repeat them if you are
still in the same python session.

.. ipython:: python

    from prody import *
    from pylab import *
    ion()


Loading the structure and initialising adaptive ANM
-------------------------------------------------------------------------------
First, we parse the structures that we want to analyse with Adaptive ANM.
For this tutorial, we will fetch chain A of PDB structures 1GRU and 1GR5. 
We import a subset containing the calpha atoms, which we use in downstream steps.

.. ipython:: python

    r_pp = parsePDB('1gruA', subset='ca')
    t = parsePDB('1gr5A', subset='ca')


Next, we create an :class:`.AdaptiveANM` object, which will be used for the calculations.

.. ipython:: python

    aanm = AdaptiveANM([r_pp, t])

This new implementation of Adaptive ANM increases the number of modes calculated 
each cycle depending on which ones were used in the previous cycles and their 
overlaps. We can give it keyword arguments for the maximum number of modes to use
(3N-6 by default for N atoms), starting number of modes (20 by default), and other 
ANM parameters. These parameters can also be presented to the functions below.

It also has 3 ways of running calculations, each of which will be demonstrated below.

Standard Adaptive ANM: Alternating Steps
-------------------------------------------------------------------------------

As in the original Adaptive ANM ([ZY09]_), we can alternate between which structure 
to use for the ANM calculations and mode selection. This is done using the class 
method :meth:`.runManyStepsAlternating`. We provide a large number of steps to ensure 
convergence is reached.

.. ipython:: python

    aanm.runManyStepsAlternating(100)

One way Adaptive ANM
-------------------------------------------------------------------------------

If one structure is very compact and unlikely to undergo much conformational change 
along ANM modes, you can just use the other structure for driving the transition with 
the class method :meth:`.runManySteps`.

.. ipython:: python

    aanm.runManySteps(100)

Alternatively, we can run one structure as far it can go and then run the rest of the steps 
with the other structure with the class method :meth:`.runManySteps`. In this case, it helps 
to provide the maximum number of modes to control how far it goes before switching structures.

.. ipython:: python

    aanm.runManyStepsFurthestEachWay(100, maxModes=30)
