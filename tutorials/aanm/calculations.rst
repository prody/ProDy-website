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
For this tutorial, we will fetch chain A of PDB structures 1GRU and 1GR5, 
which correspond to the R'' and T states, which we will call r and t.

We import a subset containing the calpha atoms, which we use in downstream steps.

.. ipython:: python

    r = parsePDB('1gruA', subset='ca')
    t = parsePDB('1gr5A', subset='ca')

In order to have the same number of atoms in both structures, we use the function 
:func:`.matchChains` to create :class:`.AtomMap` objects with the same number. 

As r has more atoms (524) than t (517), we provide t first to find the common atoms 
that are present in t:

.. ipython:: python

    t_amap, r_amap, ident, cover = matchChains(t, r)[0]

This new implementation of Adaptive ANM increases the number of modes calculated 
each cycle depending on which ones were used in the previous cycles and their 
overlaps. We can give it keyword arguments for the maximum number of modes to use
(3N-6 by default for N atoms), starting number of modes (20 by default), and other 
ANM parameters. These parameters can also be presented to the functions below.

It also has 3 ways of running calculations, each of which will be demonstrated below.

1. Standard Adaptive ANM: Alternating Steps
-------------------------------------------------------------------------------

As in the original Adaptive ANM ([ZY09]_), we can alternate between which structure 
to use for the ANM calculations and mode selection. This is done using the class 
method :meth:`.runManyStepsAlternating`. We provide a large number of steps to ensure 
convergence is reached.

Next, we create an :class:`.AdaptiveANM` object, which will be used for the calculations.

.. ipython:: python

    aanm_1 = runAdaptiveANM(r_amap, t_amap, 20)

2. One way Adaptive ANM
-------------------------------------------------------------------------------

If one structure is very compact and unlikely to undergo much conformational change 
along ANM modes, you can just use the other structure for driving the transition with 
the class method :meth:`.runManySteps`.

.. ipython:: python

    aanm_2 = runOneWayAdaptiveANM(r_amap, t_amap, 20)

3. One way at a time Adaptive ANM
-------------------------------------------------------------------------------

Alternatively, we can run one structure as far it can go and then run the rest of the steps 
with the other structure with the class method :meth:`.runManySteps`. In this case, it helps 
to provide the maximum number of modes to control how far it goes with each structure.

.. ipython:: python

    aanm_3 = runBothWaysAdaptiveANM(r_amap, t_amap, 20, maxModes=3)
