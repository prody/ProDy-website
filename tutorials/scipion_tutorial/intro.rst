Introduction
===============================================================================

This tutorial shows how to use ProDy within Scipion including experimental ensemble analysis,
NMA and ClustENMD simulations (see Figure 3 of [KJ23]_).


Required Programs
-------------------------------------------------------------------------------

Scipion_ 3 is required. How to install it can be found on its website.

.. _Scipion: https://scipion-em.github.io/docs/release-3.0.0/

Installation
-------------------------------------------------------------------------------

To install the latest release version of the scipion-em-prody plugin (3.3.0), you can use the following command::
   
   scipion3 installp -p scipion-em-prody

This will also install the latest release version (2.4.1) of ProDy into environments that Scipion can use
along with TEMPy and OpenMM.

The latest development version of the plugin can also be installed as follows::

   git clone https://github.com/scipion-em/scipion-em-prody.git
   scipion3 installp -p ./scipion-em-prody --devel

There is also the option to use the release or development version of ProDy with the plugin.
This can be achieved by using commands like the following::

   scipion3 installb prody-2.4.1
   scipion3 installb prody-github

Alternatively, you can activate the conda environments with one of these names and also the main scipion3 environment
and pip install your own ProDy there.

You can generate some example projects and test that everything is working by running the following::

   scipion3 test prody2.tests

Getting Started in Scipion
-------------------------------------------------------------------------------

You can start Scipion from the command line as follows::

   scipion3

This will open the Projects window. From here, you can create a new project or open and existing one.

.. figure:: images/scipion_projects.png
   :scale: 80%

Once you create or select a project, you will get a new window for your project as you can see for TestProDyClustENMsingle::

.. figure:: images/scipion_TestProDyClustENMsingle.png
   :scale: 80%

You can also open up the last project that was open or changed as follows::

   scipion3 last

We will look at the interface more on the next pages but you can also find more information 
at https://scipion-em.github.io/docs/release-3.0.0/docs/user/scipion-gui.html.

.. [KJ23] Krieger JM, Sorzano COS, Carazo JM.
   Scipion-EM-ProDy: A Graphical Interface for the ProDy Python Package within the Scipion Workflow Engine 
   Enabling Integration of Databases, Simulations and Cryo-Electron Microscopy Image Processing. *Int. J. Mol. Sci.* **2023** 24(18):14245.