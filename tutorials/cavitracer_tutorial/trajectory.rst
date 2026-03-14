.. _cavitracer_single:

Detection of channels in molecular dynamics (MD) trajectory
===============================================================================

Analysis of the trajectory will be performed on a short MD trajectory
containing a few frames of simulation performed for 15-lipoxygenase from
P.aeruginosa (pLoxA). This protein contains 665 residues and a catalytic center
with iron. During the analysis, the metal center is ignored.

Before analyzing the trajectory, its need to be parsed (see more details
in `Trajectory tutorial`_).

.. ipython:: python
   :verbatim:

   PDBfile = 'pLoxA2.pdb'
   DCDfile = 'pLoxA2_ev5.dcd'
   atoms = parsePDB(PDBfile)
   dcd = Trajectory(DCDfile)
   dcd.link(atoms)
   dcd.setCoords(atoms)

.. parsed-literal::

   @> 10298 atoms and 1 coordinate set(s) were parsed in 0.11s.

To detect channels in MD trajectories, we need to use
:func:`.scalcChannelsMultipleFrames` function.

.. ipython:: python
   :verbatim:

   channels4, surfaces4=calcChannelsMultipleFrames(atoms, dcd, output_path = 'channels_pLoxA_dcd', separate=True)

.. parsed-literal::

   @> Frame: 0
   @> Detected 23 channels.
   @> Saving multiple results to directory ..
   @> Frame: 1
   @> Detected 21 channels.
   @> Saving multiple results to directory ..
   @> Frame: 2
   @> Detected 20 channels.
   @> Saving multiple results to directory ..
   @> Frame: 3
   @> Detected 17 channels.
   @> Saving multiple results to directory ..
   @> Frame: 4
   @> Detected 17 channels.
   @> Saving multiple results to directory ..  

To have access to a particular frame:

.. ipython:: python
   :verbatim:

   frame2 = dcd.getFrame(2)
   frame2

.. parsed-literal::

   <Frame: 2 from pLoxA2_ev5 (10298 atoms)>

.. ipython:: python
   :verbatim:

   getChannelResidueNames(frame2, channels4[1], residues_file_name='pLoxA_DCD_res')

.. parsed-literal::

   @> 10633 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10438 atoms and 1 coordinate set(s) were parsed in 0.09s.
   @> 10793 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10618 atoms and 1 coordinate set(s) were parsed in 0.09s.
   @> 10663 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10638 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10543 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10513 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10678 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10588 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10753 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10753 atoms and 1 coordinate set(s) were parsed in 0.11s.
   @> 10463 atoms and 1 coordinate set(s) were parsed in 0.11s.
   @> 10423 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10668 atoms and 1 coordinate set(s) were parsed in 0.11s.
   @> 10753 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10633 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10368 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10513 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10433 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> 10358 atoms and 1 coordinate set(s) were parsed in 0.09s.

   ['channel0: PHE83, GLY84, LEU154, GLU155, LYS158, ASN159, GLU373, MET374, LEU378, HSD382, GLY412, THR413, ILE416, ASN417, ALA420, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, ALA551, HSD555, PHE560, TYR616, HSP617, GLY618, TYR623, ARG624, GLN625, THR626, GLY627, PHE628, VAL633, PRO680, ALA681, SER682, THR683, ASN684, ILE685',
    'channel1: THR80, PHE83, GLY84, ALA92, ARG93, GLU373, MET374, LEU378, HSD382, LEU383, GLU386, ILE416, ASN417, ALA420, ALA421, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, ILE685',
    'channel2: GLN36, ILE39, ASP40, LEU49, LEU54, PRO55, GLN56, ARG66, ARG67, VAL69, LEU70, LYS73, PHE83, GLY84, GLU373, MET374, LEU378, HSD382, ALA397, PRO398, HSE403, VAL404, ALA407, GLY412, THR413, ILE416, ASN417, ALA420, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, ASP508, ALA511, ASP512, VAL513, GLU514, ALA551, HSD555, PHE560, TYR616, HSP617, GLY618, TYR623, ARG624, VAL633, PRO680, ALA681, SER682, THR683, ASN684, ILE685',
    'channel3: PHE83, GLY84, GLU373, MET374, LEU378, HSD382, PRO408, GLY412, THR413, ILE416, ASN417, ALA420, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, ALA551, HSD555, PHE560, TYR616, HSP617, GLY618, ASP622, TYR623, ARG624, PRO632, VAL633, PHE634, THR640, LEU646, PRO680, ALA681, SER682, THR683, ASN684, ILE685',
    'channel4: PHE83, GLY84, ILE113, LEU116, PHE120, VAL189, LEU193, THR201, ALA295, PRO296, SER297, GLY298, ALA299, VAL300, LYS302, THR305, GLY306, THR307, GLY308, PHE309, ALA310, HSD340, GLU373, MET374, LEU378, HSD382, ILE416, ALA420, ILE423, LEU424, LEU425, PHE430, ILE431, ASP432, VAL433, MET434, PHE435, ALA436, ALA437, PRO438, THR442, GLN443, ALA446, MET575, PRO579, PRO598, LEU600, VAL601, GLU604, ILE608, LEU611, LEU612, ILE685',
    'channel5: GLN36, ILE39, ASP40, LEU49, PHE83, GLY84, GLU373, MET374, LEU378, HSD382, PRO398, HSE403, VAL404, GLY412, THR413, ILE416, ASN417, ALA420, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, ALA551, HSD555, PHE560, TYR616, HSP617, GLY618, TYR623, ARG624, VAL633, PRO680, ALA681, SER682, THR683, ASN684, ILE685',
    'channel6: PHE83, GLY84, ILE113, LEU116, ASN119, PHE120, SER123, VAL189, LEU193, GLN205, GLU373, MET374, LEU378, HSD382, ILE416, ALA420, ILE423, LEU424, LEU425, PHE430, ILE431, THR442, GLN443, ALA446, LEU600, LEU603, GLU604, ASN607, ILE608, LEU611, LEU612, ILE685',
    'channel7: PHE83, GLY84, GLU373, MET374, LEU378, HSD382, GLY412, THR413, ILE416, ASN417, ALA420, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, ALA551, HSD555, ASN559, PHE560, SER614, VAL615, TYR616, ARG678, PRO680, SER682, THR683, ASN684, ILE685',
    'channel8: ALA41, SER42, LEU45, LEU46, PHE83, GLY84, ALA164, GLY165, LEU168, LEU169, GLU373, MET374, LEU378, HSD382, GLY412, THR413, ILE416, ASN417, ALA420, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, ALA551, HSD555, PHE560, TYR616, HSP617, GLY618, TYR623, ARG624, PHE628, PRO629, VAL633, PRO680, ALA681, SER682, THR683, ASN684, ILE685',
    'channel9: PHE83, GLY84, ILE113, LEU116, PHE120, VAL189, LEU193, LYS194, SER197, THR201, LEU303, LEU304, THR305, GLU373, MET374, LEU378, HSD382, ILE416, ALA420, ILE423, LEU424, LEU425, PHE430, ILE431, ASP432, VAL433, MET434, THR442, GLN443, ALA446, PRO598, LEU600, VAL601, GLU604, ILE608, LEU611, LEU612, ILE685',
    'channel10: PHE83, GLY84, ILE113, LEU116, PHE120, VAL189, LEU193, THR201, ALA295, PRO296, SER297, GLY298, ALA299, VAL300, LYS302, THR305, GLY306, THR307, GLY308, PHE309, ALA310, CYS333, HSD340, PRO341, MET342, PHE343, VAL344, LEU353, GLU373, MET374, LEU378, HSD382, ILE416, ALA420, ILE423, LEU424, LEU425, PHE430, ILE431, ASP432, VAL433, MET434, PHE435, ALA436, ALA437, PRO438, THR442, GLN443, ALA446, MET575, PRO579, ALA580, PRO581, ASP582, PRO584, PRO598, LEU600, VAL601, GLU604, ILE608, LEU611, LEU612, ILE685',
    'channel11: GLN36, ILE39, ASP40, LEU49, ARG66, ARG67, VAL69, LEU70, LYS73, LYS74, PHE83, GLY84, GLU373, MET374, LEU378, HSD382, THR395, ALA397, PRO398, HSE403, VAL404, GLY412, THR413, ILE416, ASN417, ALA420, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, VAL513, GLU514, ALA517, GLU521, ALA551, HSD555, PHE560, TYR616, HSP617, GLY618, TYR623, ARG624, VAL633, PRO680, ALA681, SER682, THR683, ASN684, ILE685',
    'channel12: PHE83, GLY84, LEU106, ILE179, ALA180, LEU182, THR183, GLU373, MET374, LEU378, HSD382, PHE415, ILE416, GLY419, ALA420, ARG422, ILE423, LEU424, LEU425, ILE431, THR442, GLN443, ALA446, LEU612, ILE685',
    'channel13: ARG213, VAL224, ALA225, PHE228, TYR311, GLN358, THR362, VAL363, VAL366, TYR568, ALA569, PRO570, ALA571, ILE572, CYS573, SER576, TRP592, MET595, MET596, PRO597',
    'channel14: GLU223, VAL224, SER227, PHE228, ASP230, ASP231, GLU232, ALA233, PHE234, ALA235, TYR236, VAL239, ARG282, TYR311, LEU319, GLY320, LYS321, ASP322, ALA324, ARG325, LEU326, LEU327, GLN358, LYS361, THR362, VAL363, GLN365, VAL366, LEU473, ALA474, LEU475, PRO476, ASP477, ALA569, PRO570, ALA571, ILE572, CYS573, SER576, TRP592, MET595, MET596, PRO597, ARG669',
    'channel15: GLU223, VAL224, SER227, PHE228, ASP230, ASP231, GLU232, ALA233, PHE234, ALA235, TYR236, VAL239, GLU261, PHE264, ARG265, ARG266, MET268, GLY269, ALA270, ASP272, SER273, LEU274, ARG282, TYR311, LEU319, GLY320, LYS321, ASP322, ALA324, ARG325, LEU326, LEU327, GLN358, LYS361, THR362, VAL363, GLN365, VAL366, LEU473, ALA474, LEU475, PRO476, ASP477, ALA569, PRO570, ALA571, ILE572, CYS573, SER576, TRP592, MET595, MET596, PRO597, ARG669',
    'channel16: GLU223, VAL224, SER227, PHE228, ASP230, ASP231, GLU232, ALA233, PHE234, ALA235, TYR236, VAL239, TYR311, ARG325, LEU326, LEU327, GLN358, LYS361, THR362, VAL363, GLN365, VAL366, ALA569, PRO570, ALA571, ILE572, CYS573, SER576, TRP592, MET595, MET596, PRO597, ALA666, ARG667, ARG668, ARG669',
    'channel17: PRO479, ASP483, ILE660, ARG661, ASN664, ARG667, TYR671, GLU672, LEU674, LEU675',
    'channel18: ALA379, GLN380, THR381, VAL384, SER385, PHE388, LEU405, LEU406, PRO408, HSD409, PHE410, ASN449, ARG450, GLY452, PHE453, PHE455, ILE491, TRP494, LEU544, MET546, VAL547, ILE548, THR550, ALA551, LEU620, TYR623, PHE649',
    'channel19: PHE228, ARG229, ASP230, ASP231, PHE234, ARG238, LEU286, TYR288, PRO313, ILE314, ALA315, PHE317, PRO328, ILE331, ARG345, PRO346, TYR354, TRP357, LYS361, VAL364',
    'channel20: VAL114, ILE117, VAL118, VAL121, GLU133, ILE135, ALA136, THR137, LEU139, SER140, VAL185, LEU188']

.. _Trajectory tutorial: http://www.bahargroup.org/prody/tutorials/trajectory_analysis/