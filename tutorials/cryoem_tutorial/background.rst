.. _background:

Cryo-EM Microscopy and Analysis
==============================================================
Cryo-electron microscopy (cryo-EM) is a type of transmission 
electron microscopy where the studied sample, e.g. protein 
solution, is rapidly frozen. In structural biology, single 
particle cryo-EM is gaining more and more popularity due to 
its advantage of allowing the observation of protein under 
the native environment and improved resolution, and as 
a result, the structures of numerous mega-Dalton biomolecules 
are captured and stored in the form of cryo-EM density maps, 
with resolutions ranging from 2 to 100 Å (maybe link 
Emdatabank here: http://www.emdatabank.org/). Fitting all-atom 
structures to the density map is time-consuming, and sometimes 
not feasible due to low resolutions. In addition, for the 
purpose of studying the dynamics, MD simulations of 
supercomplexes with atomic details are extremely 
computationally expensive. It is thus desirable to apply 
coarse-grained methods, for example Anisotropic Network Model 
(ANM), to such data to study the “big motions” of molecular 
machines beyond the atomic level. ProDy can be used for 
constructing the bead-and-spring model from cryo-EM data 
(using a published algorithm*), which can be further used 
for either coarse-grained MD simulation or ANM.

