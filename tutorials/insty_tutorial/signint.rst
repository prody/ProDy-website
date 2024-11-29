.. _insty_tutorial:
=======

Signature Interactions of ensembles
===============================================================================

We have developed a pipeline that leverages BLAST, Dali, or Foldseek search to
identify and download homologs for a particular PDB structure. The
identified structures are prepared for InSty analysis by applying the following
steps: (i) download PDB (ii), extract particular chains, (iii) add
hydrogens/side chains, (iv) perform structural alignment (v) create a
folder and put their prepared structures.

As an example, we will use Aurora kinase A structure(PDB: **1OL5**).


BLAST approach
-------------------------------------------------------------------------------

We will use the PDB code and chain ID to access information from the BLAST server.
Additionally, we can also define many different types of parameters. Here, we
will define sequence identity (seqid) and the method to add missing hydrogen bonds
and side chains (fixer; it can be PDBFixer or Openbabel). We can also define a name for
a folder in which all the downloaded files will be uploaded. 

To download homologs, add missing hydrogens and align structures, we will
use :func:`.runBLAST` function:


.. ipython:: python
   :verbatim:

   PDBcode = '1ol5'
   runBLAST(PDBcode, 'A', seqid=70, fixer='openbabel', folder_name='struc_homologs_BLAST')

.. parsed-literal::

   @> PDB file is found in working directory (1ol5.pdb).
   @> 2607 atoms and 1 coordinate set(s) were parsed in 0.03s.
   @> Blast searching NCBI PDB database for "SKKRQ..."
   @> Blast search completed in 46.3s.                     
   @> Separating chains and saving into PDB file
   @> PDB code 1muo and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 1muo downloaded (1muo.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2029 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> PDB code 2j4z and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 2j4z downloaded (2j4z.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 4626 atoms and 1 coordinate set(s) were parsed in 0.14s.
   @> PDB code 7o2v and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 7o2v downloaded (7o2v.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2265 atoms and 1 coordinate set(s) were parsed in 0.10s.
   @> PDB code 6c83 and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 6c83 downloaded (6c83.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 4937 atoms and 1 coordinate set(s) were parsed in 0.16s.
   @> PDB code 6cpe and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 6cpe downloaded (6cpe.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 4024 atoms and 1 coordinate set(s) were parsed in 0.13s.
   @> PDB code 6cpf and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 6cpf downloaded (6cpf.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 4100 atoms and 1 coordinate set(s) were parsed in 0.13s.
   @> PDB code 6cpg and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 6cpg downloaded (6cpg.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 9962 atoms and 1 coordinate set(s) were parsed in 0.20s.
   @> PDB code 8sso and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 8sso downloaded (8sso.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 6151 atoms and 1 coordinate set(s) were parsed in 0.18s.
   @> PDB code 2xng and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 2xng downloaded (2xng.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2001 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> PDB code 4b0g and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 4b0g downloaded (4b0g.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2065 atoms and 1 coordinate set(s) were parsed in 0.07s.
   @> PDB code 8ssp and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 8ssp downloaded (8ssp.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 3066 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> PDB code 2x6d and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 2x6d downloaded (2x6d.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2144 atoms and 1 coordinate set(s) were parsed in 0.08s.
   @> PDB code 2x6e and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 2x6e downloaded (2x6e.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2007 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> PDB code 4byi and chain A
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 4byi downloaded (4byi.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 2049 atoms and 1 coordinate set(s) were parsed in 0.09s.
   @> PDB code 4byj and chain A
   ..
   ..
   @> Connecting wwPDB FTP server RCSB PDB (USA).
   @> Downloading PDB files via FTP failed, trying HTTP.
   @> 5k3y downloaded (5k3y.pdb.gz)
   @> PDB download via HTTP completed (1 downloaded, 0 failed).
   @> 5831 atoms and 1 coordinate set(s) were parsed in 0.19s.
   @> Adding hydrogens to the structures..
   @> Hydrogens were added to the structure. New structure is saved as addH_1muoA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_2j4zA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_7o2vA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_6c83A.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_6cpeA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_6cpfA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_6cpgA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_8ssoA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_2xngA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4b0gA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_8sspA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_2x6dA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_2x6eA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4byiA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4byjA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4zs0A.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4ztqA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4ztrA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4ztsA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_5oneA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_6graA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_6z4yA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_7ayhA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_7ayiA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_5zanA.pdb.
   ..
   ..
   @> Hydrogens were added to the structure. New structure is saved as addH_2vgoA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_2vgpA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_3ztxA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_2vrxA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4c2vA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_4c2wA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_6gr8A.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_6gr9A.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_2bfyA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_5eykA.pdb.
   @> Hydrogens were added to the structure. New structure is saved as addH_5k3yA.pdb.
   @> 209 PDBs were parsed in 8.49s.        
   @> Aligning the structures..
   @> addH_2j4zA
   @> Checking AtomGroup addH_2j4zA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_2j4zA (len=263) and Chain A from addH_1muoA (len=251):
   @> 	Match: 251 residues match with 100% sequence identity and 95% overlap.
   @> Aligning the structures..
   @> addH_7o2vA
   @> Checking AtomGroup addH_7o2vA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_7o2vA (len=264) and Chain A from addH_1muoA (len=251):
   @> 	Match: 250 residues match with 100% sequence identity and 95% overlap.
   @> Aligning the structures..
   @> addH_6c83A
   @> Checking AtomGroup addH_6c83A: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_6c83A (len=248) and Chain A from addH_1muoA (len=251):
   @> 	Match: 246 residues match with 99% sequence identity and 98% overlap.
   @> Aligning the structures..
   @> addH_6cpeA
   @> Checking AtomGroup addH_6cpeA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_6cpeA (len=256) and Chain A from addH_1muoA (len=251):
   @> 	Match: 248 residues match with 99% sequence identity and 97% overlap.
   @> Aligning the structures..
   @> addH_6cpfA
   @> Checking AtomGroup addH_6cpfA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_6cpfA (len=259) and Chain A from addH_1muoA (len=251):
   @> 	Match: 249 residues match with 99% sequence identity and 96% overlap.
   @> Aligning the structures..
   @> addH_6cpgA
   @> Checking AtomGroup addH_6cpgA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_6cpgA (len=249) and Chain A from addH_1muoA (len=251):
   @> 	Match: 248 residues match with 100% sequence identity and 99% overlap.
   @> Aligning the structures..
   @> addH_8ssoA
   @> Checking AtomGroup addH_8ssoA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_8ssoA (len=256) and Chain A from addH_1muoA (len=251):
   @> 	Match: 251 residues match with 100% sequence identity and 98% overlap.
   @> Aligning the structures..
   @> addH_2xngA
   @> Checking AtomGroup addH_2xngA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_2xngA (len=253) and Chain A from addH_1muoA (len=251):
   @> 	Match: 248 residues match with 99% sequence identity and 98% overlap.
   @> Aligning the structures..
   @> addH_4b0gA
   @> Checking AtomGroup addH_4b0gA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_4b0gA (len=247) and Chain A from addH_1muoA (len=251):
   @> 	Match: 241 residues match with 98% sequence identity and 96% overlap.
   @> Aligning the structures..
   @> addH_8sspA
   @> Checking AtomGroup addH_8sspA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_8sspA (len=264) and Chain A from addH_1muoA (len=251):
   @> 	Match: 251 residues match with 100% sequence identity and 95% overlap.
   @> Aligning the structures..
   @> addH_2x6dA
   @> Checking AtomGroup addH_2x6dA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_2x6dA (len=255) and Chain A from addH_1muoA (len=251):
   @> 	Match: 246 residues match with 98% sequence identity and 96% overlap.
   @> Aligning the structures..
   @> addH_2x6eA
   @> Checking AtomGroup addH_2x6eA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_2x6eA (len=249) and Chain A from addH_1muoA (len=251):
   @> 	Match: 245 residues match with 98% sequence identity and 98% overlap.
   @> Aligning the structures..
   @> addH_4byiA
   @> Checking AtomGroup addH_4byiA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_4byiA (len=254) and Chain A from addH_1muoA (len=251):
   @> 	Match: 249 residues match with 99% sequence identity and 98% overlap.
   @> Aligning the structures..
   @> addH_4byjA
   @> Checking AtomGroup addH_4byjA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_4byjA (len=246) and Chain A from addH_1muoA (len=251):
   @> 	Match: 244 residues match with 99% sequence identity and 97% overlap.
   @> Aligning the structures..
   @> addH_4zs0A
   @> Checking AtomGroup addH_4zs0A: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_4zs0A (len=257) and Chain A from addH_1muoA (len=251):
   @> 	Match: 248 residues match with 99% sequence identity and 96% overlap.
   ..
   ..
   @> Aligning the structures..
   @> addH_5eykA
   @> Checking AtomGroup addH_5eykA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_5eykA (len=264) and Chain A from addH_1muoA (len=251):
   @> 	Failed to match chains (seqid=5%, overlap=81%).
   @> Trying to match chains based on local sequence alignment:
   @>  Comparing Chain A from addH_5eykA (len=264) and Chain A from addH_1muoA (len=251):
   @> 	Match: 246 residues match with 70% sequence identity and 93% overlap.
   @> Aligning the structures..
   @> addH_5k3yA
   @> Checking AtomGroup addH_5k3yA: 1 chains are identified
   @> Checking AtomGroup addH_1muoA: 1 chains are identified
   @> Trying to match chains based on residue numbers and names:
   @>   Comparing Chain A from addH_5k3yA (len=268) and Chain A from addH_1muoA (len=251):
   @> 	Failed to match chains (seqid=5%, overlap=81%).
   @> Trying to match chains based on local sequence alignment:
   @>  Comparing Chain A from addH_5k3yA (len=268) and Chain A from addH_1muoA (len=251):
   @> 	Match: 246 residues match with 71% sequence identity and 92% overlap.


To compute all types of interactions for each homolog, we use
:func:`.calcSignatureInteractions` function, and we are providing the name of the
folder with the structures.
 
This function will create additional files with the prefix **'INT_'+type of interactions**
for each file. In such a file except for protein structure, we will have dummy atoms
that will correspond to the interactions. The dummy atoms will be inserted
exactly between the residue-residue pair which is interacting. 
We are computing seven types of non-covalent interactions (hydrogen bonds - HBs,
salt bridges - SBs, repulsive ionic bonding - RIB, pi-cation - PiCat,
pi-stacking - PiStack, hydrophobic interactions - HPh, and disulfide bonds - DiBs).


.. ipython:: python
   :verbatim:

   calcSignatureInteractions('struc_homologs_BLAST')

.. parsed-literal::


   @> struc_homologs_BLAST/align__addH_5l8kA.pdb
   @> 4385 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> Calculating hydrogen bonds.
   @>      DONOR (res chid atom)   <--->       ACCEPTOR (res chid atom)    Distance  Angle
   @>     LEU178    A         N_868  <--->     VAL174    A         O_808     2.8    25.5
   @>     GLU321    A        N_3203  <--->     VAL317    A        O_3141     2.8    27.8
   @>     VAL324    A        N_3257  <--->     TYR320    A        O_3187     2.9    37.4
   @>     LEU208    A        N_1397  <--->     LYS162    A         O_609     2.9    31.8
   @>     ILE360    A        N_3839  <--->     ALA356    A        O_3779     2.9    19.9
   @>     LYS224    A        N_1666  <--->     ARG220    A        O_1596     2.9    39.0
   @>     LEU263    A        N_2315  <--->     LYS271    A        O_2427     2.9    34.6
   @>     LEU318    A        N_3152  <--->     SER314    A        O_3104     2.9    35.1
   @>     GLY173    A         N_796  <--->     GLU170    A         O_754     3.0    38.1
   @>     LEU149    A         N_377  <--->     ARG137    A         O_189     3.0    28.6
   @>     ALA243    A        N_1978  <--->     GLU239    A        O_1925     3.0    23.2
   @>     ILE253    A        N_2140  <--->     SER278    A        O_2541     3.0    23.1
   @>     TYR236    A        N_1866  <--->     ARG232    A        O_1809     3.0     7.4
   @>     SER314    A        N_3099  <--->     VAL310    A        O_3033     3.0    30.1
   @>     ARG220    A        N_1591  <--->     THR217    A      OG1_1548     3.0    17.6
   @>     ILE193    A        N_1133  <--->     HIS190    A        O_1093     3.0    13.9
   @>     VAL317    A        N_3136  <--->     TRP313    A        O_3080     3.0    34.9
   @>     SER342    A        N_3556  <--->     TYR338    A        O_3475     3.0    25.6
   @>     THR337    A        N_3456  <--->     THR333    A        O_3394     3.0    24.4
   @>     TYR320    A        N_3182  <--->     GLY316    A        O_3135     3.1    39.9
   @>     LEU169    A         N_730  <--->     PHE165    A         O_666     3.1    35.0
   @>     LEU315    A        N_3110  <--->     ASP311    A        O_3049     3.1    22.9
   @>     HIS187    A        N_1028  <--->     ILE184    A         O_986     3.1    38.9
   @>     ARG362    A        N_3869  <--->     ASP358    A        O_3813     3.1    23.5
   @>     PHE133    A         N_123  <--->     LEU130    A          O_82     3.1     3.3
   @>     HIS380    A        N_4189  <--->     VAL377    A        O_4144     3.1    34.9
   @>     LEU363    A        N_3893  <--->     LEU359    A        O_3825     3.1    39.5
   @>     GLY316    A        N_3129  <--->     LEU312    A        O_3061     3.1    36.0
   @>     LYS339    A        N_3491  <--->     GLN335    A        O_3429     3.1    15.1
   @>     ILE184    A         N_981  <--->     ARG180    A         O_916     3.1    16.3
   @>     LEU359    A        N_3820  <--->     GLY355    A        O_3773     3.1    37.9
   @>     ARG179    A         N_887  <--->     GLU175    A         O_824     3.2     4.7
   @>     LEU262    A        N_2296  <--->     PRO259    A        O_2257     3.2    25.9
   @>     VAL352    A        N_3722  <--->     PRO349    A        O_3680     3.2    26.7
   @>     VAL377    A        N_4139  <--->     MET373    A        O_4069     3.2    15.7
   @>     GLU183    A         N_966  <--->     ARG179    A         O_892     3.2     6.3
   @>     TYR219    A        N_1570  <--->     PRO259    A        O_2257     3.2    39.1
   @>     LYS309    A       NZ_3024  <--->     PRO372    A        O_4054     3.2    38.5
   @>     ALA172    A         N_786  <--->     GLN168    A         O_718     3.2     6.7
   @>     LEU323    A        N_3238  <--->     CYS319    A        O_3176     3.2    13.8
   @>     GLU221    A        N_1615  <--->     THR217    A        O_1545     3.2    24.1
   @>     ALA241    A        N_1954  <--->     THR238    A        O_1911     3.3    37.2
   @>     GLU336    A        N_3441  <--->     THR333    A        O_3394     3.3    13.4
   @>     ILE301    A        N_2880  <--->     PRO297    A        O_2824     3.3    24.6
   @>     HIS248    A        N_2050  <--->     SER245    A        O_2012     3.4    38.4
   @>     ARG205    A       NE_1351  <--->     ASP202    A      OD2_1311     3.4    29.3
   @>     SER245    A        N_2007  <--->     ASN242    A        O_1969     3.4    31.2
   @>     THR333    A        N_3389  <--->     GLU336    A      OE1_3454     3.4    29.4
   @>     GLY140    A         N_241  <--->     VAL147    A         O_345     3.4    25.0
   @>     LEU225    A        N_1688  <--->     LEU222    A        O_1635     3.5    35.1
   @>     LYS166    A        NZ_699  <--->     ALA203    A        O_1317     3.5    31.4
   @>     ASN386    A      ND2_4298  <--->     TRP382    A      NE1_4233     3.5    40.0
   @>     GLY145    A         N_319  <--->     GLY142    A         O_276     3.5    32.9
   @>     LYS250    A        N_2078  <--->     TYR246    A        O_2023     3.5    14.1
   @> Number of detected hydrogen bonds: 54.
   @> Creating file with dummy atoms
   @> struc_homologs_BLAST/align__addH_3dj5A.pdb
   @> 4154 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> Calculating hydrogen bonds.
   @>      DONOR (res chid atom)   <--->       ACCEPTOR (res chid atom)    Distance  Angle
   @>     PHE241    A        N_1708  <--->     PHE335    A        O_3061     2.8    39.0
   @>     LEU162    A         N_337  <--->     ARG150    A         O_149     2.8    30.9
   @>     LEU172    A         N_518  <--->     ALA163    A         O_361     2.9    30.2
   @>     LYS179    A         N_643  <--->     THR217    A        O_1293     2.9    21.4
   @>     ASP269    A        N_2168  <--->     ASP287    A      OD2_2456     2.9    28.9
   @>     PHE335    A        N_3056  <--->     LEU331    A        O_2995     2.9    34.2
   @>     ILE222    A        N_1382  <--->     GLY211    A        O_1207     2.9    31.4
   @>     ARG384    A       NE_3875  <--->     LYS378    A        O_3774     2.9    28.1
   @>     ARG208    A        N_1137  <--->     GLU224    A      OE2_1434     2.9    11.0
   @>     GLN236    A        N_1615  <--->     TYR232    A        O_1541     2.9     8.0
   @>     GLN190    A       NE2_831  <--->     GLU194    A       OE2_915     2.9    28.4
   @>     LEU238    A        N_1654  <--->     GLU234    A        O_1586     2.9    25.0
   @>     CYS260    A        N_2007  <--->     ALA256    A        O_1951     2.9    34.8
   @>     ASP371    A        N_3646  <--->     GLU367    A        O_3595     3.0    35.9
   @>     ARG356    A        N_3405  <--->     ARG352    A        O_3332     3.0    36.5
   @>     GLU165    A         N_390  <--->     PHE170    A         O_484     3.0    39.4
   @>     GLU196    A         N_932  <--->     ARG192    A         O_858     3.0    20.9
   @>     ILE373    A        N_3677  <--->     ALA369    A        O_3617     3.0    37.6
   @>     LYS263    A        N_2046  <--->     TYR259    A        O_1991     3.0    17.9
   @>     ARG352    A        N_3327  <--->     GLN348    A        O_3265     3.0     2.1
   @>     VAL330    A        N_2974  <--->     TRP326    A        O_2918     3.0    37.0
   @>     ARG233    A        N_1557  <--->     THR230    A      OG1_1514     3.0    21.3
   @>     TYR308    A        N_2618  <--->     THR305    A        O_2578     3.0     3.8
   @>     ALA185    A         N_752  <--->     GLN181    A         O_684     3.0    35.9
   @>     ARG192    A         N_853  <--->     GLU188    A         O_790     3.0    36.0
   @>     HIS261    A        N_2018  <--->     LEU257    A        O_1961     3.1    32.6
   @>     ASP324    A        N_2882  <--->     ASP320    A        O_2822     3.1     2.2
   @>     GLU252    A        N_1888  <--->     THR248    A        O_1825     3.1    17.2
   @>     SER262    A        N_2035  <--->     SER258    A        O_1980     3.1    19.0
   @>     LEU182    A         N_696  <--->     PHE178    A         O_628     3.1    13.4
   @>     GLY158    A         N_279  <--->     GLY155    A         O_236     3.1    33.6
   @>     ARG375    A        N_3707  <--->     ASP371    A        O_3651     3.1    11.7
   @>     GLU194    A         N_901  <--->     GLN190    A         O_822     3.1    20.1
   @>     VAL323    A        N_2866  <--->     ASP320    A      OD1_2827     3.1    12.3
   @>     TYR333    A        N_3020  <--->     GLY329    A        O_2973     3.1    36.7
   @>     GLU349    A        N_3277  <--->     THR346    A      OG1_3233     3.1    22.2
   @>     LEU325    A        N_2894  <--->     GLU321    A        O_2834     3.1    15.2
   @>     ARG193    A         N_877  <--->     HIS189    A         O_805     3.1    33.8
   @>     SER327    A        N_2937  <--->     VAL323    A        O_2871     3.1    22.3
   @>     HIS393    A        N_4011  <--->     VAL390    A        O_3966     3.2    26.4
   @>     LEU372    A        N_3658  <--->     GLY368    A        O_3611     3.2    26.4
   @>     THR350    A        N_3292  <--->     THR346    A        O_3230     3.2    18.0
   @>     ILE197    A         N_947  <--->     ARG193    A         O_882     3.2    33.7
   @>     GLY329    A        N_2967  <--->     LEU325    A        O_2899     3.2    22.1
   @>     LEU376    A        N_3731  <--->     LEU372    A        O_3663     3.2    39.7
   @>     GLY186    A         N_762  <--->     GLU183    A         O_720     3.2    16.4
   @>     VAL187    A         N_769  <--->     LEU182    A         O_701     3.2     9.9
   @>     TYR232    A        N_1536  <--->     PRO272    A        O_2225     3.2    26.8
   @>     LEU336    A        N_3076  <--->     CYS332    A        O_3014     3.2    32.0
   @>     GLU389    A        N_3946  <--->     THR386    A      OG1_3911     3.2    35.4
   @>     ASN399    A        N_4117  <--->     ILE396    A        O_4071     3.3    27.4
   @>     GLU392    A        N_3996  <--->     ALA388    A        O_3941     3.4    17.0
   @>     ILE314    A        N_2718  <--->     PRO310    A        O_2662     3.4    20.0
   @>     LYS179    A        NZ_661  <--->     HIS214    A      ND1_1259     3.5    35.5
   @>     ARG384    A      NH1_3878  <--->     HIS379    A        O_3796     3.5    11.1
   @> Number of detected hydrogen bonds: 55.
   @> Creating file with dummy atoms
   @> struc_homologs_BLAST/align__addH_2x6dA.pdb
   @> 4202 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> Calculating hydrogen bonds.
   @>      DONOR (res chid atom)   <--->       ACCEPTOR (res chid atom)    Distance  Angle
   @>     PHE322    A        N_3096  <--->     LEU318    A        O_3035     2.8    15.1
   @>     LEU240    A        N_1880  <--->     TYR236    A        O_1816     2.8    33.3
   @>     SER361    A        N_3736  <--->     ARG357    A        O_3667     2.9    30.2
   @>     SER342    A        N_3434  <--->     TYR338    A        O_3353     2.9    39.1
   @>     ARG343    A        N_3445  <--->     LYS339    A        O_3374     2.9    38.2
   @>     GLN223    A        N_1594  <--->     TYR219    A        O_1520     2.9     7.0
   @>     THR238    A        N_1851  <--->     ALA234    A        O_1792     2.9    38.2
   @>     THR235    A        N_1797  <--->     GLN231    A        O_1737     2.9     9.8
   @>     CYS319    A        N_3049  <--->     LEU315    A        O_2993     2.9    37.4
   @>     LEU315    A        N_2988  <--->     ASP311    A        O_2927     2.9     9.6
   @>     ASN242    A        N_1909  <--->     THR238    A        O_1856     3.0    39.2
   @>     ARG362    A        N_3747  <--->     ASP358    A        O_3691     3.0    28.4
   @>     TYR246    A        N_1963  <--->     ASN242    A        O_1914     3.0    36.6
   @>     SER283    A        N_2549  <--->     HIS306    A        O_2845     3.0    34.3
   @>     CYS247    A        N_1984  <--->     ALA243    A        O_1928     3.0    12.8
   @>     TYR320    A        N_3060  <--->     GLY316    A        O_3013     3.0    35.7
   @>     TYR338    A        N_3348  <--->     TYR334    A        O_3286     3.0    34.3
   @>     SER186    A         N_962  <--->     VAL182    A         O_900     3.0    28.2
   @>     GLU379    A        N_4052  <--->     ARG375    A        O_3983     3.0    32.3
   @>     VAL344    A        N_3469  <--->     ILE341    A        O_3420     3.0    37.7
   @>     ILE209    A        N_1361  <--->     GLY198    A        O_1186     3.0    29.7
   @>     HIS248    A      NE2_2010  <--->     ASP311    A      OD2_2933     3.1    23.7
   @>     ASP256    A        N_2145  <--->     HIS254    A      ND1_2114     3.1    36.7
   @>     LYS271    A       NZ_2385  <--->     PRO191    A        O_1054     3.1    32.9
   @>     VAL377    A        N_4017  <--->     MET373    A        O_3947     3.1    27.4
   @>     LEU323    A        N_3116  <--->     CYS319    A        O_3054     3.1    16.5
   @>     HIS380    A        N_4067  <--->     VAL377    A        O_4022     3.2    21.1
   @>     THR233    A        N_1773  <--->     ASP229    A        O_1710     3.2    32.3
   @>     GLU336    A        N_3319  <--->     THR333    A        O_3272     3.2    26.9
   @>     HIS306    A        N_2840  <--->     SER283    A        O_2554     3.2    18.0
   @>     TRP313    A        N_2953  <--->     LYS309    A        O_2889     3.2    31.5
   @>     ALA356    A        N_3652  <--->     THR353    A      OG1_3624     3.2    16.9
   @>     ARG371    A        N_3904  <--->     PRO368    A        O_3866     3.2    19.6
   @>     VAL182    A         N_895  <--->     ARG179    A         O_837     3.2    31.6
   @>     ARG371    A      NH1_3922  <--->     HIS366    A        O_3836     3.2    39.7
   @>     SER314    A        N_2977  <--->     VAL310    A        O_2911     3.3    19.9
   @>     GLU134    A         N_143  <--->     ARG151    A         O_411     3.3    17.4
   @>     LYS224    A       NZ_1629  <--->     GLU221    A      OE1_1573     3.3    25.2
   @>     LYS250    A        N_2023  <--->     TYR246    A        O_1968     3.3    38.1
   @>     ARG179    A         N_832  <--->     GLU175    A         O_769     3.3    25.3
   @>     LEU359    A        N_3698  <--->     GLY355    A        O_3651     3.3     9.3
   @>     LEU169    A         N_730  <--->     PHE165    A         O_666     3.4     6.5
   @>     TYR219    A        N_1515  <--->     PRO259    A        O_2202     3.4    26.8
   @>     THR384    A        N_4141  <--->     PRO381    A        O_4088     3.4    39.9
   @>     GLY316    A        N_3007  <--->     LEU312    A        O_2939     3.4    18.4
   @>     LYS339    A        N_3369  <--->     GLN335    A        O_3307     3.4    18.9
   @>     HIS254    A      ND1_2114  <--->     ILE257    A        N_2157     3.4    35.9
   @>     LEU270    A        N_2348  <--->     GLU239    A      OE1_1878     3.4    36.6
   @>     LYS224    A        N_1611  <--->     GLU221    A        O_1565     3.4    32.2
   @>     LYS271    A       NZ_2385  <--->     GLU211    A      OE2_1413     3.4    35.6
   @>     ARG304    A      NH2_2820  <--->     HIS366    A      NE2_3847     3.4    27.2
   @>     HIS248    A        N_1995  <--->     SER245    A        O_1957     3.5    38.5
   @>     VAL324    A        N_3135  <--->     GLU321    A        O_3086     3.5    32.7
   @>     PHE157    A         N_517  <--->     GLU152    A         O_435     3.5    33.7
   @>     LEU378    A        N_4033  <--->     ARG375    A        O_3983     3.5    37.4
   @> Number of detected hydrogen bonds: 55.
   @> Creating file with dummy atoms
   @> struc_homologs_BLAST/align__addH_2c6eA.pdb
   @> 4141 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> Calculating hydrogen bonds.
   @>      DONOR (res chid atom)   <--->       ACCEPTOR (res chid atom)    Distance  Angle
   @>     LEU243    A        N_1964  <--->     LEU239    A        O_1916     2.7    38.9
   @>     GLN369    A        N_3804  <--->     ASN366    A        O_3770     2.8    29.7
   @>     LEU221    A        N_1606  <--->     VAL217    A        O_1535     2.9    32.6
   @>     ILE208    A        N_1392  <--->     GLY197    A        O_1217     2.9    16.0
   @>     TRP381    A      NE1_4028  <--->     GLU229    A      OE2_1762     2.9    33.1
   @>     ARG194    A        N_1147  <--->     GLU210    A      OE2_1444     2.9    40.0
   @>     VAL217    A        N_1530  <--->     LEU261    A        O_2277     2.9    11.6
   @>     TYR337    A        N_3265  <--->     TYR333    A        O_3203     2.9    33.3
   @>     ALA233    A        N_1818  <--->     GLU229    A        O_1753     2.9    39.6
   @>     LYS338    A        N_3286  <--->     GLN334    A        O_3224     3.0    25.4
   @>     ILE183    A         N_957  <--->     ARG179    A         O_892     3.0    25.3
   @>     GLU180    A         N_911  <--->     GLN176    A         O_832     3.0    14.5
   @>     VAL251    A        N_2100  <--->     CYS246    A        O_2020     3.0    13.1
   @>     ARG356    A        N_3579  <--->     THR352    A        O_3538     3.0    20.1
   @>     ALA149    A         N_372  <--->     LEU158    A         O_537     3.0    38.4
   @>     SER313    A        N_2894  <--->     VAL309    A        O_2828     3.0    38.0
   @>     ARG342    A        N_3362  <--->     LYS338    A        O_3291     3.0    22.7
   @>     ILE340    A        N_3332  <--->     THR336    A        O_3256     3.0    31.2
   @>     CYS318    A        N_2966  <--->     LEU314    A        O_2910     3.0    23.8
   @>     LEU158    A         N_532  <--->     ALA149    A         O_377     3.0    33.9
   @>     LEU168    A         N_706  <--->     PHE164    A         O_642     3.0    17.4
   @>     ALA384    A        N_4072  <--->     PRO380    A        O_4005     3.1    16.4
   @>     LYS161    A         N_580  <--->     LEU207    A        O_1378     3.1    37.4
   @>     SER341    A        N_3351  <--->     TYR337    A        O_3270     3.1    21.0
   @>     ARG219    A        N_1567  <--->     THR216    A      OG1_1524     3.1    23.9
   @>     LEU239    A        N_1911  <--->     TYR235    A        O_1847     3.1    29.2
   @>     VAL316    A        N_2931  <--->     TRP312    A        O_2875     3.1    33.5
   @>     CYS246    A        N_2015  <--->     ALA242    A        O_1959     3.1    33.9
   @>     THR232    A        N_1804  <--->     ASP228    A        O_1741     3.1    18.4
   @>     LEU362    A        N_3688  <--->     LEU358    A        O_3620     3.1    37.1
   @>     ILE359    A        N_3634  <--->     ALA355    A        O_3574     3.1    28.2
   @>     SER387    A        N_4107  <--->     THR234    A      OG1_1836     3.1    23.1
   @>     HIS379    A        N_3984  <--->     VAL376    A        O_3939     3.1     7.8
   @>     LEU224    A        N_1664  <--->     GLU220    A        O_1596     3.2     3.5
   @>     THR336    A        N_3251  <--->     THR332    A        O_3189     3.2    12.0
   @>     GLU238    A        N_1896  <--->     THR234    A        O_1833     3.2    18.5
   @>     ARG361    A        N_3664  <--->     ASP357    A        O_3608     3.2     8.1
   @>     ASN241    A        N_1940  <--->     THR237    A        O_1887     3.2    20.2
   @>     ASN260    A      ND2_2269  <--->     ASP255    A      OD1_2186     3.2    33.1
   @>     GLY315    A        N_2924  <--->     LEU311    A        O_2856     3.2    28.2
   @>     LEU358    A        N_3615  <--->     GLY354    A        O_3568     3.2    27.7
   @>     TYR319    A        N_2977  <--->     GLY315    A        O_2930     3.2    17.1
   @>     SER248    A        N_2043  <--->     SER244    A        O_1988     3.2    15.8
   @>     LEU322    A        N_3033  <--->     CYS318    A        O_2971     3.2    27.9
   @>     THR237    A        N_1882  <--->     ALA233    A        O_1823     3.2    24.3
   @>     LEU214    A        N_1490  <--->     LEU263    A        O_2315     3.2    33.9
   @>     VAL146    A         N_316  <--->     GLY139    A         O_223     3.3    31.4
   @>     ASN385    A        N_4082  <--->     TRP381    A        O_4020     3.3    37.2
   @>     GLU378    A        N_3969  <--->     ARG374    A        O_3900     3.3    27.3
   @>     GLU335    A        N_3236  <--->     THR332    A      OG1_3192     3.3    39.6
   @>     LYS249    A        N_2054  <--->     TYR245    A        O_1999     3.3    26.4
   @>     ALA128    A          N_43  <--->     ASP131    A        OD1_97     3.4     4.6
   @>     GLU220    A        N_1591  <--->     THR216    A        O_1521     3.4    37.4
   @>     VAL309    A        N_2823  <--->     ASP306    A        O_2779     3.4    19.5
   @>     TYR218    A        N_1546  <--->     PRO258    A        O_2233     3.4    28.7
   @>     GLY197    A        N_1211  <--->     ILE208    A        O_1397     3.4    20.4
   @>     ARG254    A        N_2152  <--->     ASP310    A      OD1_2849     3.4    20.2
   @>     ASP310    A        N_2839  <--->     GLU307    A        O_2791     3.4    27.3
   @>     ALA171    A         N_762  <--->     GLN167    A         O_694     3.4    23.6
   @>     VAL181    A         N_926  <--->     ARG178    A         O_868     3.4    34.7
   @>     GLY215    A        N_1509  <--->     ALA212    A        O_1471     3.4     6.4
   @>     TRP312    A        N_2870  <--->     LYS308    A        O_2806     3.5    22.5
   @>     GLN184    A       NE2_990  <--->     ASP273    A      OD1_2459     3.5    17.3
   @>     VAL376    A        N_3934  <--->     MET372    A        O_3864     3.5    26.9
   @> Number of detected hydrogen bonds: 64.
   @> Creating file with dummy atoms
   ..
   ..


Now, we need to go to the newly created folder struc_homologs_BLAST and use
:func:`.findClusterCenters` for each type of interactions to detect clusters
of interactions. We cannot do it in an automatic way becasue each system is
different and often default parameters of *numC* and *distC* should be tuned.

To compute fingerprint interactions of hydrogen bonds use the command below.
This function uses the prefix 'INT_HBs_' to select analyzed interactions.

.. ipython:: python
   :verbatim:

   findClusterCenters('INT_HBs_*.pdb', selection = 'resname DUM')

.. parsed-literal::

   @> 4362 atoms and 1 coordinate set(s) were parsed in 0.11s.
   @> 4433 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4122 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4383 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4568 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4594 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4269 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4249 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4205 atoms and 1 coordinate set(s) were parsed in 0.05s.
   ..
   ..
   @> 4177 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4428 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4245 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4269 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4417 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4610 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> Results are saved in clusters_INT_HBs_.pdb.


To compute fingerprint interactions of salt bridges use the prefix
'INT_SBs_':

.. ipython:: python
   :verbatim:

   findClusterCenters('INT_SBs_*.pdb', selection = 'resname DUM')

.. parsed-literal::

   @> 4306 atoms and 1 coordinate set(s) were parsed in 0.11s.
   @> 4383 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4082 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4337 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4532 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4558 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4228 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4213 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4154 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4303 atoms and 1 coordinate set(s) were parsed in 0.06s.
   @> 4350 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4107 atoms and 1 coordinate set(s) were parsed in 0.05s.
   ..
   ..
   @> 4199 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4233 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4369 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4569 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> Results are saved in clusters_INT_SBs_.pdb.


To compute fingerprint interactions of repulsive ionic bonding use the prefix
'INT_RIB_':

.. ipython:: python
   :verbatim:

   findClusterCenters('INT_RIB_*.pdb', selection = 'resname DUM')

.. parsed-literal::

   @> 4367 atoms and 1 coordinate set(s) were parsed in 0.11s.
   @> 4322 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4510 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4537 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4218 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4199 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4142 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4286 atoms and 1 coordinate set(s) were parsed in 0.05s.
   ..
   ..
   @> 4353 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4182 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4215 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4351 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4548 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> Results are saved in clusters_INT_RIB_.pdb.


To compute fingerprint of pi-stacking interactions use the prefix
'INT_PiStack_':

.. ipython:: python
   :verbatim:

   findClusterCenters('INT_PiStack_*.pdb', selection = 'resname DUM')

.. parsed-literal::

   @> 4292 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4368 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4066 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4322 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4510 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4536 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4220 atoms and 1 coordinate set(s) were parsed in 0.04s.
   ..
   ..
   @> 4119 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4354 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4182 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4352 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4548 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> Results are saved in clusters_INT_PiStack_.pdb.


To compute fingerprint of pi-stacking interactions use the prefix
'INT_PiCat_':

.. ipython:: python
   :verbatim:

   findClusterCenters('INT_PiCat_*.pdb', selection = 'resname DUM')

.. parsed-literal::

   @> 4295 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4369 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4068 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4323 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4516 atoms and 1 coordinate set(s) were parsed in 0.05s.
   ..
   ..
   @> 4184 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4219 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4356 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4551 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> Results are saved in clusters_INT_PiCat_.pdb.


To compute fingerprint of hydrophobic interactions use the prefix
'INT_HPh_':

.. ipython:: python
   :verbatim:

   findClusterCenters('INT_HPh_*.pdb', selection = 'resname DUM')

.. parsed-literal::

   @> 4374 atoms and 1 coordinate set(s) were parsed in 0.12s.
   @> 4446 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4142 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4398 atoms and 1 coordinate set(s) were parsed in 0.05s.
   @> 4594 atoms and 1 coordinate set(s) were parsed in 0.05s.
   ..
   ..
   @> 4427 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4260 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4295 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4426 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> 4625 atoms and 1 coordinate set(s) were parsed in 0.04s.
   @> Results are saved in clusters_INT_HPh_.pdb.

We can further visualize those results in a graphical program like VMD_.
The obtained results are with default numC and distC parameters, and it means
that we were identifying interactions that were within 0.3 Ang. from each other 
(distC) in at least three structures (numC). To see more information about those 
parameters, see the WatFinder tutorial_.

Visualization of hydrogen bonds clusters is the following:

.. figure:: images/blast_hbs.png
   :scale: 60 %

   
Visualization of salt bridges clusters is the following:

.. figure:: images/blast_sbs.png
   :scale: 60 %
   
   
Visualization of repulsive ionic bonding clusters is the following:

.. figure:: images/blast_rib.png
   :scale: 60 %
   
   
Visualization of pi-cation clusters is the following:

.. figure:: images/blast_picat.png
   :scale: 60 %
   
   
Visualization of pi-stacking clusters is the following:

.. figure:: images/blast_pistack.png
   :scale: 60 %


Visualization of hydrophobic interactions clusters is the following:
   
.. figure:: images/blast_hph.png
   :scale: 60 %


.. _WatFinder tutorial: http://www.bahargroup.org/prody/tutorials/watfinder_tutorial
