Main scripts/functions:

- run\_manopt\_rsom\_genproc.m\
Main executable that runs the simulated experiments on a graph generated through testNetwork

- run\_som\_adjmat.m\
Main executable that runs the simulated experiments on a user-defined graph (need for its adjacency matrix)

- do\_rsom\_procrustes\_genproc.m\
Function that does the main job called in run\_manopt\_rsom\_genproc.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

- align2d.m\
Function that enables recovery of the 2-degree node solutions after performing RS

- check\_Rb\_ambiguity.m\
Checking 2-deg node recovery procedure step-by-step (all the formulas leading to that)

- check\_recovery\_error.m\
Checking 2-deg node recovery error on a 6 node graph case

- check\_Tij1j2\_init.m\
Checking whether the Tij1j2 and Tij1j2_tilde vectors are being formed accordingly (regards the 2-deg node recovery procedure)

- compute\_hrt.m\
Computes the H_{RT} part of the Hessian

- compute\_htr.m\
Computes the H_{TR} part of the Hessian

- cost\_genproc.m\
Calls rsom\_cost\_base

- edge\_diffs\_2\_T.m\
Recover "global" translations from the edge-wise differences; technically serves as the inverse of make\_T\_edges

- eigencheck\_hessian\_genproc.m\
Checks if the Eigenvalue-Eigenvector couple is correct for the genproc Hessian

- find\_Rb\_initguess.m\
Find initial guess for Rb recovery: DEPRECATED and substituted with recoverRitilde.m

- grad\_genproc.m\
Compute the generalized procrustes gradient

- hess\_genproc.m\
Compute the generalized procrustes hessian

- hess\_genproc\_R.m\
Compute the rotation part of the generalized procrustes hessian

- hess\_genproc\_shifted.m\
Compute the shifted generalized procrustes hessian

- hess\_genproc\_R.m\
Compute the translation part of the generalized procrustes hessian

- make\_p\_last.m\
Return matrix that extracts the last p-d rows of matrix A when multiplied on left side of A i.e., P_last_out * A

- make\_T\_edges.m\
Make edge-wise translation differences vector starting from a translation "pose" vector; technically serves as the inverse of edge\_diffs\_2\_T

- make\_Tij1j2s.m\
Make Tij1j2 and Tij1j2_tilde vectors

- POCdegree2Align.m\
Original/base/raw version of recoverRitilde.m

- POCRotateToMinimizeLastEntries.m\
Finds the Qx matrix such that 3+ degree nodes are correctly recovered after RS

- pim\_function\_genproc.m\
Power iteration method able to run on a struct hessian that involves both
a .R and a .T component

- procrustes\_Rb.m\
Part of the recoverRitilde() pipeline

- Qc\_recovery\_Rb\_initguess.m\
DEPRECATED Part of the old 2-deg node recovery

- qcrb\_geodFun.m\
Generate geodesics for a 2x2 and 4x4 orthogonal matrix pair

- qcrb\_geodFun.m\
Generate geodesics starting from a random point for a 2x2 and 4x4 orthogonal matrix pair

- recoverRitilde.m\
2-deg node recovery procedure to recover the rotation equivalent to the RS Stiefel output

- rsom\_cost\_base.m\
Compute the cost function with the basic/original formulation: 
\sum_{\ijE}\norm{R_i\iframe{}T_{ij}-T_j+T_i}^2 

- rsom\_genproc.m\
Run the RSOM pipeline with the optimization based on the generalized
Procrustes product manifold

- rsom\_pim\_hessian\_genproc.m\
Power iteration method for a generalized Procrustes Hessian 
that involves both a .R and a .T component

- solve\_gauge\_symmetry\_1rot.m\
Script containing the pipeline able to solve the "high"-degree recovery -> example running on a single rotation

- solve\_gauge\_symmetry.m\
Script containing the pipeline able to solve the "high"-degree recovery -> example running on N=6 rotations

- test\_cost\_grad\_hess\_genproc\_manopt.m\
Runs Manopt's checkgradient() and checkhessian() tests on the genproc cost

- test\_fail\_recovery\_zerocost.m\
Testing whether fail recoveries start from zero-cost (i.e., RS has succesfully stopped in a global minimum)

- test\_full\_recovery.m\
Test that runs the recovery procedure on a sample case with both low and high deg nodes

- test\_make\_T\_edges\_2ways.m\
Test whether the functions for doing the T -> T\_edges (T\_diffs) -> T transition function correctly

- test\_POCdegree2Align\_gt.m\
Test POCdegree2Align procedure directly on ground truth (padded) rotations and translations

- test\_POCdegree2Align\_som.m\
Test POCdegree2Align procedure directly on a case coming from an actual RSOM output

- test\_POCdegree2Align.m\
Test POCdegree2Align function

- test\_POCRotateToMinimizeLastEntries.m\
Test recovery procedure for 3+ degree nodes, all translation differences

- test\_recovery\_recap.m\
Testing one by one the formulas leading to the 2-deg node recovery procedure

- test\_transl\_recovery.m\
Testing whether the translation recovery (based on edge differences) actually works

- testNetwork\_adj.m\
Version of testNetwork that generates the network given an adjacency matrix as input

- testNetwork\_params.m\
Version of testNetwork that is able to take as input the following parameters: number of nodes, mode (full/banded graph), max node degree

- testNetwork\_params.csv, noise\_test\_params.csv, mus.txt, sigmas.txt\
Test params files

- remake\_plots.m\
Self-explanatory, for paper