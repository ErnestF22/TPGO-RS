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

- find\_Rb\_initguess.m\
Find initial guess for Rb recovery: DEPRECATED and substituted with recoverRitilde.m

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