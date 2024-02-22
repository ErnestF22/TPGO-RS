- run\_manopt\_rsom\_genproc.m\\
Main executable that runs the simulated experiments.

- do\_rsom\_procrustes\_genproc.m\\
Function that does the main job called in run\_manopt\_rsom\_genproc.m\

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

- compute\_hrt.m\\
Computes the H_{RT} part of the Hessian

- compute\_htr.m\\
Computes the H_{TR} part of the Hessian

- cost\_genproc.m\\
Calls rsom\_cost\_base

- eigencheck\_hessian\_genproc.m\\
Checks if the Eigenvalue-Eigenvector couple is correct for the genproc Hessian

- grad\_genproc.m\\
Compute the generalized procrustes gradient

- hess\_genproc.m\\
Compute the generalized procrustes hessian

- hess\_genproc\_R.m\\
Compute the rotation part of the generalized procrustes hessian

- hess\_genproc\_shifted.m\\
Compute the shifted generalized procrustes hessian

- hess\_genproc\_R.m\\
Compute the translation part of the generalized procrustes hessian

- pim\_function\_genproc.m\\
Power iteration method able to run on a struct hessian that involves both
a .R and a .T component

- rsom\_cost\_base.m\\
Compute the cost function with the basic/original formulation: 
\sum_{\ijE}\norm{R_i\iframe{}T_{ij}-T_j+T_i}^2 

- rsom\_genproc.m\\
Run the RSOM pipeline with the optimization based on the generalized
Procrustes product manifold

- rsom\_pim\_hessian\_genproc.m\\
Power iteration method for a generalized Procrustes Hessian 
that involves both a .R and a .T component

- test\_cost\_grad\_hess\_genproc\_manopt.m\\
Runs Manopt's checkgradient() and checkhessian() tests on the genproc cost




