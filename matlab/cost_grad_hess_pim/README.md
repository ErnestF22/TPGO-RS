- compute\_hrt.m\
Compute one of the cross terms of the Hessian

- compute\_htr.m\
Compute the other cross term of the Hessian

- cost\_genproc.m\
Return cost for function defined on genproc manifold

- eigencheck\_hessian\_genproc.m\
Checks if the Eigenvalue-Eigenvector couple is correct for the genproc Hessian

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

- pim\_function\_genproc.m\
Power iteration method able to run on a struct hessian that involves both
a .R and a .T component

- rsom\_cost\_base.m\
Compute the cost function with the basic/original formulation: 
\sum_{\ijE}\norm{R_i\iframe{}T_{ij}-T_j+T_i}^2 

- rsom\_pim\_hessian\_genproc.m\
Power iteration method for a generalized Procrustes Hessian 
that involves both a .R and a .T component