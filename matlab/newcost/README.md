- do\_rsom\_procrustes\_rs.m\
Function that runs an instance of the run\_manopt\_rsom.m tests

- rsom\_manopt.m\
Run the Manopt-based step1 RSOM optimization (composed of step1, step2)

- rsom\_pim\_hessian.m\
Function that runs the power iteration method for the step1 RSOM Hessian

- rsom\_rs.m\
Run the RSOM Riemannian Staircase without genproc

- rsom\_step1.m\
Run the Manopt-based step1 RSOM optimization (rotation estimation)

- rsom\_step2.m\
Run the Manopt-based step1 RSOM optimization (translation estimation)

- run\_manopt\_rsom.m\
Run the repetitive simulated tests comparing Manopt, ICP-Manopt and Procrustes

- save\_data\_script.m\
Useful script to save matrices/terms that appear in the cost function 
when stopping test execution at a desired breakpoint (P, frct, LR, PR, BR);
supports also saving of the "next" terms, with next referring to the
following staircase step

- test\_rsom\_pim\_hessian.m\
Test the RSOM PIM Hessian functioning on given cost terms/matrices (that
are being loaded from the data/ folder)

- test\_rsom\_pim\_hessian\_fun.m\
Test the RSOM PIM Hessian functioning on given cost terms/matrices (that
are being loaded from the data/ folder), using a direct call to the
rsom\_pim\_hessian function



