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