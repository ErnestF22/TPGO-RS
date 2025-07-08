# TPGO-RS
Translation-Only Pose Graph Optimization through a Riemannian Staircase

MATLAB Code for Shape of Motion project.

PDF of Report "Hessians on Stiefel for a bilinear cost function" -> [link](https://drive.google.com/file/d/1jg5BSRPsQMLih3ln2P9feqrpJWpfVzHv/view?usp=share_link)

Dependencies:
- Manopt -> [link](https://www.manopt.org/)
- CVX -> [link](http://cvxr.com/)

\
Brief description of the scripts/functions present in this folder (check subfolders README for more infos on those):

- run\_manopt\_rsom\_genproc.m\
Main executable that runs the simulated experiments on a graph generated through testNetwork

- run\_som\_adjmat.m\
Main executable that runs the simulated experiments on a user-defined graph (need for its adjacency matrix)

- run\_final\_tests.m\
Script that runs the most relevant tests

- do\_rsom\_procrustes\_genproc.m\
Function that does the main job called in run\_manopt\_rsom\_genproc.m

- remake\_plots.m\
Self-explanatory, for paper

- rsom\_genproc.m\
Run the TPGO-RS pipeline with the optimization based on the generalized
Procrustes product manifold, and the Riemannian Staircase adaptation

