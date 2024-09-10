# som
Shape of Motion
MATLAB Code for Shape of Motion project.

PDF of Report "Hessians on Stiefel for a bilinear cost function" -> [link](https://drive.google.com/file/d/1jg5BSRPsQMLih3ln2P9feqrpJWpfVzHv/view?usp=share_link)

Dependencies:
- Manopt -> [link](https://www.manopt.org/)
- CVX -> [link](http://cvxr.com/)

\
Brief description of the scripts/functions present in this folder (check subfolders README for more infos on those):\
- check\_prev\_cost\_script.m\
Used to check whether RS padding impacts cost (it should not)

- compare\_procrustes\_manopt\_genproc.m\
Compare Procrustes VS Manopt-ICP VS Manopt Genproc (SOM); also called "cpm" in the following

- compare\_procrustes\_manopt\_genproc.m\
Compare Procrustes VS Manopt-ICP VS Manopt SOM (iterative RS, no Generalized Procrustes) 

- compare\_procrustes\_manopt\_genproc.m\
Compare Procrustes VS Manopt-ICP

- cp2m\_noiseinit.m\
Compare Procrustes VS the two Manopt versions, with noisy initial guess

- cp2m\_test\_end\_workspace.mat\
WS saved at the end of running cp2m\_noiseinit.m

The do\_[...] functions are called from the compare\_[...] scripts and do the whole procedure with the same [...], while the compare\_[...] scripts setup the test and gather results.\
"Two modes" referring to gradients or Hessians refer to manual VS automated (through Manopt library AD) compatation.\

- make\_cost\_matrix.m\
Used to understand SE-Sync cost matrix formulation

- make\_input\_subset.m\
Make smaller subsets of the input graph, counting only 2, 3, ... nodes

- manopt\_first\_tutorial.m\
Self explanatory, copied from Manopt library documentation

- old\_compare\_procrustes\_manopt.m\
First (old, DEPRECATED) version of cpm

- plot\_results.m\
Function that aids result plotting (for the different versions of the tests)

- sesync\_simple\_stiefel\_onlyrot.m\
Running RTR on a simple, rotation-only formulation of SE-Sync's cost

- ShapeOfMotion.m\
Draft class for running things C++ style

- som\_testnet.m\
Example of calling testNetwork() function -> problem case generator 

- testNetwork\_params.csv, noise\_test\_params.csv, mus.txt, sigmas.txt
Test params files
