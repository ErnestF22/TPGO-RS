From "PROGRESSIVE TESTS (CODEMETA’S TESTNETWORK INPUTS)" section on project's overleaf report

- test\_procrustes\_dummy.m\
Test that contains a trivial set/setup of points that can be translated and rotated ad libitum in order to test the correctness of any implemented Procrustes method (i.e. Kabsch, Umeyama, SoM, ...)

- test\_rotation\_cost.m\
The cost function of step 1, when rotations are given and equal to ground truth, must be equal to 0

- test\_translation\_cost.m\ 
The cost function of step 2, when translation are given and equal to ground truth, must be equal to 0.\
Without noise, and with known rotations (e.g. setting all rotations equal to gt, with all rotations equal to $I$ in the first try), $A T - B = \mathbf{0}$ must be true, with $A$, $B$ as in \eqref{eq:Lij_step1} and \eqref{eq:Pij_step1} and $T$ filled with ground truth translation values.

- test\_make\_LT\_PT.m\
Form L(T) and P(T) with two separate methods (matricial and loop) and check if they come out as equal.   

- test\_make\_A\_b.m\
Equivalent of $test\\_make\\_LT\\_PT.m$

- test\_make\_M\_N.m\
Equivalent of $test\\_make\\_LT\\_PT.m$

- funCheckDer\_test.m\
Test taken from codemeta, to be adapted and reused for testing gradients

- Manopt gradients \
Test correctness of translation estimation gradients hand-written expression, using an adapted version of codemeta's funCheckDer \
In particular, when it comes to gradients, we test:
    - Euclidean gradient wrt R\
    (in "test\_manopt\_step1\_gradients.m")
    - Riemannian gradient wrt R\
    (in "test\_manopt\_step1\_gradients.m")
    - Euclidean gradient wrt T\
    (in "test\_manopt\_step2\_gradients.m")
Without noise, they must be equal to 0\
In code, $g$ indicates Euclidean gradient, $h$ indicates Riemannian gradient.

- test\_procrustes\_fullgraph\_groundtruth.m\
Test SoM Procrustes considering that all data given is correct and that the graph is full (i.e. all edges are available)

- test\_cpm\_single\_nonoise.m\
Compare Procrustes and Manopt on a single run of a test without noise\
"Note: cpm is an acronym for Compare Procrustes and Manopt}

- test\_cpm\_single.m\
Compare Procrustes and Manopt on a single run of a noise test (fixed low $\sigma$ i.e. $\sigma = 0.1$, null $\mu$ "added" inside Codemeta's testNetwork)

- test\_manopt\_sep\_gen.m\
Compare transf out results obtained with the two Manopt pipelines: separated, generalized

- test\_procrustes\_randinit.m\ 
Test created to verify Procrustes pipeline performance when dealing with increasingly noisy inputs (which should lead the estimation to fail)

- compare\_manopt\_maxiter.m\
Test the improvement obtained on Manopt results by increasing the maximum number of ICP iterations

- compare\_manopt\_riemgrad.m\
Compare the results obtained through Manopt when using auto-computed Riemannian gradient vs the hand-computed one; the comparison includes also execution times; graphical results from this test can be found in figure \ref{fig:manopt_analytical_vs_ad_grad}

- compare\_manopt\_hessian\_ad\_analytic.m\
Test Manopt results and execution times obtained when we use Manopt with automatic differentiation for computing Hessian, VS analytically pre-computed Hessian; graphical results from this test can be found in figure \ref{fig:manopt_analytical_vs_ad_hessian}

- test\_make\_input\_subset.m\
Test that checks if the "make\_input\_subset()" function works. Said function returns a partial input of a full testdata input generated through the testNetwork methods; it returns all data associated to the $cameras\_ids$ passed as input, extracting it by
removing the data associated to other cameras.

- test\_manopt\_randinit.m\
Test that starts running Manopt on a limited number of cameras, without any initial guess given, and then using the results obtained as initial guesses for those cameras, progressively/iteratively increasing the number of cameras considered. \
OBS. Even if the filename is similar, the functioning/reasoning behind this test is in fact quite different.

- test\_manopt\_hessian\_rotation.m\
Test the rotation step analytically computed Hessian (and gradients) using Manopt "checkHessian" (and "checkGradient") functions

- test\_manopt\_hessian\_translation.m\
Test the translation step analytically computed Hessian (and gradient) using Manopt "checkHessian" (and "checkGradient") functions

- test\_hessian\_manopt\_genproc.m\
Test hessians of translation and rotation inside generalized Procrustes Manopt

- test\_rotation\_hessian\_geodesic.m\
Run a test similar to "test\_manopt\_step1\_gradients.m", this time checking also the Hessian

- test\_rotation\_hessian\_geodesic\_P.m\
Run a test similar to "test\_rotation\_hessian\_geodesic.m", only considering the $tr(R^T P)$ part of the cost

- test\_translation\_hessian\_geodesic.m\
Run a test similar to "test\_manopt\_step2\_gradients.m", this time checking also the Hessian

