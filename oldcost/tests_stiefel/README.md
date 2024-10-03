From "PROGRESSIVE TESTS (CODEMETA’S TESTNETWORK INPUTS)" section on project's overleaf report:

- test\_manopt\_stiefel\_vs\_som.m\
Test whether Manopt gives same results optimizing on $d \times d \times N$ Stiefel manifold, against $SO(d)^N$

- Manopt gradients\
Test correctness of translation estimation gradients hand-written expression, using an adapted version of codemeta's funCheckDer \
In particular, when it comes to gradients, we test:
    - Euclidean gradient wrt R\
    (in "test\_step1\_gradients\_geodesic\_stiefel.m)
    - Riemannian gradient wrt R\
    (in "test\_step1\_gradients\_geodesic\_stiefel.m)
    - Euclidean gradient wrt T\
    (in "test\_step2\_gradients\_geodesic\_stiefel.m)
    - Hessian wrt R\
    (in "test\_step1\_gradients\_geodesic\_stiefel.m)
    - Hessian wrt T\
    (in "test\_step2\_gradients\_geodesic\_stiefel.m)\
Without noise, they must be equal to 0.\
In code, $g$ indicates Euclidean gradient, $h$ indicates Riemannian gradient.

- test\_cpm\_single\_nonoise\_stiefel.m\
Compare Procrustes and Manopt on a single run of a test without noise\ Note: cpm is an acronym for Compare Procrustes and Manopt

- test\_hessian\_step1\_geodesic.m\
Test expression for Hessian on a pseudo-geodesic curve
$ \ddot{Y}(t)=\ddot{R_L}(t) Y_0 R_R(t) +  2 \dot{R_L}(t) Y_0 \dot{R_R}(t) + R_L(t) Y_0 \ddot{R_R}(t) $
this test also includes a check on if the symmetry property is respected

- test\_hessian\_step1\_symmetry.m\
Checks whether Hessian is symmetric i.e., for $x \in St(d,p)^m$ ~ and ~ $v1, v2 \in \bot_x St(d,p)^m$ ~ whether
$ <v1, Hess(x)[v2]> == <v2, Hess(x)[v1]> $

- test\_cpm\_single\_stiefel.m\
Compare Procrustes and Manopt on a single run of a noise test (fixed low $\sigma$ i.e. $\sigma = 0.1$, null $\mu$ "added" inside Codemeta's testNetwork)

- test\_manopt\_sep\_gen\_stiefel.m\
Compare transf out results obtained with the two Manopt pipelines: separated, generalized

- compare\_manopt\_maxiter\_stiefel.m\
Test the improvement obtained on Manopt results by increasing the maximum number of ICP iterations

- compare\_manopt\_riemgrad\_stiefel.m\
Compare the results obtained through Manopt when using auto-computed Riemannian gradient vs the hand-computed one; the comparison includes also execution times

- test\_manopt\_randinit\_stiefel.m\
Test that starts running Manopt on a limited number of cameras, without any initial guess given, and then using the results obtained as initial guesses for those cameras, progressively/iteratively increasing the number of cameras considered.\
OBS. Even if the filename is similar, the functioning/reasoning behind this test is in fact quite different.

- test\_manopt\_hessian\_ad\_analytic\_stiefel.m\
Test Manopt results and execution times obtained when we use Manopt with automatic differentiation for computing Hessian, VS analytically pre-computed Hessian --- TODO: create test file and everything

- test\_make\_LT\_PT\_stiefel.m\
Form "Stiefel-expanded version" of L(T) and P(T) with two separate methods (matricial and loop) and check if they come out as equal. Brief comparison also with the non-Stiefel versions of L, P.

- test\_make\_A\_b\_stiefel.m\
Equivalent of "test\_make\_LT\_PT\_stiefel.m 

- test\_riemannian\_staircase.m\
Run a single test of the Riemannian staircase pipeline on generated dataset (through Codemeta's testNetwork)

- test\_stiefel\_cost\_step1.m\
Test whether step 1 (ref. eq. \eqref{eq:coord_desc_first_step}) costs defined with R on $SO(3)^N$ and, respectively, Stiefel manifold, are the same, for a randomly generated R.

- test\_stiefel\_cost\_step2.m\
Test whether step 2 (ref. eq. \eqref{eq:coord_desc_second_step}) costs defined with R on $SO(3)^N$ and, respectively, Stiefel manifold, are the same, for a randomly generated R.

- test\_manopt\_step1\_gradients\_stiefel.m\
Test Manopt step 1 gradients and Hessians through Manopt built-in "check...()" functions with $num\\_rows\\_stiefel = d+1$

- test\_manopt\_step2\_gradients\_stiefel.m\
Test Manopt step 2 gradients and Hessians through Manopt built-in "check...()" functions $num\\_rows\\_stiefel = d+1$

- test\_stiefel\_randn.m\
Just to explore how to make "stiefel\_randn()" and "stiefel\_randGeodFun()" work


All of these serve to form the following functions:

- som\_manopt\_genproc\_stiefel.m

- som\_manopt\_stiefel.m

- som\_stepone\_manopt\_stiefel.m

- som\_steptwo\_manopt\_stiefel.m 

