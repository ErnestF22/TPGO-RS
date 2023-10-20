From "PROGRESSIVE TESTS (CODEMETA’S TESTNETWORK INPUTS)" section on project's overleaf report

\item \textit{test\_procrustes\_dummy.m} test that contains a trivial set/setup of points that can be translated and rotated ad libitum in order to test the correctness of any implemented Procrustes method (i.e. Kabsch, Umeyama, SoM, ...)
\item \textit{test\_rotation\_cost.m} \\ the cost function of step 1, when rotations are given and equal to ground truth, must be equal to 0
\item \textit{test\_translation\_cost.m} \\ the cost function of step 2, when translation are given and equal to ground truth, must be equal to 0 \\
Without noise, and with known rotations (e.g. setting all rotations equal to gt, with all rotations equal to $I$ in the first try), $A T - B = \mathbf{0}$ must be true, with $A$, $B$ as in \eqref{eq:Lij_step1} and \eqref{eq:Pij_step1} and $T$ filled with ground truth translation values.
\item \textit{test\_make\_LT\_PT.m} \\ Form L(T) and P(T) with two separate methods (matricial and loop) and check if they come out as equal.   
\item \textit{test\_make\_A\_b.m} \\ Equivalent of \textit{test\_make\_LT\_PT.m}
\item \textit{test\_make\_M\_N.m} \\ Equivalent of \textit{test\_make\_LT\_PT.m}
\item \textit{funCheckDer\_test.m} \\ test taken from codemeta, to be adapted and reused for testing gradients
\item Manopt gradients \\ Test correctness of translation estimation gradients hand-written expression, using an adapted version of codemeta's funCheckDer \\
In particular, when it comes to gradients, we test:
\begin{itemize}
    \item Euclidean gradient wrt R \\ (in \textit{test\_manopt\_step1\_gradients.m})
    \item Riemannian gradient wrt R \\ (in \textit{test\_manopt\_step1\_gradients.m})
    \item Euclidean gradient wrt T \\ (in \textit{test\_manopt\_step2\_gradients.m})
\end{itemize}
Without noise, they must be equal to 0\\
In code, $g$ indicates Euclidean gradient, $h$ indicates Riemannian gradient.
\item \textit{test\_procrustes\_fullgraph\_groundtruth.m} \\ Test SoM Procrustes considering that all data given is correct and that the graph is full (i.e. all edges are available)
\item \textit{test\_cpm\_single\_nonoise.m} \\ Compare Procrustes and Manopt on a single run of a test without noise\\\textit{Note: cpm is an acronym for Compare Procrustes and Manopt}
\item \textit{test\_cpm\_single.m} \\ Compare Procrustes and Manopt on a single run of a noise test (fixed low $\sigma$ i.e. $\sigma = 0.1$, null $\mu$ "added" inside Codemeta's testNetwork)
\item \textit{test\_manopt\_sep\_gen.m} \\ Compare transf out results obtained with the two Manopt pipelines: separated, generalized
\item \textit{test\_procrustes\_randinit.m} \\ Test created to verify Procrustes pipeline performance when dealing with increasingly noisy inputs (which should lead the estimation to fail)
\item \textit{compare\_manopt\_maxiter.m} \\
Test the improvement obtained on Manopt results by increasing the maximum number of ICP iterations
\item \textit{compare\_manopt\_riemgrad.m} \\
Compare the results obtained through Manopt when using auto-computed Riemannian gradient vs the hand-computed one; the comparison includes also execution times; graphical results from this test can be found in figure \ref{fig:manopt_analytical_vs_ad_grad}
\item \textit{compare\_manopt\_hessian\_ad\_analytic.m}\\Test Manopt results and execution times obtained when we use Manopt with automatic differentiation for computing Hessian, VS analytically pre-computed Hessian; graphical results from this test can be found in figure \ref{fig:manopt_analytical_vs_ad_hessian}
\item \textit{test\_make\_input\_subset.m} \\
Test that checks if the $make\_input\_subset()$ function works. Said function returns a partial input of a full testdata input generated through the testNetwork methods; it returns all data associated to the $cameras\_ids$ passed as input, extracting it by
removing the data associated to other cameras.
\item \textit{test\_manopt\_randinit.m}\\Test that starts running Manopt on a limited number of cameras, without any initial guess given, and then using the results obtained as initial guesses for those cameras, progressively/iteratively increasing the number of cameras considered. 
\\OBS. Even if the filename is similar, the functioning/reasoning behind this test is in fact quite different.
\item \textit{test\_manopt\_hessian\_rotation.m}\\Test the rotation step analytically computed Hessian (and gradients) using Manopt \textit{checkHessian} (and \textit{checkGradient}) functions
\item \textit{test\_manopt\_hessian\_translation.m}\\Test the translation step analytically computed Hessian (and gradient) using Manopt \textit{checkHessian} (and \textit{checkGradient}) functions
\item \textit{test\_hessian\_manopt\_genproc.m}\\Test hessians of translation and rotation inside generalized Procrustes Manopt
\item \textit{test\_rotation\_hessian\_geodesic.m}\\Run a test similar to \textit{test\_manopt\_step1\_gradients.m}, this time checking also the Hessian
\item \textit{test\_rotation\_hessian\_geodesic\_P.m}\\Run a test similar to \textit{test\_rotation\_hessian\_geodesic.m}, only considering the $\trace (R\transpose P)$ part of the cost
\item \textit{test\_translation\_hessian\_geodesic.m}\\Run a test similar to \textit{test\_manopt\_step2\_gradients.m}, this time checking also the Hessian