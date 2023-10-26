\item \textit{eucl\_tangentProj\_test.m}\\Test projection on tangent space of the Euclidean space; this is $eucl\_proj(a,b) = b$;
\item \textit{test\_double\_stiefel\_proj\_composition.m}\\Check whether $Proj_x[u] \circ Proj_x[u] == Proj_x[u]$, as stated in Boumal's Intro book, page 111
\item \textit{test\_problem.m}\\Set up a \textit{problem} struct with cost, egrad, rgrad, ehess, rhess associated functions and L, P matrices fit for dimensions $sz_{problem} = [ num\_rows\_stiefel=4, d=3, N=5 ]$;
\item \textit{test\_problem\_curve.m}\\Create a pseudo-geodesic curve starting from $R^N Y_0 R^P$ with $R^N \in SO(N), R^P \in SO(P)$ and $Y_0 = stiefel_\_randn(sz_{problem})$;
\item \textit{test\_problem\_curve\_test.m}\\Check that the pseodo-geodesic curve is properly built, along with its first and second derivatives
\item \textit{test\_problem\_test\_egrad.m}\\Check that the expression for the Euclidean gradient is correct
\item \textit{test\_problem\_test\_ehess.m}\\Check that the expression for the Euclidean Hessian is correct
\item \textit{test\_problem\_test\_rgrad.m}\\Check that the expression for the Riemannian gradient is correct
\item \textit{test\_problem\_test\_rhess.m}\\Check that the expression for the Riemannian Hessian (computed with the formula reported in eq. \eqref{eq:riem_hess_proj_diff_proj_perp}) is correct
\item \textit{test\_problem\_test\_rhess\_canonical.m}\\Check that the expression for the Riemannian Hessian, built using the canonical metric instead of the Euclidean one, is correct
\item \textit{test\_problem\_test\_rhess\_rconn.m}\\Check that the expression for the Riemannian Hessian (computed through applying the Riemannian connection operator to $grad ~ f$) is correct
\item \textit{test\_rotGeodFun\_SoM.m}\\Test first and second derivatives of random geodesic built through $rot\_geodFun()$.
 