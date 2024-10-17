- eucl\_tangentProj\_test.m\
Test projection on tangent space of the Euclidean space; this is $eucl\\_proj(a,b) = b$;

- make\_rand\_stiefel\_3d\_array.m\
Make random Stiefel 3D array

- plot\_pim.m\
Plot functions starting from PIM output, towards found "eigenvector" direction

- riemannian\_connection\_stiefel.m\
Apply the Riemannian Connection operator for an array on the 3D-Stiefel manifold

- test\_double\_stiefel\_proj\_composition.m\
Check whether $Proj_x[u] \circ Proj_x[u] == Proj_x[u]$, as stated in Boumal's Intro book, page 111

- test\_problem.m\
Set up a "problem" struct with cost, egrad, rgrad, ehess, rhess associated functions and L, P matrices fit for dimensions $sz_{problem} = [ num\\_rows\\_stiefel=4, d=3, N=5 ]$;

- test\_problem\_curve.m\
Create a pseudo-geodesic curve starting from $R^N Y_0 R^P$ with $R^N \in SO(N), R^P \in SO(P)$ and $Y_0 = stiefel\_{randn} (sz_{problem})$;

- test\_problem\_curve\_test.m\
Check that the pseodo-geodesic curve is properly built, along with its first and second derivatives

- test\_problem\_test\_egrad.m\
Check that the expression for the Euclidean gradient is correct

- test\_problem\_test\_ehess.m\
Check that the expression for the Euclidean Hessian is correct

- test\_problem\_test\_hess\_antidiag.m\
Check whether the antidiagonal part of the Hessian, computed as $\frac{\partial^2 R}{\partial T^2}$ is the same as the one computed by $\frac{\partial^2 T}{\partial R^2}$

- test\_problem\_test\_pim\_nextstep.m\
Test curvature of what PIM Hessian is able to extract

- test\_problem\_test\_rgrad.m\
Check that the expression for the Riemannian gradient is correct

- test\_problem\_test\_rhess.m\
Check that the expression for the Riemannian Hessian (computed with the formula reported in eq. \eqref{eq:riem_hess_proj_diff_proj_perp}) is correct

- test\_problem\_test\_rhess\_canonical.m\
Check that the expression for the Riemannian Hessian, built using the canonical metric instead of the Euclidean one, is correct

- test\_problem\_test\_rhess\_manopt.m\
Check that the expression for the Riemannian Hessian (computed through applying the Riemannian connection operator to $grad ~ f$) is correct through Manopt's $checkHessian()$

- test\_rand\_stiefel\_generator.m\
Test $make\_rand\_stiefel\_array()$ function

- test\_rotGeodFun\_SoM.m\
Test first and second derivatives of random geodesic built through $rot\\_geodFun()$.
 
