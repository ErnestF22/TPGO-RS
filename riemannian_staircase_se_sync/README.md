- test\_newdescentdir\_geod.m 
Geodesic plot of cost function in the $f(x(t)) = x_0 + tv$ form, checking whether newly found descent direction actually decreases the cost (after an initial stale point);

- test\_yplus\_grad.m
At the end of every staircase step, Manopt finds us a point s.t. $grad_Y f(Y) = 0 \rightarrow $ check whether $grad_{Y_{plus}} f(Y_{plus}) == 0$, with $Y_{plus}$ computed as described above. This assures us that the critical points do not drastically change when going onto the next staircase step;

- test\_stiefel\_randgeodfun.m
Usage Example of $stiefel\_randgeodfun()$ Matlab function, checking whether all submats of the initial point of the generated geodesic are actually on the Stiefel manifold

- stiefel\_rand\_geod\_fun\_new\_test.m
Test/Usage example for new version of \textit{stiefel\_rand\_geod\_fun\_new.m} that substitutes old, partial versions of the Stiefel geodesic generating functions that were inside Codemeta and/or Manopt

- test\_stiefel\_der.m
Test derivative of a curve $f(t) = R_n(t) Y_0 R_p(t)$ on the Stiefel manifold, with $R_n \in SO(n), R_p \in SO(p)$, with $R_n, R_p$ as randomly generated geodesics in the respective manifolds