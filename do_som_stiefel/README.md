The $do\_som()$ functions have their equivalent on the Stiefel manifold, on files whose name is obtained by adding the $\_stiefel$ suffix.
Generally, the Manopt manifold constructed has been of size $num\_rows\_stiefel \times d \times N$. $num\_rows\_stiefel$ is passed as an argument to the SoM-Stiefel functions.
\\
Note that, for the moment, only step one of Manopt has been considered/moved on Stiefel.   
%
\subsubsection{riemannian\_staircase\_se\_sync folder README file}
\begin{itemize}
    \item \textit{check\_yplus\_grad.m}\\Returns boolean that states whether $max(Yplus, [~], 'all')$ is close enough to zero (i.e., $< eps$)
    \item \textit{cpm.m}\\Compare Procrustes, Manopt and Manopt with the SE-Sync Riemannian Staircase addition over the course of multiple noisy tests, and averaging the rotation estimation, translation estimation and execution time results for each sigma
    \item \textit{create\_stiefel\_geod\_example.m}\\Usage example for creating a geodesic on Stiefel manifold
    \item \textit{lastNRowsAllZeros.m}\\Util function that checks whether all last N rows in a written point on the 3D-$Stiefel^N$ manifold are equal to 0; this means that the first $d$ rows of each submatrix can be directly taken as the $SO(d)$ projection of that submatrix ($\pm$ determinant sign to be precise)
    \item \textit{make\_rand\_stiefel\_mat.m}\\Function that creates a random matrix on a 3D Stiefel manifold (i.e., every matrix of the third dimension is on the Stiefel manifold specified in the function's input parameters)
    \item \textit{pam.m}\\Function that returns an eigenvector associated to the maximum eigenvalue of a squared matrix passed as input; tested in \textit{test\_pam.m}
    \item \textit{pam\_hessian.m}\\Function that returns an the minimum eigenvalue of a Hessian calculated on point $x\in St(d,p)^m$ manifold; tested in \textit{test\_pam\_hessian.m} 
    \item \textit{pam\_hessian\_sketch.m}\\Sketch code that useful as guidance for implementation of \textit{pam\_hessian()} pipeline
    \item \textit{riemannian\_staircase\_se\_sync.m}\\Function called from \textit{riemannian\_staircase\_se\_sync\_main.m} that does the work for the SoM variation of the SE-Sync Riemannian Staircase
    \item \textit{riemannian\_staircase\_se\_sync\_main.m}\\Run a single instance of the SoM variation of the SE-Sync Riemannian Staircase
    \item \textit{riemannian\_staircase\_se\_sync\_rotonly.m}\\Function called from \textit{riemannian\_staircase\_se\_sync\_rotonly\_main.m} that does the work for the SoM variation of the SE-Sync Riemannian Staircase, considering only the rotation estimation part of SoM (translation estimation is left to other tests)
    \item \textit{riemannian\_staircase\_se\_sync\_rotonly\_main.m}\\Run a single instance of the SoM variation of the SE-Sync Riemannian Staircase, considering only the rotation estimation part of SoM
    \item \textit{startup.m}\\Used to import Manopt from inside SE-Sync, as it appears that newer versions of Manopt cause some conflicts when running the SE-Sync Riemannian Staircase
    \item \textit{test\_newdescentdir\_geod.m}\\Tests whether $Y_{opt}$ is a critical (minimum) point, and the new direction calculated from $Y_{plus}$ leads instead to decreasing the cost, generating appropriate geodesics starting from those points and progressing onto the old/new directions
    \item \textit{test\_pam.m}\\Test power augmentation method for computing minimum eigenvalue of a square matrix
    \item \textit{test\_pam\_hessian.m}\\Sketch test for running power augmentation method and computing minimum eigenvalue of a linear map for computing the Hessian on a given point $x$
    \item \textit{test\_stiefel\_randgeodfun.m}\\Example usage for the $stiefel\_randgeodfun()$ function, that also tests it correct functioning
    \item \textit{test\_yplus\_grad.m}\\Test gradient at $Y_{plus}$ i.e., the initial point of descent for the new staircase step, in order to check whether that is still 0, after only augmenting it with zeros, before calculating new descent direction
\end{itemize}