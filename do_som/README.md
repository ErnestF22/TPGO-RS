Many tests/comparisons required a $do\\_som()$ function that generally calls the Shape of Motion pipeline with 2 separate methods.
So, multiple $do\\_som...()$ with different appendices in the filename have been implemented:

- DO\_SOM \ Function that executes the Shape of Motion algorithms through Manopt and Procrustes pipelines, returning rotation and translation errors, as well as the mean execution time; the SoM algorithms are run on Gaussian noisy input data, with sigma variance and mu mean (that are being passed as arguments to the DO\_SOM function) \This version uses the two steps pipeline for Manopt, where optimization is done repetitively first on rotations and then on translations (in an ICP style pipeline).

- DO\_SOM\_GENPROC\ Function that executes the Shape of Motion algorithms through Manopt and Procrustes pipelines, returning rotation and translation errors, as well as the mean execution time; the SoM algorithms are run on Gaussian noisy input data, with sigma variance and mu mean (that are being passed as arguments to the DO\_SOM function) \This version uses the generalized Procrustes pipeline for Manopt, which optimizes at the same time for rotations and translations, opposite to the do\_som() original two-steps pipeline.

- DO\_SOM\_RGRAD\_TWOMODES\ Function that executes the Shape of Motion algorithms through Manopt pipelines, returning rotation and translation errors, as well as the mean execution time.\This $do\\_som()$ is for comparing the precision, accuracy and exectime using the automated Manopt Riemannian gradient estimation, against the computed-by-hand gradient.

- DO\_SOM\_MANOPT\_MAXITER\ Function that executes the Shape of Motion algorithms through Manopt pipelines, returning rotation and translation errors, as well as the mean execution time.\This $do\\_som()$ is for comparing the precision, accuracy and exectime when increasing the maximum potential number of ICP iterations.

- DO\_SOM\_CP2M\_NOISEINIT\ Function that executes the Shape of Motion algorithms through both Manopt pipelines as well as the Procrustes one, returning rotation and translation errors, as well as the mean execution time. CP2M stands for Compare Procrustes VS 2 Manopts.

