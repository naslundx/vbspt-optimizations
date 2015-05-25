# vbspt-optimizations
Optimization efforts for the vbSPT software for the Master's thesis "High Performance Computing aspects of Single Particle Machine Learning" written during the Spring of 2015 by Marcus Näslund.

There are seven files present:
* **VB3_VBEMiterator.m**: The most optimized version of the EM algorithm.
* **VB3_preprocess.m**: The new version of the preprocessing algorithm. Can be parallelized if necessary.
* **HMM_multiForwardBackward.c**: An optimized version of the FB algorithm with no functional changes. 
* **HMM_multiForwardBackward_openmp.c**: A version of the FB algorithm parallelized using OpenMP. Remember to turn off the parallelization in MATLAB before using this.
* **HMM_multiForwardBackward_blas.c**: A version of the FB algorithm using BLAS routines. Remember to use the VB3_VBEMiterator_blas.m and VB3_preprocess_blas.m instead of the original files, so that data is ordered in the correct way.
* **VB3_VBEMiterator_blas.m**: See above.
* **VB3_preprocess_blas.m**: See above.

For details on performance please read the thesis report (link coming in June).

Comments can be sent to marcus.naslund.7959 (at) student.uu.se.
