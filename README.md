# LSQR

LSQR was first introduced in by Piage and Saunders to solve linear linear systems $Ax=b$ and least squares problems $\min ||Ax-b||_2^2$. 

In our work we present a matrix-oriented version of the algorithm that works wit low rank approximations. We implement the truncation function such that it takes advantage of previous factorizations, leading to 
a cheaper step that otherwise would have been memory and cost requiring.

If you use any of the code in the repository, please cite [1].

In order to run the test, download first the dataset MNIST. 

# REFERENCES
[1] L. Piccinini, and V. Simoncini. Truncated LSQR for matrix least squares problems, 
Computational Optimization and Applications 91 (2), 905-932.
DOI: https://doi.org/10.1007/s10589-024-00629-w

[2] C. C. Paige and M. A. Saunders, LSQR: An algorithm for sparse linear equations and sparse
least squares, ACM Transactions on Mathematical Software (TOMS), 8 (1982), pp. 43–71.
