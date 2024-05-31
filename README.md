# REFERENCES
[1] C. C. Paige and M. A. Saunders, LSQR: An algorithm for sparse linear equations and sparse
least squares, ACM Transactions on Mathematical Software (TOMS), 8 (1982), pp. 43–71.

[2] Valeria Simoncini and Lorenzo Piccinini. TRUNCATED LSQR FOR
MATRIX LEAST SQUARES PROBLEMS AND APPLICATION TO
DICTIONARY LEARNING *. working paper or preprint, February
2024. https://hal.science/hal-04437719/


# LSQR

LSQR was first introduced in by Piage and Saunders to solve linear linear systems $Ax=b$ and least squares problems $\min ||Ax-b||_2^2$. 

In our work we present a matrix-oriented version of the algorithm that works wit low rank approximations. We implement the truncation function such that it takes advantage of previous factorizations, leading to 
a cheaper step that otherwise would have been memory and cost requiring.


In order to run the test, rememeber to download first the dataset MNIST. 
