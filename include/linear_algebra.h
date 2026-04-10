#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

// Gauss-Jordan elimination with full pivoting
// to solve Ax=b with A a matrix of size n
//
// After use A->identity and b->permutation of x
void GaussJordan_fullpivot(int n, double **A, double *b, double *x);


#endif

