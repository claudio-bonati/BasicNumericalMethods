#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include"../include/linear_algebra.h"

// vector assignement: out=in
inline void equal_vec(int n, double *out, double *in)
  {
  int i;

  for(i=0; i<n; i++)
     {
     out[i]=in[i];
     }
  } 

// vector -=: a-=b
inline void minuseq_vec(int n, double *a, double *b)
  {
  int i;

  for(i=0; i<n; i++)
     {
     a[i]-=b[i];
     }
  } 

// matrix assignement: out=in
inline void equal_mat(int n, double **out, double **in) 
  {
  int i, j;

  for(i=0; i<n; i++) 
     {
     for(j=0; j<n; j++) 
        {
        out[i][j] = in[i][j];
        }
     }
  }

// scalar product of two vectors
inline double scalprod(int n, double *v1, double *v2) 
  {
  int i;
  double ris = 0.0;

  for(i=0; i<n; i++)
      ris += v1[i] * v2[i];
  return ris;
  }

// matrix-vector multiplication: ris = A * x
inline void matvec_mult(int n, double *ris, double **A, double *x) 
  {
  int i, j;

  for(i=0; i<n; i++) 
     {
     ris[i] = 0.0;
     for(j=0; j<n; j++) 
        {
        ris[i] += A[i][j] * x[j];
        }
     }
  }

// matrix-matrix multiplication: ris = A * B
inline void matmat_mult(int n, double **ris, double **A,  double **B) 
  {
  int i, j, k;
  double tmp;

  for(i=0; i<n; i++) 
     {
     for(j=0; j<n; j++) 
        {
        tmp=0.0;
        for(k=0; k<n; k++)
           { 
           tmp+=A[i][k]*B[k][j];
           }
        ris[i][j]=tmp;
        }
     }
  }

// matrix-matrix multiplication: ris = A^{t} * B
inline void matdagmat_mult(int n, double **ris, double **A, double **B) 
  {
  int i, j, k;
  double tmp;

  for(i=0; i<n; i++) 
     {
     for(j=0; j<n; j++) 
        {
        tmp=0.0;
        for(k=0; k<n; k++)
           { 
           tmp+=A[k][i]*B[k][j];
           }
        ris[i][j]=tmp;
        }
     }
  }


// Gauss-Jordan elimination with full pivoting
// to solve Ax=b with A a matrix of size n
//
// After use A->identity and b->permutation of x
void GaussJordan_fullpivot(int n,     // size of the matrix
                           double **A, 
                           double *b, 
                           double *x);


// Gauss-Seidel iterative solution of Ax=b with A a matrix of size n
//
// A and b stay constant in this function
void GaussSeidel(int n,       // size of the matrix
                 double **A, 
                 double *b, 
                 double *x,
                 double accuracy, // elementwise accuracy
                 int maxiter); // maximum number of iterations;

#endif

