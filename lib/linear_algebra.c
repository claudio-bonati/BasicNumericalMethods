#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>


#include"../include/linear_algebra.h"

//#define DEBUG

// inline functions defined in linear_algebra.h

// vector assignement: out=in
void equal_vec(int n, double *out, double *in);

// vector -=: a-=b
void minuseq_vec(int n, double *a, double *b);

// matrix assignement: out=in
void equal_mat(int n, double **out, double **in); 

// scalar product of two vectors
double scalprod(int n, double *v1, double *v2);

// matrix-vector multiplication: ris = A * x
void matvec_mult(int n, double *ris, double **A, double *x); 

// -----------------------

// Gauss-Jordan elimination with full pivoting
// to solve Ax=b with A a matrix of size n
//
// After use A->identity and b->permutation of x
void GaussJordan_fullpivot(int n,     // size of the matrix
                           double **A, 
                           double *b, 
                           double *x)
  {
  int i, j, k, pivotrow, pivotcol, tmpindex;
  double max, tmp, pivot;
  const double epsilon=1.0e-12;
  int *colindex; // to track column swaps

  colindex=(int *)malloc((unsigned long int)(n)*sizeof(int));
  if(colindex == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialized to no-swap
  for(i=0; i<n; i++) 
     {
     colindex[i] = i;
     }
  
   for(k=0; k<n; k++)  // loop on columns
      {
      // find the pivot (max abs value) in the submatrix with indices >=k
      max=0.0;
      pivotrow=k;
      pivotcol=k;
  
      for(i=k; i<n; i++) 
         {
         for(j=k; j<n; j++) 
            {
            if(fabs(A[i][j]) > max) 
              {
              max = fabs(A[i][j]);
              pivotrow = i;
              pivotcol = j;
              }
            }  
         }
  
      if(fabs(max) < epsilon)
        {
        fprintf(stderr, "Matrix is (nearly) singular (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
  
      // swap rows
      if(pivotrow != k)
        {
        for(j = 0; j < n; j++)
           {
           tmp = A[k][j];
           A[k][j] = A[pivotrow][j];
           A[pivotrow][j] = tmp;
           }
        tmp = b[k];
        b[k] = b[pivotrow];
        b[pivotrow] = tmp;
        }
  
      // swap columns
      if(pivotcol != k) 
        {
        for(i=0; i<n; i++) 
           {
           tmp = A[i][k];
           A[i][k] = A[i][pivotcol];
           A[i][pivotcol] = tmp;
          }
        tmpindex = colindex[k];
        colindex[k] = colindex[pivotcol];
        colindex[pivotcol] = tmpindex;
        }

      // now the pivot is in position [k][k]  

      // normalize pivot row
      pivot = A[k][k];
      for(j=0; j<n; j++) 
         {
         A[k][j] /= pivot;
         }
      b[k] /= pivot;
  
      // eliminate other rows
      for(i=0; i<n; i++) 
         {
         if(i != k)  // i = row index different from that of the pivot
           {
           tmp = A[i][k];
           for(j = 0; j < n; j++) 
              {
              A[i][j] -= tmp * A[k][j];
              }
           b[i] -= tmp * b[k];
           }
         }
      } // <- end of the loop on columns
  
  // reorder the solution due to column swaps
  for(i=0; i<n; i++) 
     {
     x[colindex[i]] = b[i];
     }

  free(colindex);
  }


// Gauss-Seidel iterative solution of Ax=b with A a matrix of size n
// (work for sure if A is symmetric and positive defined)
//
// A and b stay constant in this function
void GaussSeidel(int n,       // size of the matrix
                 double **A, 
                 double *b, 
                 double *x,
                 double accuracy, // elementwise accuracy
                 int maxiter) // maximum number of iterations;
  {
  int i, j, iter;
  double sum, err, diff;
  double *x_old;
  
  x_old=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(x_old == NULL)
    {
    fprintf(stderr, "GS: allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  for(i=0; i<n; i++)
     {
     x[i]=0.0; // step zero, anything would work
     }
  
  err=accuracy+1.0;
  for(iter=0; iter<maxiter && err>accuracy; iter++) 
     {
     // store previous iteration
     for(i=0; i<n; i++)
        {
        x_old[i]=x[i];
        }  

     // iteration
     for(i=0; i<n; i++) 
        {
        sum = 0.0;
        for(j=0; j<n; j++) 
           {
           if(j!=i)
             {
             sum += A[i][j] * x[j];
             }
           }
  
        x[i] = (b[i] - sum) / A[i][i];

        //// This is the modification needed for the Succesive Overrelaxation algorithm
        //double omega = 1.2; // typically between 1 and 2
        //x[i] = (1 - omega) * x[i] + omega * (b[i] - sum) / A[i][i];
        }
  
     // check convergence (componentwise)
     err = 0.0;
     for(i=0; i<n; i++) 
        {
        diff = fabs(x[i] - x_old[i]);
        if(diff>err)
          {
          err=diff;
          }
        }

     #ifdef DEBUG
     printf("GS: iter=%d err=%g\n", iter, err);
     #endif
     }
  
  free(x_old);

  if(iter==maxiter)
    {
    fprintf(stderr, "GS: reached maxiter=%d in (%s, %d)\n", maxiter, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// conjugate gradient solver of Ax=b with A a positive defined matrix of size n
//
// A and b stay constant in this function
void conjugate_gradient(int n,      // size of the matrix
                        double **A, 
                        double *b, 
                        double *x, 
                        double accuracy,  // accuracy of the solution: ||r_k||<accuracy
                        int maxiter)      // maximun number of iterations
  {
  int i, iter;
  double rs_new, rs_old, alpha;
  double *r, *p, *Ap;
  
  r=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(r == NULL)
    {
    fprintf(stderr, "CG: allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  p=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(p == NULL)
    {
    fprintf(stderr, "CG: allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  Ap=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(Ap == NULL)
    {
    fprintf(stderr, "CG: allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  matvec_mult(n, r, A, x);  // r=A*x
  for(i=0; i<n;i++) 
     {
     r[i]=b[i]-r[i];   
     p[i]=r[i];
     }
  // now p = r = b - A*x
 
  rs_old = scalprod(n, r, r);
  #ifdef DEBUG
  printf("CG: iter=-1, rs=%g\n", sqrt(rs_old));
  #endif
  
  for(iter=0; iter<maxiter && sqrt(rs_old)>accuracy; iter++) 
     {
     matvec_mult(n, Ap, A, p);  // Ap=A*p

     alpha = rs_old / scalprod(n, p, Ap);
  
     // x = x + alpha * p
     for(i=0; i<n; i++)
        {
        x[i] += alpha * p[i];
        }
  
     // r = r - alpha * Ap
     for(i=0; i<n; i++)
        {
        r[i] -= alpha * Ap[i];
        }
  
     rs_new = scalprod(n, r, r);
     #ifdef DEBUG
     printf("CG: iter=%d, rs=%g\n", iter, sqrt(rs_new));
     #endif
  
     // p = r + (rs_new / rs_old) * p
     for(i=0; i<n; i++)
        {
        p[i] = r[i] + (rs_new / rs_old) * p[i];
        }
  
     rs_old = rs_new;
     }

  free(r);
  free(p);
  free(Ap);

  if(iter==maxiter)
    {
    fprintf(stderr, "CG: eached maxiter=%d in (%s, %d)\n", maxiter, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// ------------------------

// determine the eigenvalue with maximum absolute value and its the corresponding eigenvector
// using the power method (A stays constant)
void maxeig_power(int n,      // size of the matrix
                  double **A, 
                  double *eigvalue,
                  double *eigvector,
                  double accuracy,  // accuracy of the eigenvalue
                  int maxiter)      // maximun number of iterations
  {
  int i, iter;  
  double norm, err, old, new;
  double *Ax;

  Ax=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(Ax == NULL)
    {
    fprintf(stderr, "MaxEigPow: allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize the random number generator with the time
  srand((unsigned int)time(NULL));
 
  // initialize eigvector with random numbers in [-10,10]
  for(i=0; i<n; i++)
     {
     eigvector[i]=(double)((rand() % 20)-10);
     }
  norm=sqrt(scalprod(n, eigvector, eigvector));
  for(i=0; i<n; i++)
     {
     eigvector[i]/=norm;
     }
  // now \|eigvector\|=1

  iter=0;
  do
     {
     matvec_mult(n, Ax, A, eigvector);   // Ax=A*eigvector

     new=scalprod(n, eigvector, Ax); // new=<eigvector, Ax>
     #ifdef DEBUG
     printf("MaxEigPow: iter=%d ; eigvalue=%.12lf\n", iter, new);
     #endif
  
     if(iter==0)
       {
       err=1.0;
       }
     else
       {
       err=fabs(new-old);
       }

     norm=sqrt(scalprod(n, Ax, Ax));
     for(i=0; i<n; i++)
        {
        eigvector[i]=Ax[i]/norm;  // now eigvector = Ax / \|Ax\|
        }

     old=new;
     iter++;
     }
  while(iter<maxiter && err>accuracy);

  *eigvalue=new; 

  free(Ax);
  
  if(iter==maxiter)
    {
    fprintf(stderr, "MaxEigPow: reached maxiter=%d in (%s, %d)\n", maxiter, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


// determine the eigenvalue with minimum absolute value and its the corresponding eigenvector
// using the inverse power method (A stays constant)
void mineig_power(int n,      // size of the matrix
                  double **A, 
                  double *eigvalue,
                  double *eigvector,
                  double accuracy,  // accuracy of the eigenvalue
                  int maxiter)      // maximun number of iterations
  {
  int i, iter;  
  double norm, err, old, new;
  double *Ax;

  Ax=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(Ax == NULL)
    {
    fprintf(stderr, "MinEigPow: allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize the random number generator with the time
  srand((unsigned int)time(NULL));
 
  // initialize eigvector with random numbers in [-10,10]
  for(i=0; i<n; i++)
     {
     eigvector[i]=(double)((rand() % 20)-10);
     }
  norm=sqrt(scalprod(n, eigvector, eigvector));
  for(i=0; i<n; i++)
     {
     eigvector[i]/=norm;
     }
  // now \|eigvector\|=1

  iter=0;
  do
     {
     conjugate_gradient(n, A, eigvector, Ax, 1.0e-10, 1000); // Ax = A^{-1}*eigvector
     // note that the starting point of the solution
     // is the solution at the previous iteration, 
     // which makes CG converge faster

     new=1.0/scalprod(n, eigvector, Ax); // new=<eigvector, Ax>
     #ifdef DEBUG
     printf("MinEigPow: iter=%d ; eigvalue=%.12lf\n", iter, new);
     #endif
  
     if(iter==0)
       {
       err=1.0;
       }
     else
       {
       err=fabs(new-old);
       }

     norm=sqrt(scalprod(n, Ax, Ax));
     for(i=0; i<n; i++)
        {
        eigvector[i]=Ax[i]/norm;  // now eigvector = Ax / \|Ax\|
        }

     old=new;
     iter++;
     }
  while(iter<maxiter && err>accuracy);

  *eigvalue=new; 

  free(Ax);
 
  if(iter==maxiter)
    {
    fprintf(stderr, "MinEigPow: reached maxiter=%d in (%s, %d)\n", maxiter, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  }


#undef DEBUG
