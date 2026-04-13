#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include"../include/linear_algebra.h"

//#define DEBUG

/*
Test for the linear solvers
*/

int main(void)
  {
  int i, j;
  const int n = 10;
  double **A, *b, *x;
  double **Acopy, *bcopy, *auxvec;
  double test; 

  // initialize the random number generator with the time
  srand((unsigned int)time(NULL));

  // allocate A
  A=(double **)malloc((unsigned long int)(n)*sizeof(double*));
  if(A == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(i=0; i<n; i++)
     {
     A[i]=(double *)malloc((unsigned long int)(n)*sizeof(double));
     if(A[i] == NULL)
       {
       fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // allocate Acopy
  Acopy=(double **)malloc((unsigned long int)(n)*sizeof(double*));
  if(Acopy == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(i=0; i<n; i++)
     {
     Acopy[i]=(double *)malloc((unsigned long int)(n)*sizeof(double));
     if(Acopy[i] == NULL)
       {
       fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // allocate b, bcopy and x
  b=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(b == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  bcopy=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(bcopy == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  x=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(x == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  auxvec=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(auxvec == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  printf("\n");

  // ---------------------------------------
  printf("Gauss-Jordan solver for random matrix\n");

  // initialize everything
  for(i=0; i<n; i++)
     {
     for(j=0; j<n; j++)
        {
        A[i][j]=(double)((rand() % 20)-10); // integer random number in [-10, 10] converted in double
        }
     b[i]=(double)((rand() % 20)-10);
     }
  equal_mat(n, Acopy, A);  // Acopy=A
  equal_vec(n, bcopy, b);  // bcopy=b

  #ifdef DEBUG
  printf("Matrix A\n");
  for(i=0; i<n; i++)
     { 
     for(j=0; j<n; j++)
        {
        printf("%+3.0lf ", A[i][j]);
        }
     printf("\n");
     }
  printf("Vector b\n"); 
  for(i=0; i<n; i++) printf("%+3.0lf\n", b[i]);
  #endif

  // Gauss Jordan solution of Ax=b
  GaussJordan_fullpivot(n, A, b, x);

  // test the solution
  matvec_mult(n, auxvec, Acopy, x); // auxvec=Acopy*x
  minuseq_vec(n, auxvec, bcopy);    // auxvec-=bcopy
  test=sqrt(scalprod(n, auxvec, auxvec)); // test = ||auxvec||
  printf("Test of Gauss-Jordan solution\n");
  printf("|Ax-b|_2 = %g\n", test);
  printf("\n\n"); 

  // ---------------------------------------
  printf("Gauss-Seidel solver for tridiagonal matrix\n");

  // to test Gauss-Seidel we use a specific matrix
  for(i=0; i<n; i++)
     {
     A[i][i]=2.0;
     if(i+1<n)
       {
       A[i+1][i]=-1.0;
       A[i][i+1]=-1.0;
       }
     }

  // and for b a random vector
  for(i=0; i<n; i++)
     {
     b[i]=(double)((rand() % 20)-10);
     }

  #ifdef DEBUG
  printf("Matrice A\n");
  for(i=0; i<n; i++)
     { 
     for(j=0; j<n; j++)
        {
        printf("%+3.0lf ", A[i][j]);
        }
     printf("\n");
     }
  printf("Vettore b\n"); 
  for(i=0; i<n; i++) printf("%+3.0lf\n", b[i]);
  #endif

  // Gauss-Seidel solution of Ax=b
  GaussSeidel(n, A, b, x, 1.0e-7, 1000); 
 
  // test the solution
  matvec_mult(n, auxvec, A, x); // auxvec=A*x
  minuseq_vec(n, auxvec, b);    // aux-=b
  test=sqrt(scalprod(n, auxvec, auxvec)); // test = ||auxvec||
  printf("Test of Gauss-Seidel solution\n");
  printf("|Ax-b|_2 = %g\n", test);
  printf("\n\n"); 

  // ---------------------------------------
  printf("Conjugate gradient solver for tridiagonal matrix\n");

  // same matrix as for Gauss-Seidel

  #ifdef DEBUG
  printf("Matrice A\n");
  for(i=0; i<n; i++)
     { 
     for(j=0; j<n; j++)
        {
        printf("%+3.0lf ", A[i][j]);
        }
     printf("\n");
     }
  printf("Vettore b\n"); 
  for(i=0; i<n; i++) printf("%+3.0lf\n", b[i]);
  #endif

  conjugate_gradient(n, A, b, x, 1.0e-7, 200); 

  // test the solution
  matvec_mult(n, auxvec, A, x); // auxvec=A*x
  minuseq_vec(n, auxvec, b);    // aux-=b
  test=sqrt(scalprod(n, auxvec, auxvec)); // test = ||auxvec||
  printf("Test of conjugate gradient solution\n");
  printf("|Ax-b|_2 = %g\n", test);
  printf("\n\n"); 

  // ------------------------------- 
  // deallocate everything
  for(i=0; i<n; i++)
     {
     free(A[i]);
     free(Acopy[i]);
     }
  free(A);
  free(Acopy);
  free(b);
  free(bcopy);
  free(x);
  free(auxvec);

  return EXIT_SUCCESS;
  }

