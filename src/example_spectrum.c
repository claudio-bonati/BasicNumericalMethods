#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include"../include/linear_algebra.h"

//#define DEBUG

/*
Test for eigensolvers
*/

int main(void)
  {
  int i, j;
  const int n = 30;
  double **A, *x;
  double eigvalue; 
  const double pi=3.141592653589793238462643383279502884197169399375105820;

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

  // allocate x 
  x=(double *)malloc((unsigned long int)(n)*sizeof(double));
  if(x == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize A with a tridiagonal matrix
  for(i=0; i<n; i++)
     {
     for(j=0; j<n; j++)
        {
        A[i][j]=0.0;
        }
     }
  for(i=0; i<n; i++)
     {
     A[i][i]=2.0;
     if(i+1<n)
       {
       A[i+1][i]=-1.0;
       A[i][i+1]=-1.0;
       }
     }

  printf("\n");

  // ---------------------------------------
  printf("Maximum egenvalue using power iteration for tridiagonal matrix\n");

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
  #endif

  // maximum eigenvalue using power iterations
  maxeig_power(n, A, &eigvalue, x, 1.0e-7, 1000);

  // test the solution
  printf("eig      = %.12lf\n", eigvalue);
  printf("analytic = %.12lf\n", 4.0*pow(sin(pi*(double)n/(2.0*(double)(n+1))),2.0));
  printf("\n\n"); 

  // all the eigenvalues of the discrete second derivative 
  // are 4*sin^2(pi j/[2*(n+1)])
  // see, e.g., D'Elia "Lezioni di meccanica classica" sez. 4.4.1
  // or
  // https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors_of_the_second_derivative

  // ---------------------------------------
  printf("Minimum egenvalue using inverse power iteration for tridiagonal matrix\n");

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
  #endif

  // maximum eigenvalue using power iterations
  mineig_power(n, A, &eigvalue, x, 1.0e-7, 1000);

  printf("eig      = %.12lf\n", eigvalue);
  printf("analytic = %.12lf\n", 4.0*pow(sin(pi/(2.0*(double)(n+1))),2.0));
  printf("\n\n"); 


  // ------------------------------- 
  // deallocate everything
  for(i=0; i<n; i++)
     {
     free(A[i]);
     }
  free(A);
  free(x);

  return EXIT_SUCCESS;
  }

