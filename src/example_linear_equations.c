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
  double **Acopy, *bcopy;
  double tmp, test; 

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
  bcopy=(double *)malloc((unsigned long int)(n)*sizeof(double));
  x=(double *)malloc((unsigned long int)(n)*sizeof(double));

  // initialize everything
  for(i=0; i<n; i++)
     {
     for(j=0; j<n; j++)
        {
        A[i][j]=(double)((rand() % 20)-10); // integer random number in [-10, 10] converted in double
        Acopy[i][j]=A[i][j];
        }
     b[i]=(double)((rand() % 20)-10);
     bcopy[i]=b[i];
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


  // Gauss Jordan solution of Ax=b
  GaussJordan_fullpivot(n, A, b, x);

  // test the solution
  test=0.0;
  for(i=0; i<n; i++)
     {
     tmp=0.0;
     for(int j=0; j<n; j++)
        {
        tmp+=Acopy[i][j]*x[j];
        }
     tmp-=bcopy[i];
     test+=fabs(tmp);
     }
  printf("Test of Gauss-Jordan solution\n");
  printf("|Ax-b|_1 = %lg\n", test);

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

  return EXIT_SUCCESS;
  }

