#include<stdio.h>
#include<stdlib.h>

#include"../include/linear_algebra.h"
#include"../include/ode_init.h"  // <--- DIM is defined here as a macro

// eigenvalue problem on [a, b] for a SCALAR equation (single component)
// -y'' + q(x)y(x) = lambda y(x)
// with y(a)=0, y(b)=0
//
// eigvals are sorted in such a way that eigvals[0]<=eigvals[1]<=eigvals[2]<=...
// and eigvects[:][i] is the eigenvector corresponfing to eigvals[i]
void eigen_fdm(double (*q)(double), 
               double a,
               double b,
               int nsteps,  // numbert of points of the discretization
               double *eigvals,
               double **eigvects, 
               double accuracy,    // accuracy and maxiter for Jacobi diagonalization
               int maxiter)
  {
  int i;
  double **A;   // finite difference matrix
  const double h=(b-a)/(double)(nsteps+1);

  // allocate A
  A=(double **)malloc((unsigned long int)(nsteps)*sizeof(double*));
  if(A == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(i=0; i<nsteps; i++)
     {
     A[i]=(double *)malloc((unsigned long int)(nsteps)*sizeof(double));
     if(A[i] == NULL)
       {
       fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  // initialize variables
  for(i=0; i<nsteps; i++)
     {
     A[i][i]=2.0/h/h + q(a+h*(double)(i+1));
     if(i+1<nsteps)
       { 
       A[i+1][i]=-1.0/h/h;
       A[i][i+1]=-1.0/h/h;
       }
     }

  // diagonalize A
  Jacobi_diag(nsteps, A, eigvals, eigvects, accuracy, maxiter);

  for(i=0; i<nsteps; i++)
     {
     free(A[i]);
     }
  free(A);
  }


