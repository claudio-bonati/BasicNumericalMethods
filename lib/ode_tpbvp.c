#include<stdio.h>
#include<stdlib.h>

#include"../include/linear_algebra.h"
#include"../include/ode_init.h"  // <--- DIM is defined here as a macro


// two point boundary value problem on [a, b] for a SCALAR equation (single component)
// -y'' + q(x)y(x) = f(x)
// with y(a)=0, y(b)=0
void tpbvp_fdm(double (*q)(double), 
               double (*f)(double),
               double a,
               double b,
               int nsteps,  // numbert of points of the discretization
               double *sol, // nsteps vector: sol[i]=solution in x=a+(i+1)h, with h=(b-a)/(nsteps+1)
                            // sol[0] is the solution in a+h
                            // sol[nsteps-1] is the solution in a+nsteps*(b-a)/(nsteps+1)=b-h
               double accuracy,    // accuracy and maxiter for Conjugate Gradient solver
               int maxiter)
  {
  int i;
  double **A;   // finite difference matrix
  double *fvec; // values of f stored in a vector
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

  // allocate fvec
  fvec=(double *)malloc((unsigned long int)(nsteps)*sizeof(double));
  if(fvec == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
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

     fvec[i]=f(a+h*(double)(i+1));
     }

  // solve A*sol=fvec
  conjugate_gradient(nsteps, A, fvec, sol, accuracy, maxiter); 

  for(i=0; i<nsteps; i++)
     {
     free(A[i]);
     }
  free(A);
  free(fvec);
  }


