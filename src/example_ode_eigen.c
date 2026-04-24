#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define STRING_LENGTH 50

#include"../include/ode_init.h"  // <---- DIM is #defined here!
#include"../include/ode_eigen.h" 
#include"../include/root_find.h"

// two point boundary value problem on [a, b] for a SCALAR equation (single component)
// -y'' + q(x)y(x) = f(x)
// with y(a)=0, y(b)=0

double q(double x)
  {
  return x*x;
  }


// ------------------- FOR SHOOTING -------------------------

//  we write 
// -y'' + q(x)y(x) = lambda y(x)
// as
// y[0]'=y[1]
// y[1]'= q(x) y[0] - lambda y[0]
void for_rhs(double x, double y[DIM], double dydx[DIM], double lambda)
  {
  #if DIM==2
    dydx[0]=y[1];
    dydx[1]=q(x)*y[0]-lambda*y[0]; 
  #else
    fprintf(stderr, "Function valid only for DIM=2 and DIM=%d (%s, %d)\n", DIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  #endif 
  }


// perform the shoot
double func_for_shooting(double lambda, 
                         double a,
                         double b,
                         int nsteps,
                         double **sol)  // sol[nsteps+1][2] solution of the problem if sol[nsteps][0]==0
  {
  double y0[2];

  y0[0]=0.0;
  y0[1]=(b-a)/(double)(nsteps+1); // derivative 1  

  RK4_param(for_rhs, lambda, a, y0, b, nsteps, sol); 

  return sol[nsteps][0];
  }


typedef struct {
  double a;
  double b;
  int nsteps; 
  double **sol;
} Params;


// function to be used when only y0prime is important
double func_for_secant_with_context(double lambda, void *ctx)
  {
  Params *p=ctx;
  
  return func_for_shooting(lambda, p->a, p->b, p->nsteps, p->sol);
  }


// ---------------------------


int main(int argc, char **argv)
  {
  int i, j, nsteps_fde, nsteps_shoot; 
  const double a=-10;
  const double b=10;
  double h, norm, lambda;
  char datafile[STRING_LENGTH];
  FILE *fp;
  Params ctx;

  double *eigvals_fde;   //  eigvals_fde[steps] eigenvalues 
  double **eigvects_fde; //  eigvects_fde[steps][steps] eigenvectors
                         //  eigvects_fde[i][j] is the eigfunction corresponding to eigvals_fde[j]
                         //  evaluated at x=a+(i+1)h, with h=(b-a)/(nsteps_fde+1)

  double **sol_shoot; // sol_shoot[nsteps_fde+1][2]

  // check input from command line
  if(argc != 3)
    {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s nsteps_fde nsteps_shoot\n\n", argv[0]);
    fprintf(stdout, "  nsteps_fde = number of steps for FDE\n");
    fprintf(stdout, "  nsteps_shoot = number of integration steps in shooting\n");
    return EXIT_SUCCESS;
    }
  else
    {  
    // read input value
    nsteps_fde=atoi(argv[1]);
    nsteps_shoot=atoi(argv[2]);
    }

  // allocate eigvals_fde and eigvects_fde
  eigvals_fde=(double *)malloc((unsigned long int)(nsteps_fde)*sizeof(double));
  if(eigvals_fde == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  eigvects_fde=(double **)malloc((unsigned long int)(nsteps_fde)*sizeof(double*));
  if(eigvects_fde == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  for(i=0; i<nsteps_fde; i++)
     {
     eigvects_fde[i]=(double *)malloc((unsigned long int)(nsteps_fde)*sizeof(double));
     if(eigvects_fde[i] == NULL)
       {
       fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
       return EXIT_FAILURE;
       }
     }

  // allocate sol_shoot 
  sol_shoot=(double **)malloc((unsigned long int)(nsteps_shoot+1)*sizeof(double*));
  if(sol_shoot == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  for(i=0; i<=nsteps_shoot; i++)
     {
     sol_shoot[i]=(double *)malloc((unsigned long int)(DIM)*sizeof(double));
     if(sol_shoot[i] == NULL)
       {
       fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
       return EXIT_FAILURE;
       }
     }

  // ----------- FINITE DIFFERENCE METHOD ------------ 

  // initialize everything to zero
  for(i=0; i<nsteps_fde; i++)
     {
     eigvals_fde[i]=0.0;
     for(j=0; j<nsteps_fde; j++)
        {
        eigvects_fde[i][j]=0.0;
        }
     }

  // finite difference solution of the two point boundary value problem
  eigen_fdm(q, a, b, nsteps_fde, eigvals_fde, eigvects_fde, 1.0e-6, 100000);

  h=(b-a)/(double)(nsteps_fde+1);

  // normalize eigenvectors
  for(i=0; i<nsteps_fde; i++)
     {
     norm=0.0;
     for(j=0; j<nsteps_fde; j++)
        {
        norm+=pow(eigvects_fde[j][i],2.0)*h;
        }
     norm=sqrt(norm);
     for(j=0; j<nsteps_fde; j++)
        {
        eigvects_fde[j][i]/=norm;
        }
     }

  //print results
  sprintf(datafile, "ris_fdm_%d.dat", nsteps_fde);

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  fprintf(fp,"#%d %.12lf %.12lf\n", nsteps_fde, eigvals_fde[0], eigvals_fde[2]);

  for(i=0; i<nsteps_fde; i++)
     {
     fprintf(fp,"%.12lf %.12lf %.12lf\n", a+h*(double)(i+1), eigvects_fde[i][0], eigvects_fde[i][1]);
     }

  fclose(fp);

  // ----------- SHOOTING METHOD ------------ 

  lambda=1.0; // eigvals_fde[0];

  ctx.a=a;
  ctx.b=b;
  ctx.nsteps=nsteps_shoot;
  ctx.sol=sol_shoot;

  lambda=secant_ctx(func_for_secant_with_context,
             (void*)&ctx,
             lambda,
             lambda*1.2,
             1.0e-8,
             1000); 

  // normalize the solution
  h=(b-a)/(double)(nsteps_shoot+1);

  norm=0.0;
  for(i=1; i<nsteps_shoot; i++)
     {
     norm+=pow(sol_shoot[i][0],2.0)*h;
     }
  norm=sqrt(norm);
  for(i=1; i<nsteps_shoot; i++)
     {
     sol_shoot[i][0]/=norm;
     sol_shoot[i][1]/=norm;
     }

  //print results
  sprintf(datafile, "ris_shoot_%d.dat", nsteps_shoot);

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  fprintf(fp,"#%d %.12lf \n", nsteps_shoot, lambda);

  for(i=0; i<=nsteps_shoot; i++)
     {
     fprintf(fp,"%.12lf %.12lf\n", a+h*(double)i, sol_shoot[i][0]);
     }

  fclose(fp);

  // ------------------------------
 
  free(eigvals_fde);
  for(i=0; i<nsteps_fde; i++)
     {
     free(eigvects_fde[i]);
     }
  free(eigvects_fde);

  for(i=0; i<=nsteps_shoot; i++)
     {
     free(sol_shoot[i]);
     }
  free(sol_shoot);

  return EXIT_SUCCESS;
  }


