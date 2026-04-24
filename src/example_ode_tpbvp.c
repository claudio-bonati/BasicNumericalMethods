#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define STRING_LENGTH 50

#include"../include/ode_init.h"  // <---- DIM is #defined here!
#include"../include/ode_tpbvp.h" 
#include"../include/root_find.h"

// two point boundary value problem on [a, b] for a SCALAR equation (single component)
// -y'' + q(x)y(x) = f(x)
// with y(a)=0, y(b)=0

// sin(pi*x) is the exact solution of 
// -y''+x y(x)=(pi^2+x)sin(pi*x) 
// with y(0)=y(1)=0

// y(x)=x(1-x) is the exact solution of
// -y''+x y(x)=2+x^2(1-x)

double q(double x)
  {
  return x;
  }

double f(double x)
  {
  //const double pi=3.141592653589793238462643383279502;

  return 2.0+x*x*(1.0-x);
  }


// ------------------- FOR SHOOTING -------------------------

//  we write 
// -y'' + q(x)y(x) = f(x)
// as
// y[0]'=y[1]
// y[1]'= q(x) y[0] -f(x)
void for_rhs(double x, double y[DIM], double dydx[DIM])
  {
  #if DIM==2
    dydx[0]=y[1];
    dydx[1]=q(x)*y[0]-f(x); 
  #else
    fprintf(stderr, "Function valid only for DIM=2 and DIM=%d (%s, %d)\n", DIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  #endif 
  }


// perform the shoot
double func_for_shooting(double y0prime, // y'(a) 
                         double a,
                         double b,
                         int nsteps,
                         double **sol)  // sol[nsteps+1][2] solution of the problem if sol[nsteps][0]==0
  {
  double y0[2];

  y0[0]=0.0;
  y0[1]=y0prime;

  RK4(for_rhs, a, y0, b, nsteps, sol);                

  return sol[nsteps][0];
  }


typedef struct {
  double a;
  double b;
  int nsteps; 
  double **sol;
} Params;


// function to be used when only y0prime is important
double func_for_secant_with_context(double y0prime, void *ctx)
  {
  Params *p=ctx;
  
  return func_for_shooting(y0prime, p->a, p->b, p->nsteps, p->sol);
  }


// ---------------------------


int main(int argc, char **argv)
  {
  int i, nsteps; 
  const double a=0.0;
  const double b=1.0;
  double h, y0prime;
  char datafile[STRING_LENGTH];
  FILE *fp;
  Params ctx;

  double *sol_fde; // sol_fde[nsteps] solution
  double **sol_shoot; // sol_shoot[nsteps+1][2]

  // check input from command line
  if(argc != 2)
    {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s nsteps\n\n", argv[0]);
    fprintf(stdout, "  nsteps = number of steps in FDE or integration steps\n");

    return EXIT_SUCCESS;
    }
  else
    {  
    // read input value
    nsteps=atoi(argv[1]);
    }

  // allocate sol_fde
  sol_fde=(double *)malloc((unsigned long int)(nsteps)*sizeof(double));
  if(sol_fde == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  // allocate sol_shoot 
  sol_shoot=(double **)malloc((unsigned long int)(nsteps+1)*sizeof(double*));
  if(sol_shoot == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  for(i=0; i<=nsteps; i++)
     {
     sol_shoot[i]=(double *)malloc((unsigned long int)(DIM)*sizeof(double));
     if(sol_shoot[i] == NULL)
       {
       fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
       return EXIT_FAILURE;
       }
     }

  // ----------- FINITE DIFFERENCE METHOD ------------ 

  // initialize the solution to zero
  for(i=0; i<nsteps; i++)
     {
     sol_fde[i]=0.0;
     }

  // finite difference solution of the two point boundary value problem
  tpbvp_fdm(q, f, a, b, nsteps, sol_fde, 1.0e-8, 1000); 

  //print results
  sprintf(datafile, "ris_fdm_%d.dat", nsteps);

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  h=(b-a)/(double)(nsteps+1);
  for(i=0; i<nsteps; i++)
     {
     fprintf(fp,"%.12lf %.12lf\n", a+h*(double)(i+1), sol_fde[i]);
     }

  fclose(fp);

  // ----------- SHOOTING METHOD ------------ 

  y0prime=1;

  ctx.a=a;
  ctx.b=b;
  ctx.nsteps=nsteps;
  ctx.sol=sol_shoot;

  y0prime=secant_ctx(func_for_secant_with_context,
                     (void*)&ctx,
                     y0prime,
                     y0prime*1.2,
                     1.0e-8,
                     1000); 

  //print results
  sprintf(datafile, "ris_shoot_%d.dat", nsteps);

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  h=(b-a)/(double)(nsteps);
  for(i=0; i<nsteps; i++)
     {
     fprintf(fp,"%.12lf %.12lf\n", a+h*(double)i, sol_shoot[i][0]);
     }

  fclose(fp);

  // ------------------------------
 
  free(sol_fde);
  for(i=0; i<=nsteps; i++)
     {
     free(sol_shoot[i]);
     }
  free(sol_shoot);

  return EXIT_SUCCESS;
  }


