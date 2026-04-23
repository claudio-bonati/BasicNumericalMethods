#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define STRING_LENGTH 50

#include"../include/ode_tpbvp.h"   // <---- DIM is #defined here!

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


int main(int argc, char **argv)
  {
  int i, nsteps; 
  const double a=0.0;
  const double b=1.0;
  double h;
  char datafile[STRING_LENGTH];
  FILE *fp;

  double *sol; // sol[nsteps] solution

  // check input from command line
  if(argc != 2)
    {
    fprintf(stdout, "How to use this program:\n");
    fprintf(stdout, "  %s nsteps\n\n", argv[0]);
    fprintf(stdout, "  nsteps = number of integration steps\n");

    return EXIT_SUCCESS;
    }
  else
    {  
    // read input value
    nsteps=atoi(argv[1]);
    }

  // allocate y
  sol=(double *)malloc((unsigned long int)(nsteps)*sizeof(double));
  if(sol == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  // ----------- FINITE DIFFERENCE METHOD ------------ 

  h=(b-a)/(double)(nsteps+1);

  // finite difference solution of the two point boundary value problem
  tpbvp_fdm(q, f, a, b, nsteps, sol); 

  //print results
  sprintf(datafile, "ris_fdm_%d.dat", nsteps);

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }


  for(i=0; i<nsteps; i++)
     {
     fprintf(fp,"%.12lf %.12lf\n", a+h*(double)(i+1), sol[i]);
     }

  fclose(fp);

  // ------------------------

  free(sol);

  return EXIT_SUCCESS;
  }


