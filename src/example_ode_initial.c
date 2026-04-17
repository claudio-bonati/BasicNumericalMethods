#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define STRING_LENGTH 50

#include"../include/ode.h"


// Define the system: dy/dt = f(t, y)
//
// Example system for DIM==2: 
// dy0/dt = y1
// dy1/dt = -y0   
//
// y0=p  ;  y1=x
// H=p^2/2 + x^2/2  
void func(double t,         // in case the r.h.s also depends on time
          double y[DIM],    // input
          double dydt[DIM]) // output
  {
  (void)t;  // just to avoid compile time warning 
  #if DIM==2
    dydt[0] = y[1];
    dydt[1] = -y[0];
  #else
    fprintf(stderr, "Function valid only for DIM=2 and DIM=%d (%s, %d)\n", DIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  #endif 
  }


// Define dp/dt = fp(t, y) for Hamiltonian system, 
// with p=y[0,...,DIM/2-1], q=y[DIM/2, ..., DIM-1]
//
// Example system: 
// dy0/dt = y1
//
// y0=p  ;  y1=x
// H=p^2/2 + x^2/2  
void funcp(double t,         // in case the r.h.s also depends on time
           double y[DIM],    // input
           double dydt[DIM]) // output
  {
  (void)t;  // just to avoid compile time warning 
  #if DIM==2
    dydt[0] = y[1];
  #else
    fprintf(stderr, "Function valid only for DIM=2 and DIM=%d (%s, %d)\n", DIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  #endif
  }


// Define dq/dt = fp(t, y) for Hamiltonian system, 
// with p=y[0,...,DIM/2-1], q=y[DIM/2, ..., DIM-1]
//
// Example system: 
// dy1/dt = -y0
//
// y0=p  ;  y1=x
// H=p^2/2 + x^2/2  
void funcq(double t,         // in case the r.h.s also depends on time
           double y[DIM],    // input
           double dydt[DIM]) // output
  {
  (void)t;  // just to avoid compile time warning 
  #if DIM==2
    dydt[1] = -y[0];
  #else
    fprintf(stderr, "Function valid only for DIM=2 and DIM=%d (%s, %d)\n", DIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  #endif
  }

// -------------------------------------------------------

int main(void) 
  {
  double t0=0.0;
  double y0[DIM] = {1.0, 0.0};  // initial condition: y0=1, y1=0
  double time=20;
  double energy;
  int step, nsteps = 100; 
  char datafile[STRING_LENGTH];
  FILE *fp;

  double **y; // y[nsteps+1][DIM] solution

  // allocate y
  y=(double **)malloc((unsigned long int)(nsteps+1)*sizeof(double*));
  if(y == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    return EXIT_FAILURE;
    }
  for(step=0; step<=nsteps; step++)
     {
     y[step]=(double *)malloc((unsigned long int)(DIM)*sizeof(double));
     if(y[step] == NULL)
       {
       fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
       return EXIT_FAILURE;
       }
     }

  // ----------------------------

  // solve the system with Euler integration scheme
  Euler(&func, t0, y0, time, nsteps, y); 

  //print results
  strcpy(datafile, "ris_file_euler.dat");

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  for(step=0; step<=nsteps; step++)
     {
     energy=(y[step][0]*y[step][0]+y[step][1]*y[step][1])/2.0;

     fprintf(fp,"%.12lf ", t0+(double)step *(time-t0)/(double)nsteps);
     fprintf(fp,"%.12lf %.12lf ", y[step][0], y[step][1]);
     fprintf(fp,"%.12lf\n", energy);
     }

  fclose(fp);

  // ------------------------------

  // solve the system with Runge-Kutta 4 scheme
  RK4(&func, t0, y0, time, nsteps, y); 

  //print results
  strcpy(datafile, "ris_file_rk4.dat");

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  for(step=0; step<=nsteps; step++)
     {
     energy=(y[step][0]*y[step][0]+y[step][1]*y[step][1])/2.0;

     fprintf(fp,"%.12lf ", t0+(double)step *(time-t0)/(double)nsteps);
     fprintf(fp,"%.12lf %.12lf ", y[step][0], y[step][1]);
     fprintf(fp,"%.12lf\n", energy);
     }

  fclose(fp);

  // ------------------------------

  // solve the system with symplectic Euler integration scheme
  SympEuler(&funcp, &funcq, t0, y0, time, nsteps, y); 

  //print results
  strcpy(datafile, "ris_file_symp_euler.dat");

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  for(step=0; step<=nsteps; step++)
     {
     energy=(y[step][0]*y[step][0]+y[step][1]*y[step][1])/2.0;

     fprintf(fp,"%.12lf ", t0+(double)step *(time-t0)/(double)nsteps);
     fprintf(fp,"%.12lf %.12lf ", y[step][0], y[step][1]);
     fprintf(fp,"%.12lf\n", energy);
     }

  fclose(fp);

  // ------------------------------

  // solve the system with leapfrog
  leapfrog(&funcp, &funcq, t0, y0, time, nsteps, y); 

  //print results
  strcpy(datafile, "ris_file_symp_leapfrog.dat");

  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    return EXIT_FAILURE;
    }

  for(step=0; step<=nsteps; step++)
     {
     energy=(y[step][0]*y[step][0]+y[step][1]*y[step][1])/2.0;

     fprintf(fp,"%.12lf ", t0+(double)step *(time-t0)/(double)nsteps);
     fprintf(fp,"%.12lf %.12lf ", y[step][0], y[step][1]);
     fprintf(fp,"%.12lf\n", energy);
     }

  fclose(fp);



  // ------------------------

  // free y
  for(step=0; step<=nsteps; step++)
     {
     free(y[step]);
     }
  free(y);

  return EXIT_SUCCESS;
  }


