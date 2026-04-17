#include<stdio.h>
#include<stdlib.h>

#include"../include/ode.h"  // <--- DIM is defined here as a macro

//#define DEBUG   // if defined the energy of the harmonic oscillator is also printed

// Euler algorith to solve 
// the system: dy/dt = func(t, y)
void Euler(void (*func)(double, double*, double*), // r.h.s 
           double t0,        // initial time
           double y0[DIM],   // initial position
           double tend,      // integrate from t0 to tend
           int nsteps,       // using n steps of identical size
           char *datafile)  // name of the output file
  {
  int i, step;
  double y[DIM];     // value of the solution at time t
  double dydt[DIM];  // value of the r.h.s at time t
  double t=t0;
  double h=(tend-t0)/(double)nsteps; // step size
  FILE *fp;

  // open data file
  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  for(i=0; i<DIM; i++)
     {
     y[i]=y0[i];
     }

  #ifdef DEBUG
    fprintf(fp, "%.12lf %.12lf %.12lf %.12lf\n", t, y[0], y[1], (y[0]*y[0]+y[1]*y[1])/2.0);
  #else
    fprintf(fp, "%.12lf %.12lf %.12lf\n", t, y[0], y[1]);
  #endif

  for(step=0; step<nsteps; step++)
     {
     // evaluate the r.h.s
     func(t, y, dydt);

     // update the components
     for(i=0; i<DIM; i++) 
        {
        y[i] = y[i] + h * dydt[i];
        }

     // increase time
     t=t0+(double)(step+1)*h;

     #ifdef DEBUG
       fprintf(fp, "%.12lf %.12lf %.12lf %.12lf\n", t, y[0], y[1], (y[0]*y[0]+y[1]*y[1])/2.0);
     #else
       fprintf(fp, "%.12lf %.12lf %.12lf\n", t, y[0], y[1]);
     #endif
     }

  // close datafile
  fclose(fp);
  }



// Symplectic Euler algorith to solve 
// the Hamiltonian system: dy/dt = func(t, y)
// with p=y[0,...,DIM/2-1], q=y[DIM/2, ..., DIM-1]
void SympEuler(void (*funcp)(double, double*, double*), // r.h.s for p 
               void (*funcq)(double, double*, double*), // r.h.s for q
               double t0,       // initial time
               double y0[DIM],  // initial position
               double tend,     // integrate from t0 to tend
               int nsteps,      // using n steps of identical size
               char *datafile)  // name of the output file
  {
  int i, step;
  double y[DIM];     // value of the solution at time t
  double dydt[DIM];  // value of the r.h.s at time t
  double t=t0;
  double h=(tend-t0)/(double)nsteps; // step size
  FILE *fp;

  // open data file
  fp=fopen(datafile, "w");
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", datafile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  for(i=0; i<DIM; i++)
     {
     y[i]=y0[i];
     }

  #ifdef DEBUG
    fprintf(fp, "%.12lf %.12lf %.12lf %.12lf\n", t, y[0], y[1], (y[0]*y[0]+y[1]*y[1])/2.0);
  #else
    fprintf(fp, "%.12lf %.12lf %.12lf\n", t, y[0], y[1]);
  #endif

  for(step=0; step<nsteps; step++)
     {
     // evaluate the r.h.s for p
     funcp(t, y, dydt);

     // update the p components
     for(i=0; i<DIM/2; i++) 
        {
        y[i] = y[i] + h * dydt[i];
        }

     // evaluate the r.h.s for q
     funcq(t, y, dydt);

     // update the q components
     for(i=DIM/2; i<DIM; i++) 
        {
        y[i] = y[i] + h * dydt[i];
        }

     // increase time
     t=t0+(double)(step+1)*h;

     #ifdef DEBUG
       fprintf(fp, "%.12lf %.12lf %.12lf %.12lf\n", t, y[0], y[1], (y[0]*y[0]+y[1]*y[1])/2.0);
     #else
       fprintf(fp, "%.12lf %.12lf %.12lf\n", t, y[0], y[1]);
     #endif
     }

  // close datafile
  fclose(fp);
  }

