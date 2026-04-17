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
void f(double t,         // in case the r.h.s also depends on time
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
void fp(double t,         // in case the r.h.s also depends on time
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
void fq(double t,         // in case the r.h.s also depends on time
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
  int nsteps = 500; 

  char datafile[STRING_LENGTH];

  // solve the system with Euler integration scheme
  strcpy(datafile, "ris_file_euler.dat");
  Euler(&f, t0, y0, time, nsteps, datafile); 

  // solve the system with symplectic Euler integration scheme
  strcpy(datafile, "ris_file_symp_euler.dat");
  SympEuler(&fp, &fq, t0, y0, time, nsteps, datafile); 

  return EXIT_SUCCESS;
  }


