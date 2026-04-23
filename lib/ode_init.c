#include<stdio.h>
#include<stdlib.h>

#include"../include/ode_init.h"  // <--- DIM is defined here as a macro


// Euler algorith to solve 
// the system: dy/dt = f(t, y)
void Euler(void (*func)(double, double*, double*), // func(t,y,dy/dx) assigns dy/dx=f(t,y)
           double t0,        // initial time
           double y0[DIM],   // initial position
           double tend,      // integrate from t0 to tend
           int nsteps,       // using n steps of identical size
           double **y)       // y[steps+1][DIM] solution
  {
  int i, step;
  double dydt[DIM];  // value of the r.h.s at time t
  double t=t0;
  double h=(tend-t0)/(double)nsteps; // step size

  for(i=0; i<DIM; i++)
     {
     y[0][i]=y0[i];
     }

  for(step=0; step<nsteps; step++)
     {
     // evaluate the r.h.s
     func(t, y[step], dydt);

     // update the components
     for(i=0; i<DIM; i++) 
        {
        y[step+1][i] = y[step][i] + h * dydt[i];
        }

     // increase time
     t=t0+(double)(step+1)*h;
     }
  }


// Runge-Kutta 4 algorith to solve 
// the system: dy/dt = f(t, y)
void RK4(void (*func)(double, double*, double*), // func(t,y,dy/dx) assigns dy/dx=f(t,y)
         double t0,        // initial time
         double y0[DIM],   // initial position
         double tend,      // integrate from t0 to tend
         int nsteps,       // using n steps of identical size
         double **y)       // y[steps+1][DIM] solution
  {
  int i, step;
  double dydt1[DIM], dydt2[DIM], dydt3[DIM], dydt4[DIM], yaux[DIM];
  double t=t0;
  double h=(tend-t0)/(double)nsteps; // step size

  for(i=0; i<DIM; i++)
     {
     y[0][i]=y0[i];
     }

  for(step=0; step<nsteps; step++)
     {
     // evaluate the r.h.s
     func(t, y[step], dydt1);  // f(t, y)
   
     for(i=0; i<DIM; i++)
        {
        yaux[i]=y[step][i]+h*dydt1[i]/2.0;
        }
     func(t+h/2.0, yaux, dydt2);  // f(t+h/2, y+k1*h/2)
    
     for(i=0; i<DIM; i++)
        {
        yaux[i]=y[step][i]+h*dydt2[i]/2.0;
        }
     func(t+h/2.0, yaux, dydt3);  // f(t+h/2, y+k2*h/2)
 
     for(i=0; i<DIM; i++)
        {
        yaux[i]=y[step][i]+h*dydt3[i];
        }
     func(t+h, yaux, dydt4);  // f(t+h, y+k3*h)
 
     // update the components
     for(i=0; i<DIM; i++) 
        {
        y[step+1][i] = y[step][i] + h*(dydt1[i]+2.0*dydt2[i]+2.0*dydt3[i]+dydt4[i])/6.0;
        }

     // increase time
     t=t0+(double)(step+1)*h;
     }
  }


// Symplectic Euler algorith to solve 
// the Hamiltonian system: dy/dt = f(y)
// with p=y[0,...,DIM/2-1], q=y[DIM/2, ..., DIM-1]
void SympEuler(void (*funcp)(double*, double*), // funcp(t,y,dy/dx) assigns dp/dt=f(y)
               void (*funcq)(double*, double*), // funcq(t,y,dy/dx) assigns dq/dt=f(y)
               double t0,       // initial time
               double y0[DIM],  // initial position
               double tend,     // integrate from t0 to tend
               int nsteps,      // using n steps of identical size
               double **y)       // y[steps+1][DIM] solution
  {
  int i, step;
  double dydt[DIM];  // value of the r.h.s at time t
  double h=(tend-t0)/(double)nsteps; // step size

  for(i=0; i<DIM; i++)
     {
     y[0][i]=y0[i];
     }

  for(step=0; step<nsteps; step++)
     {
     // evaluate the r.h.s for p
     funcp(y[step], dydt);

     // update the p components
     for(i=0; i<DIM/2; i++) 
        {
        y[step+1][i] = y[step][i] + h * dydt[i];
        }

     // evaluate the r.h.s for q
     funcq(y[step+1], dydt);

     // update the q components
     for(i=DIM/2; i<DIM; i++) 
        {
        y[step+1][i] = y[step][i] + h * dydt[i];
        }
     }
  }


// Leapfrog (aka Verlet) algorith to solve 
// the Hamiltonian system: dy/dt = func(y)
// with p=y[0,...,DIM/2-1], q=y[DIM/2, ..., DIM-1]
void leapfrog(void (*funcp)(double*, double*), //  funcp(t,y,dy/dx) assigns dp/dt=f(y)
              void (*funcq)(double*, double*), //  funcq(t,y,dy/dx) assigns dq/dt=f(y)
              double t0,       // initial time
              double y0[DIM],  // initial position
              double tend,     // integrate from t0 to tend
              int nsteps,      // using n steps of identical size
              double **y)      // y[steps+1][DIM] solution
  {
  int i, step;
  double dydt[DIM];  // value of the r.h.s at time t
  double h=(tend-t0)/(double)nsteps; // step size

  for(i=0; i<DIM; i++)
     {
     y[0][i]=y0[i];
     }

  for(step=0; step<nsteps; step++)
     {
     // evaluate the r.h.s for p
     funcp(y[step], dydt);

     // update the p components
     for(i=0; i<DIM/2; i++) 
        {
        y[step+1][i] = y[step][i] + h*dydt[i]/2.0;
        }

     // evaluate the r.h.s for q
     funcq(y[step+1], dydt);

     // update the q components
     for(i=DIM/2; i<DIM; i++) 
        {
        y[step+1][i] = y[step][i] + h*dydt[i];
        }

     // evaluate the r.h.s for p
     funcp(y[step+1], dydt);

     // update the p components
     for(i=0; i<DIM/2; i++) 
        {
        y[step+1][i] = y[step+1][i] + h*dydt[i]/2.0;
        }
     }
  }


