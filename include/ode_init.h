#ifndef ODE_INIT_H
#define ODE_INIT_H

#define DIM 2  // Dimension of the system

// Euler algorith to solve 
// the system: dy/dt = func(t, y)
void Euler(void (*func)(double, double*, double*), // r.h.s 
           double t0,        // initial time
           double y0[DIM],   // initial position
           double tend,      // integrate from t0 to tend
           int nsteps,       // using n steps of identical size
           double **y);      // y[steps+1][DIM] solution


// Runge-Kutta 4 algorith to solve 
// the system: dy/dt = func(t, y)
void RK4(void (*func)(double, double*, double*), // r.h.s 
         double t0,        // initial time
         double y0[DIM],   // initial position
         double tend,      // integrate from t0 to tend
         int nsteps,       // using n steps of identical size
         double **y);      // y[steps+1][DIM] solution


// Symplectic Euler algorith to solve 
// the Hamiltonian system: dy/dt = func(y)
void SympEuler(void (*funcp)(double*, double*), // r.h.s for p 
               void (*funcq)(double*, double*), // r.h.s for q
               double t0,       // initial time
               double y0[DIM],  // initial position
               double tend,     // integrate from t0 to tend
               int nsteps,      // using n steps of identical size
               double **y);     // y[steps+1][DIM] solution


// Leapfrog (aka Verlet) algorith to solve 
// the Hamiltonian system: dy/dt = func(y)
// with p=y[0,...,DIM/2-1], q=y[DIM/2, ..., DIM-1]
void leapfrog(void (*funcp)(double*, double*), // r.h.s for p 
              void (*funcq)(double*, double*), // r.h.s for q
              double t0,       // initial time
              double y0[DIM],  // initial position
              double tend,     // integrate from t0 to tend
              int nsteps,      // using n steps of identical size
              double **y);     // y[steps+1][DIM] solution
 

#endif
