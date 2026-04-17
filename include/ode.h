#ifndef ODE_H
#define ODE_H

#define DIM 2  // Dimension of the system

// Euler algorith to solve 
// the system: dy/dt = func(t, y)
void Euler(void (*func)(double, double*, double*), // r.h.s 
           double t0,        // initial time
           double y0[DIM],   // initial position
           double tend,      // integrate from t0 to tend
           int nsteps,       // using n steps of identical size
           char *datafile);  // name of the output file


// Symplectic Euler algorith to solve 
// the Hamiltonian system: dy/dt = func(t, y)
void SympEuler(void (*funcp)(double, double*, double*), // r.h.s for p 
               void (*funcq)(double, double*, double*), // r.h.s for q
               double t0,       // initial time
               double y0[DIM],  // initial position
               double tend,     // integrate from t0 to tend
               int nsteps,      // using n steps of identical size
               char *datafile); // name of the output file


#endif
