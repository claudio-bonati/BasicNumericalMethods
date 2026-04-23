#ifndef ODE_INIT_H
#define ODE_INIT_H

#define DIM 2  // Dimension of the system

// Euler algorith to solve 
// the system: dy/dt = func(t, y)
void Euler(void (*func)(double, double*, double*), // func(t,y,dy/dx) assigns dy/dx=f(t,y)
           double t0,        // initial time
           double y0[DIM],   // initial position
           double tend,      // integrate from t0 to tend
           int nsteps,       // using n steps of identical size
           double **y);      // y[steps+1][DIM] solution (y[0]=initial value)


// Runge-Kutta 4 algorith to solve 
// the system: dy/dt = func(t, y)
void RK4(void (*func)(double, double*, double*), // func(t,y,dy/dx) assigns dy/dx=f(t,y)
         double t0,        // initial time
         double y0[DIM],   // initial position
         double tend,      // integrate from t0 to tend
         int nsteps,       // using n steps of identical size
         double **y);      // y[steps+1][DIM] solution (y[0]=initial value)



// Symplectic Euler algorith to solve 
// the Hamiltonian system: dy/dt = func(y)
void SympEuler(void (*funcp)(double*, double*), // funcp(t,y,dy/dx) assigns dp/dt=f(y)
               void (*funcq)(double*, double*), // funcq(t,y,dy/dx) assigns dq/dt=f(y)
               double t0,       // initial time
               double y0[DIM],  // initial position
               double tend,     // integrate from t0 to tend
               int nsteps,      // using n steps of identical size
               double **y);     // y[steps+1][DIM] solution (y[0]=initial value)



// Leapfrog (aka Verlet) algorith to solve 
// the Hamiltonian system: dy/dt = func(y)
// with p=y[0,...,DIM/2-1], q=y[DIM/2, ..., DIM-1]
void leapfrog(void (*funcp)(double*, double*), // funcp(t,y,dy/dx) assigns dp/dt=f(y)
              void (*funcq)(double*, double*), // funcq(t,y,dy/dx) assigns dq/dt=f(y)
              double t0,       // initial time
              double y0[DIM],  // initial position
              double tend,     // integrate from t0 to tend
              int nsteps,      // using n steps of identical size
              double **y);     // y[steps+1][DIM] solution (y[0]=initial value)


#endif
