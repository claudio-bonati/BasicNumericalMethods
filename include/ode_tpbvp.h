#ifndef ODE_TPBVP_H
#define ODE_TPBVP_H


// two point boundary value problem on [a, b] for a SCALAR equation (single component)
// -y'' + q(x)y(x) = f(x)
// with y(a)=0, y(b)=0
void tpbvp_fdm(double (*q)(double), 
               double (*f)(double),
               double a,
               double b,
               int nsteps,   // numbert of points of the discretization
               double *sol); // nsteps vector: sol[i]=solution in x=a+(i+1)h, with h=(b-a)/(nsteps+1)
                             // sol[0] is the solution in a+h
                             // sol[nsteps-1] is the solution in a+nsteps*(b-a)/(nsteps+1)=b-h


#endif
