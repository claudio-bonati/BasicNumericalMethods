#ifndef ROOT_FIND_H
#define ROOT_FIND_H

// root finding by bisection
double bisection(double (*func)(double),  // we want to solve func(x)=0, the function return the solution 
                 double x1,      
                 double x2,   // interval [x1, x2] such that func(x1)*func(x2)<0 
                 double xacc, // absolute accuracy on the position of the solution
                 int maxiter); // maximum number of iterations

double secant(double (*func)(double),  // we want to solve func(x)=0, the function return the solution 
              double x1,      
              double x2,   // x1, x2 are the intitial two points
              double xacc, // absolute accuracy on the position of the solution
              int maxiter); // maximum number of iterations
 
#endif
