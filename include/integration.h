#ifndef INTEGRATION_H
#define INTEGRATION_H

// integration function with trapezoidal rule and fixed number of points
// elementary step
// int_c^d f(x)dx=(d-c)[f(d)-f(c)]/2
// asymptotic error usin N discretization points
// \int_a^b f(x)dx - trap = - (b-a)^2[f'(b)-f'(a)]/(12 N^2)
double trapezoidal(double (*func)(double),  // function to be integrated
                   double a, 
                   double b,  // integration interval [a,b] 
                   int n);     // 2^n points used for the discretization, n>=1

// integration function using trapezoidal rule and self-consistent number of points
double trapezoidal_self(double (*func)(double), // function to be integrated
                        double a, 
                        double b,  // integration interval [a,b]
                        double accuracy, // accuracy of the result 
                        int maxiter);  // 2^{maxiter} is the maximum number of points to be used

// Richardson extrapolation to improve the order of convergence
// of the integration algorithm
double richardson_for_int(double (*func)(double),  // function to be integrated
                          double a, 
                          double b,  // integration interval [a,b] 
                          int n,     // 2^n and 2^(n+1) points used for the discretization, n>=1
                          double (*intalg)(double (*)(double), double, double, int), // integration algorithm, e.g. trapezoidal
                          int order);  // order of the integration algorirhtm

#endif
