#include<math.h>
#include<stdio.h>
#include<stdlib.h>

//#define DEBUG

// integration function with trapezoidal rule and fixed number of points
// elementary step
// int_c^d f(x)dx=(d-c)[f(d)+f(c)]/2
// asymptotic error usin N discretization points
// \int_a^b f(x)dx - trap \approx - (b-a)^2[f'(b)-f'(a)]/(12 N^2)
double trapezoidal(double (*func)(double),  // function to be integrated
                   double a, 
                   double b,  // integration interval [a,b] 
                   int n)     // 2^n points used for the discretization, n>=1
  {
  double x, sum, delta;
  static double s;
  int numpoints, i;
  
  if(n<1)
    {
    fprintf(stderr, "n>0 has to be used for the trapezoidal rule (%s, %d)\n", __FILE__, __LINE__);
    return 1.0/0.0;
    }
  
  numpoints=1;
  for(i=0; i<n; i++) numpoints*=2;

  delta=(b-a)/(double) (numpoints-1);

  x=a;
  sum=0.5*(*func)(x);
  for(i=1; i<numpoints-1; i++) 
     {
     x+=delta;
     sum += (*func)(x);
     }
  sum+=0.5*(*func)(b);

  s=delta*sum;

  return s;
  }


// integration function using trapezoidal rule and self-consistent number of points
double trapezoidal_self(double (*func)(double), // function to be integrated
                        double a, 
                        double b,  // integration interval [a,b]
                        double accuracy, // accuracy of the result 
                        int maxiter)  // 2^{maxiter} is the maximum number of points to be used
  {
  int i;
  double s, olds;
  
  olds=trapezoidal(func, a, b, 4);

  for(i=5; i<=maxiter; i++)
     {
     s=trapezoidal(func, a, b, i);

     #ifdef DEBUG
     printf("i=%d ; s=%.15lf\n", i, s);
     #endif
    
     if(fabs(s-olds) < accuracy) return s;

     olds=s;
     }

  fprintf(stderr, "maxiter=%d reached (%s, %d)\n", maxiter, __FILE__, __LINE__);
  return 1.0/0.0;
  }


// Richardson extrapolation for an integration algorithm of order "order", i.e. with an error 
// which scales as 1/N^{order}, with N the number of discretization points
//
// The used combination cancels the leading term of the error
double richardson_for_int(double (*func)(double),  // function to be integrated
                          double a, 
                          double b,  // integration interval [a,b] 
                          int n,     // 2^n and 2^(n+1) points used for the discretization, n>=1
                          double (*intalg)(double (*)(double), double, double, int), // integration algorithm, e.g. trapezoidal
                          int order)  // order of the integration algorirhtm
  {
  double ris1, ris2, ris;

  ris1=(*intalg)(func, a, b, n);
  ris2=(*intalg)(func, a, b, n+1);

  ris=(pow(2.0, order)*ris2-ris1) / (pow(2.0, order) - 1.0);

  #ifdef DEBUG
  printf("n1=%d ; n2=%d\n", n, n+1);
  printf(" ris1=%.15lf\n", ris1); 
  printf(" ris2=%.15lf\n", ris2); 
  printf("  ris=%.15lf\n", ris); 
  #endif

  return ris;
  }


#undef DEBUG
