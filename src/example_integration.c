#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/integration.h"

/*
Test for integration algorithms
*/

double func(double x)
  {
  return x*x*x;
  }

int main(void)
  {
  double x1, x2, ris, accuracy;
  int n;

  const int maxiter=50;

  printf("\n");
  printf("Integration using trapezoidal rule with fixed number of points\n");
  x1=0.0;
  x2=1.0;

  n=5; // 2^(5) points
  ris=trapezoidal(&func, x1, x2, n);
  printf("n=%d ; ris=%.15lf ; |ris-trueris|=%.15lf\n", n, ris, fabs(ris-1.0/4.0));
  n=6; // 2^(6) points
  ris=trapezoidal(&func, x1, x2, n);
  printf("n=%d ; ris=%.15lf ; |ris-trueris|=%.15lf\n", n, ris, fabs(ris-1.0/4.0));
  n=10; // 2^(10) points
  ris=trapezoidal(&func, x1, x2, n);
  printf("n=%d ; ris=%.15lf ; |ris-trueris|=%.15lf\n", n, ris, fabs(ris-1.0/4.0));
  n=11; // 2^(11) points
  ris=trapezoidal(&func, x1, x2, n);
  printf("n=%d ; ris=%.15lf ; |ris-trueris|=%.15lf\n", n, ris, fabs(ris-1.0/4.0));
  n=20; // 2^(20) points
  ris=trapezoidal(&func, x1, x2, n);
  printf("n=%d ; ris=%.15lf ; |ris-trueris|=%.15lf\n", n, ris, fabs(ris-1.0/4.0));

  printf("\n\n");

  printf("Integration using trapezoidal rule with a self-consistent number of points \n");

  accuracy=1.0e-6;
  ris=trapezoidal_self(&func, x1, x2, accuracy, maxiter); 
  printf("accuracy=%g  ris=%.15lf ; |ris-trueris|=%.15lf\n", accuracy, ris, fabs(ris-1.0/4.0));

  accuracy=1.0e-12;
  ris=trapezoidal_self(&func, x1, x2, accuracy, maxiter); 
  printf("accuracy=%g  ris=%.15lf ; |ris-trueris|=%.15lf\n", accuracy, ris, fabs(ris-1.0/4.0));

  printf("\n\n");

  printf("Integration using Richardson extrapolation of the trapezoidal rule with fixed number of points\n");
  n=5; // 2^(10) points
  ris=richardson_for_int(&func, x1, x2, n, &trapezoidal, 2);
  printf("n=%d ; ris=%.15lf ; |ris-trueris|=%.15lf\n", n, ris, fabs(ris-1.0/4.0));
  n=10; // 2^(15) points
  ris=richardson_for_int(&func, x1, x2, n, &trapezoidal, 2);
  printf("n=%d ; ris=%.15lf ; |ris-trueris|=%.15lf\n", n, ris, fabs(ris-1.0/4.0));

  return EXIT_SUCCESS;
  }

