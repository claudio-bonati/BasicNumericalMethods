#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/root_find.h"

/*
Test to check root finding algorithms
*/

double func(double x)
  {
  return -(x-1.0)*(x-2.0)*(x+3.0);
  }

int main(void)
  {
  double x1, x2, ris;
  const int maxiter=50;
  const double xacc=1.0e-12;

  printf("solution with the bisection method\n");
  x1=0.8;
  x2=1.5;
  ris=bisection(&func, x1, x2, xacc, maxiter);
  printf("ris=%.15lf ; |ris-trueris|=%.15lf\n", ris, fabs(ris-1.0));

  printf("\n\n");

  printf("solution with the secant method\n");
  x1=0.8;
  x2=1.5;
  ris=secant(&func, x1, x2, xacc, maxiter);
  printf("ris=%.15lf ; |ris-trueris|=%.15lf\n", ris, fabs(ris-1.0));


  return EXIT_SUCCESS;
  }

