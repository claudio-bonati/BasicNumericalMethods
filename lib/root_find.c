#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define DEBUG

// root finding by bisection
double bisection(double (*func)(double),  // we want to solve func(x)=0, the function return the solution 
                 double x1,      
                 double x2,   // interval [x1, x2] such that func(x1)*func(x2)<0 
                 double xacc, // absolute accuracy on the position of the solution
                 int maxiter) // maximum number of iterations
  {
  int i;
  double f,fmid, xmid, dx, negx;
  
  f=(*func)(x1);
  fmid=(*func)(x2);

  if(f*fmid >= 0.0) 
    {
    fprintf(stderr, "Root must be bracketed for bisection to work! (%s, %d)\n", __FILE__, __LINE__);
    return 1.0/0.0;
    }

  // negx is the side of the integral where func<0
  if(f<0.0)
    {
    negx=x1;
    dx=x2-x1;
    }
  else
    {
    negx=x2;
    dx=x1-x2;
    }
 
  for(i=0; i<maxiter; i++)
     {
     dx*=0.5;
     xmid=negx+dx;
     fmid=(*func)(xmid);
     if(fmid <= 0.0) negx=xmid;
 
     #ifdef DEBUG
     printf("i=%d ; xmid=%g ; dx=%g ; fmid=%g\n", i, xmid, dx, fmid);
     #endif
     
     if(fabs(dx) < xacc || fmid == 0.0) return xmid;  // final result
     }

  fprintf(stderr, "Too many bissections in (%s, %d)\n", __FILE__, __LINE__);
  fprintf(stderr, "maxiter=%d reached\n", maxiter);

  return 1.0/0.0;
  }


// root finding by secant method
double secant(double (*func)(double),  // we want to solve func(x)=0, the function return the solution 
              double x1,      
              double x2,   // x1, x2 are the intitial two points
              double xacc, // absolute accuracy on the position of the solution
              int maxiter) // maximum number of iterations
  {
  int i;
  double  f1, f2, xnew, xold, dx, fnew, fold;

  f1=(*func)(x1);
  f2=(*func)(x2);

  //Pick the bound with the smaller function value as the most recent guess.
  if(fabs(f1) < fabs(f2)) 
    {
    xnew=x1;
    xold=x2;
    fnew=f1;
    fold=f2;
    } 
  else 
    {
    xnew=x2;
    xold=x1;
    fnew=f2;
    fold=f1;
    }

  for(i=0; i<maxiter; i++) 
     {
     dx=(xold-xnew)*fnew/(fnew-fold); // xnew+dx is solution of fnew + (fold-fnew)/(xold-xnew)(x-x_new)=0
     xold = xnew;
     xnew = xnew + dx;
     fold = fnew;
     fnew=(*func)(xnew);

     #ifdef DEBUG
     printf("i=%d ; xnew=%g ; dx=%g ; fnew=%g\n", i, xnew, dx, fnew);
     #endif

     if (fabs(dx) < xacc || fnew == 0.0) return xnew;
     }

  fprintf(stderr, "Too many iterations in (%s, %d)\n", __FILE__, __LINE__);
  fprintf(stderr, "maxiter=%d reached\n", maxiter);

  return 1.0/0.0;
  }


// root finding by secant method with context pointer
double secant_ctx(double (*func)(double, void*),  // we want to solve func(x)=0, the function return the solution 
                  void *ctx,   // the contex pointer used to pass other parameters to func
                  double x1,      
                  double x2,   // x1, x2 are the intitial two points
                  double xacc, // absolute accuracy on the position of the solution
                  int maxiter) // maximum number of iterations
  {
  int i;
  double  f1, f2, xnew, xold, dx, fnew, fold;

  f1=(*func)(x1, ctx);
  f2=(*func)(x2, ctx);

  printf("%.12lf %.12lf %.12lf %.12lf\n", x1, f1, x2, f2);

  //Pick the bound with the smaller function value as the most recent guess.
  if(fabs(f1) < fabs(f2)) 
    {
    xnew=x1;
    xold=x2;
    fnew=f1;
    fold=f2;
    } 
  else 
    {
    xnew=x2;
    xold=x1;
    fnew=f2;
    fold=f1;
    }

  for(i=0; i<maxiter; i++) 
     {
     dx=(xold-xnew)*fnew/(fnew-fold); // xnew+dx is solution of fnew + (fold-fnew)/(xold-xnew)(x-x_new)=0
     xold = xnew;
     xnew = xnew + dx;
     fold = fnew;
     fnew=(*func)(xnew, ctx);

     #ifdef DEBUG
     printf("i=%d ; xnew=%.12g ; dx=%.12g ; fnew=%.12g\n", i, xnew, dx, fnew);
     #endif

     if (fabs(dx) < xacc || fnew == 0.0) return xnew;
     }

  fprintf(stderr, "Too many iterations in (%s, %d)\n", __FILE__, __LINE__);
  fprintf(stderr, "maxiter=%d reached\n", maxiter);

  return 1.0/0.0;
  }


#undef DEBUG
