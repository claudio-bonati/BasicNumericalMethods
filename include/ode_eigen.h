#ifndef ODE_EIGEN_H
#define ODE_EIGEN_H

// eigenvalue problem on [a, b] for a SCALAR equation (single component)
// -y'' + q(x)y(x) = lambda y(x)
// with y(a)=0, y(b)=0
//
// eigvals are sorted in such a way that eigvals[0]<=eigvals[1]<=eigvals[2]<=...
// and eigvects[:][i] is the eigenvector corresponfing to eigvals[i]
void eigen_fdm(double (*q)(double), 
               double a,
               double b,
               int nsteps,  // numbert of points of the discretization
               double *eigvals,
               double **eigvects, 
               double accuracy,    // accuracy and maxiter for Jacobi diagonalization
               int maxiter);

#endif
