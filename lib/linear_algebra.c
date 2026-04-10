#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Gauss-Jordan elimination with full pivoting
// to solve Ax=b with A a matrix of size n
//
// After use A->identity and b->permutation of x
void GaussJordan_fullpivot(int n, double **A, double *b, double *x)
  {
  int i, j, k, pivotrow, pivotcol, tmpindex;
  double max, tmp, pivot;
  const double epsilon=1.0e-12;
  int *colindex; // to track column swaps

  colindex=(int *)malloc((unsigned long int)(n)*sizeof(int));
  if(colindex == NULL)
    {
    fprintf(stderr, "allocation problem (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialized to no-swap
  for(i=0; i<n; i++) 
     {
     colindex[i] = i;
     }
  
   for(k=0; k<n; k++)  // loop on columns
      {
      // find the pivot (max abs value) in the submatrix with indices >=k
      max=0.0;
      pivotrow=k;
      pivotcol=k;
  
      for(i=k; i<n; i++) 
         {
         for(j=k; j<n; j++) 
            {
            if(fabs(A[i][j]) > max) 
              {
              max = fabs(A[i][j]);
              pivotrow = i;
              pivotcol = j;
              }
            }  
         }
  
      if(fabs(max) < epsilon)
        {
        fprintf(stderr, "Matrix is (nearly) singular (%s, %d)\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        }
  
      // swap rows
      if(pivotrow != k)
        {
        for(j = 0; j < n; j++)
           {
           tmp = A[k][j];
           A[k][j] = A[pivotrow][j];
           A[pivotrow][j] = tmp;
           }
        tmp = b[k];
        b[k] = b[pivotrow];
        b[pivotrow] = tmp;
        }
  
      // swap columns
      if(pivotcol != k) 
        {
        for(i=0; i<n; i++) 
           {
           tmp = A[i][k];
           A[i][k] = A[i][pivotcol];
           A[i][pivotcol] = tmp;
          }
        tmpindex = colindex[k];
        colindex[k] = colindex[pivotcol];
        colindex[pivotcol] = tmpindex;
        }

      // now the pivot is in position [k][k]  

      // normalize pivot row
      pivot = A[k][k];
      for(j=0; j<n; j++) 
         {
         A[k][j] /= pivot;
         }
      b[k] /= pivot;
  
      // eliminate other rows
      for(i=0; i<n; i++) 
         {
         if(i != k)  // i = row index different from that of the pivot
           {
           tmp = A[i][k];
           for(j = 0; j < n; j++) 
              {
              A[i][j] -= tmp * A[k][j];
              }
           b[i] -= tmp * b[k];
           }
         }
      } // <- end of the loop on columns
  
  // reorder the solution due to column swaps
  for(i=0; i<n; i++) 
     {
     x[colindex[i]] = b[i];
     }

  free(colindex);
  }

