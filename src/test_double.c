#include<stdio.h>
#include<stdlib.h>

int main(void)
    {
    double x, y;
    int iter;

    printf("Test 1: 1/2^n > 0\n");
    iter=0;
    x=1.0;
    while(x>0)
         {
         iter++;
         x=x/2.0;
         printf("%d %.16lf %.16g\n", iter, x, x);
         }

    printf("\n\n\n");

    printf("Test 2: 1 + 1/2^n > 1\n");
    iter=0;
    x=1.0;
    y=2.0;
    while(y>1.0)
         {
         iter++;
         x=x/2.0;
         y=1.0+x;
         printf("%d %.16lf %.16g\n", iter, y, y);
         }

    return EXIT_SUCCESS;
    }

