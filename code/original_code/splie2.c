/*Modified from Numerical Recipes in C (Press et al. 1988)*/

#include "header.h"


void splie2(double x1a[], double x2a[], double **ya, int m, int n, int hs_start,  double **y2a)
{


	int j; 

	for (j =hs_start; j < hs_start + m; j++)
		spline(x2a, ya[j], n, 1.0e30, 1.0e30,  y2a[j]);

//	printf("1)%i %f %f %f %f\n", hs1, ya[hs1+5][6], y2a[hs1+1][1], y2a[hs1+2][1], y2a[hs1+m-1][40]); 

}

