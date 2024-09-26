/*Modified from Numerical Recipes in C (Press et al. 1988)*/

#include "header.h"

void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, int hs_start,  double x1, double x2, double *y)
{
	int j; 
	double *ytmp, *yytmp; 
	
	ytmp = (double *)malloc(sizeof(double)*m);
	yytmp = (double *)malloc(sizeof(double)*m);

	for (j =hs_start; j < hs_start+ m; j++){
		splint(x2a, ya[j], y2a[j], n, x2, &yytmp[j-hs_start]);
	}
	
	spline(x1a, yytmp, m, 1.0e30, 1.0e30, ytmp);
	splint(x1a, yytmp, ytmp, m, x1, y);

	free(ytmp); 
	free(yytmp); 
}
