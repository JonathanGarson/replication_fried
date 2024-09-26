/*Modified from Numerical Recipes in C (Press et al. 1988)*/

#include "header.h"

extern double mmax, hpmin; 
double amotry_r(double p[3][2], double y[], double psum[], int ndim, 
		double (*funk)(double[], double m, double hs, double ah,  int zi, double rho, double Omega, double  w,
			double hra_hr, double ph, int own,  int qi, int na, int pr),
		int ihi, double fac, double m, double hs,  int zi, double rho, double Omega, double w, double hra_hr,
		double ph, double ah, int own, int qi, int na,  int pr)
{
	int j; 
	double fac1, fac2, ytry, *ptry, lb[ndim], ub[ndim], rc, xprime;
	rc = 1.0e15;
	ptry = (double *)malloc(sizeof(double)*ndim); 
	fac1 = (1.0 - fac)/ndim; 
	fac2 = fac1 - fac; 

	lb[0] =0.0;
	lb[1] = 0.0; 

	//printf("aahere1\n"); 
	for (j =0; j <ndim; j++) {
		ptry[j]= psum[j]*fac1 - p[ihi][j]*fac2;
		ptry[j] = ptry[j]<lb[j] ? lb[j] : ptry[j];
	}
	ytry = (*funk)(ptry, m, hs, ah, zi, rho, Omega,  w, hra_hr, ph, own,  qi, na,   pr); 
	if (ytry <y[ihi]) {
		y[ihi] = ytry; 
		for (j = 0; j <ndim; j++){
			psum[j] +=ptry[j] -p[ihi][j];
			p[ihi][j] = ptry[j]; 
		}


	}
	free(ptry); 
	return ytry; 
}
