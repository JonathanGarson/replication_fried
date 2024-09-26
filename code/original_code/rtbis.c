/*Modified from Numerical Recipes in C (Press et al. 1988)*/

#include "header.h"
#define JMAX 400


double rtbis(double (*func)(double ka_kp, double Omega, double rho, double tfp), 
		double x1, double x2, double  xacc, double Omega, double rho, double tfp)
{
	int j; 
	double dx, f, fmid, xmid, rtb; 

	f = (*func)(x1, Omega, rho, tfp); 
	fmid = (*func)(x2, Omega, rho, tfp); 
	if (f*fmid >=0.0) printf("Root is not bracketed\n"); 
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 -x2, x2);
	for (j = 0; j< JMAX; j++) { 
		fmid = (*func)(xmid = rtb +(dx *=0.5), Omega, rho, tfp); 
		if (fmid <= 0.0) rtb = xmid; 
		if (fabs(dx) < xacc || fmid ==0.0) return rtb;
	}

	printf("Error: too many bisections\n"); 
	return 0.0; 

}
