/*Modified from Numerical Recipes in C (Press et al. 1988)*/

#include "header.h"
#define TINY 1.0e-14
#define NMAX 5000
#define GET_PSUM \
	for(j =0; j<ndim;j++) {\
		for (sum = 0.0, i = 0; i<mpts; i++) sum +=p[i][j];\
		psum[j]=sum; }
#define SWAP(a,b) {swap = (a); (a)=(b); b = swap;}

extern double mmax, hpmin;

void amoeba_r(double p[3][2], double y[], int ndim, double ftol, 
		double (*funk)(double[], double m, double hs, double ah,  int zi, double rho, double Omega, double w,
			double hra_hr, double ph,  int own, int qi, int na, int pr),
		int *nfunk,  double m, double hs, int zi, double rho, double Omega, double w, double hra_hr, double ph, double ah,
		int own,int qi, int na,  int pr)
{


	int i, ihi, ilo, inhi, j, mpts=ndim+1; 
	double rtol, sum, swap, ysave, ytry, *psum, lb[ndim], ub[ndim], rc, xprime;

	psum = (double *)malloc(sizeof(double)*ndim); 
	*nfunk =0;

	lb[0] =0.0; 
	lb[1] =0.0; 

	GET_PSUM
		for (;;) {
		
			ilo = 0;
			ihi = y[0]>y[1] ? (inhi = 1,0) : (inhi = 0, 1); 
			for (i = 0; i <mpts; i++){
				if (y[i] <= y[ilo]) ilo = i; 
				if (y[i] > y[ihi]){
					inhi = ihi; 
					ihi =i; 
				} else if (y[i] > y[inhi] && (i !=  ihi)) inhi = i; 
			}

			rtol = 2.0*fabs(y[ihi] - y[ilo])/(fabs(y[ihi]) + fabs(y[ilo]) + TINY); 
			if ((rtol < ftol)|| (*nfunk > NMAX)) {
				SWAP(y[0], y[ilo])
					for (i = 0; i < ndim; i++) SWAP(p[0][i], p[ilo][i])
						break; 
			}

			*nfunk+=2; 
			ytry = amotry_r(p, y, psum, ndim, funk, ihi, -1.0, m, hs, zi, rho, Omega,  w, hra_hr, ph, ah, own, qi,na,  pr);
			if(ytry <=y[ilo])
				ytry = amotry_r(p, y, psum, ndim, funk, ihi, 2.0, m, hs, zi, rho, Omega,w, hra_hr, ph, ah, own, qi, na, pr); 
			else if (ytry >= y[inhi]){
				ysave = y[ihi];
				ytry = amotry_r(p, y, psum, ndim, funk, ihi, 0.5, m, hs, zi,rho, Omega, w, hra_hr, ph, ah, own, qi,na,  pr); 
				if (ytry >=ysave) {
					for ( i =0; i< mpts; i++){
						if (i != ilo) {
							for (j = 0; j<ndim; j++){
								p[i][j]= 0.5*(p[i][j] + p[ilo][j]);
								p[i][j]= p[i][j] < lb[j] ?  lb[j] : p[i][j];
								psum[j] = p[i][j]; 
							}
							y[i] = (*funk)(psum, m, hs, ah, zi, rho, Omega, w, hra_hr, ph, own, qi, na,  pr);
						}
					}
					*nfunk +=ndim; 
					GET_PSUM
				}
			} else --(*nfunk);

		}
	free(psum); 
}
