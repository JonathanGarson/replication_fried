/*Modified from Numerical Recipes in C (Press et al. 1988)*/

#include "header.h"
#define R 0.61803399
#define C (1.0 - R)
#define SHIFT2(a,b,c) (a)=(b); (b)=(c);
#define SHIFT3(a,b,c,d) (a)=(b); (b)=(c); (c)=(d); 

void  golden (double rho, double ax, double bx, double cx, 
	double(*f)(double rho, double x, int hpi, int hai,  int ai, int zi, int ub), 
	double tol, double *xmin, int hpi, int hai, int ai, int zi)
{
	double f1, f2, x0, x1, x2, x3; 

	x0 = ax; 
	x3 = cx; 
	if(fabs(cx-bx) > fabs(bx-ax)){
		x1=bx; 
		x2=bx+C*(cx-bx); 
	} else {
		x2=bx; 
		x1 = bx-C*(bx-ax); 
	}
	f1 = (*f)(rho, x1, hpi, hai, ai, zi, 0); 
	f2 = (*f)(rho, x2, hpi, hai, ai, zi, 0); 
	
	while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))){
		if (f2<f1){
			SHIFT3(x0, x1, x2, R*x1 + C*x3)
			SHIFT2(f1, f2, (*f)(rho, x2, hpi, hai, ai, zi, 0))
		} else {
			SHIFT3(x3, x2, x1, R*x2 + C*x0)
			SHIFT2(f2, f1, (*f)(rho, x1, hpi, hai, ai, zi, 0))
		}
	}
	if(f1 < f2){
		*xmin = x1; 
		//return f1; 
	}else {
		*xmin = x2; 
		//return f2; 
	}
}

