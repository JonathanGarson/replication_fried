/*Modified from Numerical Recipes in C (Press et al. 1988)*/

# include "header.h"

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
	int klo, khi, k; 
	double h, b, a; 

	klo =0; 
	khi = n-1;

	while(khi - klo > 1) {
		k = (khi + klo) >>1; 
		if (xa[k] > x) khi = k; 
		else klo = k; 
	}

	h = xa[khi] - xa[klo]; 
	if (h==0.0) printf("problem: h=0i\n");
	a = (xa[khi]-x)/h; 
	b =(x - xa[klo])/h; 
	*y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	
//	printf("%f\n", y2a[n-1]); 
			
 //	 printf("%i %i %f %f %f %f %f %f\n", klo, khi, xa[klo], x, xa[khi], ya[klo], *y, ya[khi]); 
}


