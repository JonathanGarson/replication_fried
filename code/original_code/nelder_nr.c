/*
Program: nelderi_nr.c
Purpose: Shell program for numerical recipies nelder-meade routine. For the no-risk specification for homeowners. One less choice variable because homeowners do not need to puchase insurance in this specification 
Calls: ucont.c amoeba_nr.c amoeba_r.c 
*/

#include "header.h"
				
extern double ftol,  theta,   psi_h, eta, hpmin, theta, *mgrid, *hsgrid; 
				
void nelder_nr(double prime1[], double prime2[],  double ytemp, int ndim, double *hpprime, double *haprime, double *aprime, 
		double *xprime, double *hsrprime, double *ystar, double m, double hs, int zi, double rho, double Omega, 
		double w, double ar, double ph, double ah,  int own, int qi, int na, int pr) 

{
	int nfunc,j, i, k, counter, hpsign[8], hasign[8], asign[8],  new_vals, pr_ucont,  pr_amoeba, max_loops;
	double ytemp1, hpprime1, hpprime2, haprime1, haprime2, aprime1, aprime2, xprime1, xprime3,  xprime2,  d[ndim], y[ndim+1], 
	       p[ndim+1][ndim], hpprime_temp, haprime_temp, aprime_temp, F,   hptprime, c, hpprime3, 
	       haprime3, aprime3,  hsrprime1, hsrprime2, hsrprime_temp, xprime_temp;

	new_vals =0; 

	hpsign[0] = 1;  hpsign[1] = -1; hpsign[2] = 1; hpsign[3] = 1; hpsign[4]= -1; hpsign[5] = -1; hpsign[6] =1; hpsign[7] =-1;
	hasign[0] = 1; hasign[1] = 1; hasign[2] =-1; hasign[3] = 1; hasign[4]= -1; hasign[5] = 1; hasign[6] =-1; hasign[7] =-1;
	asign[0] = 1; asign[1] = 1; asign[2] = 1; asign[3] = -1; asign[4]= 1; asign[5] = -1; asign[6] =-1; asign[7] =-1;

	if (own ==1){
		hpprime1 = prime1[0]; haprime1 = prime1[1]; aprime1 = prime1[2];   
		hpprime2 = prime2[0]; haprime2 = prime2[1]; aprime2 = prime2[2];  
	}else {
		hsrprime1 = prime1[0]; hsrprime2 = prime2[0]; 
		aprime1 = prime1[1]; aprime2 = prime2[1]; 
	}
	
	pr_amoeba =0; 
	pr_ucont =0; 
	max_loops =0; 
	if(own==1){
		for (counter =0; counter <1; counter++){

			p[0][0] = hpprime1; 
			p[0][1] = haprime1; 
			p[0][2] = aprime1;
			p[1][0] = hpprime2; 
			p[1][1] = haprime1; 
			p[1][2] = aprime1;
			p[2][0] = hpprime1; 
			p[2][1] = haprime2; 
			p[2][2] = aprime1;
			p[3][0] = hpprime1; 
			p[3][1] = haprime1; 
			p[3][2] = aprime2;

			//check that starting simplex is reasonable 
			for (i =0; i < ndim +1; i++) {

				hpprime3= p[i][0]; 
				haprime3 = p[i][1]; 
				aprime3 = p[i][2]; 
				
				if (hpprime3 < hpmin){
					hpprime3 = hpmin + 0.01*i*hpprime3; 
					p[i][0] = hpprime3; 
				}
				
			}	

			for (i = 0; i <ndim+1; i ++){
				for (k = 0; k < ndim; k++){
					d[k] = p[i][k];
				}
				y[i] =  ucont(d,  m, hs, ah,  zi, rho,  Omega , w,  ar, ph, own, qi, na,  pr_ucont);
			}

			amoeba_nr(p, y, ndim, ftol, ucont, &nfunc, m , hs, zi, rho, Omega, w, ar, ph, ah, own, qi,  na, pr_amoeba); 
			
			if (counter ==0 && y[0] >= ytemp) {
				break;
			}
			if (y[0] < ytemp){
				ytemp = y[0]; 
				hpprime_temp = p[0][0]; 
				haprime_temp = p[0][1];
				aprime_temp = p[0][2];
				xprime_temp = 0.0;

				new_vals = 1; 
			}
			hpprime1 = hpprime_temp; 
			haprime1 = haprime_temp; 
			aprime1 = aprime_temp; 

			if (counter < 3){
				hpprime2 = hpprime1 + hpsign[counter]*.1*hpprime1 + 0.1;  
				haprime2 = haprime1 + hasign[counter]*.1*haprime1+ 0.005; 
				aprime2 = aprime1 + asign[counter]*.1*aprime1 + 0.1;  
			}

			max_loops +=1; 
		}
		if (new_vals ==1){
			*ystar = ytemp; 
			*hpprime = hpprime_temp; 
			*haprime = haprime_temp; 
			*aprime  = aprime_temp;
			*xprime = xprime_temp; 
			*hsrprime= 0.0; 
		}

	} else  { //solve renters problem

		for (counter =0; counter < 1; counter++){
			
			if (aprime1 <0.0) aprime1 = 0.01; 
			if(aprime2 <0.0) aprime2 =0.015; 
			if(m- aprime1 <= 1.0e-10) aprime1 =0.9*m; 
			if (m - aprime2 <= 1.0e-10) aprime2 = 0.8*m; 

			p[0][0] = hsrprime1; 
			p[0][1] = aprime1; 
			p[1][0] = hsrprime2; 
			p[1][1] = aprime1; 
			p[2][0] = hsrprime1; 
			p[2][1] = aprime2; 

			for(i = 0; i <ndim+1; i ++){
				for (k = 0; k < ndim; k++){
					d[k] = p[i][k];
				}
				y[i] = ucont(d,  m, hs, ah,  zi, rho,  Omega , w, ar, ph,  own, qi, na,  pr_ucont); 
			}


			amoeba_r(p, y, ndim, ftol, ucont, &nfunc, m , hs, zi, rho, Omega, w, ar, ph, ah, own,qi,  na,  pr_amoeba); 
			
			if (y[0] < ytemp){
				ytemp = y[0]; 
				hsrprime_temp = p[0][0]; 
				aprime_temp = p[0][1];
				if(fabs(ytemp - y[0]) > 1.0e-10) counter =0;
				new_vals = 1; 
			}
			hsrprime1 = hsrprime_temp; 
			aprime1 = aprime_temp; 

			if (counter < 3){
				hsrprime2 = hsrprime1 + hpsign[counter]*.1*hsrprime1 + 0.1;  
				aprime2 = aprime1 + asign[counter]*.1*aprime1 + 0.1;  
			}
		}


		if (new_vals ==1){
			*ystar = ytemp; 
			*hsrprime = hsrprime_temp; 
			*aprime  = aprime_temp; 
			*xprime = 0.0;
			*haprime =0.0; 
			*hpprime = 0.0; 
		}




	}
}

