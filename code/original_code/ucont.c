/*
Program: ucont.c
Purpose: Computes u + beta*E(V), interpolates using the cubic splines
Calls: splin2.c
*/
#include "header.h"

extern double r,beta, zeta, lambda, sigma, eta, tau, theta,  *hsgrid, *mgrid, zgrid[10],
	psi_h,  psi_r, deltah, deltaa,  *v0, tmat[10][10], phi, Ah,  **ya, **y2a;
extern int nm, nhs, nz, ne, nq;

double ucont(double choice[],  double m, double hs, double ah,  int zi, double rho, double  Omega, double  w, double ar, 
		double ph, int own,  int qi, int na, int pr) 
{

	int i, j,  hi_l, hi_h, mi_l, mi_h, hi, mi, hil_mat[ne], hih_mat[ne],
	z2i,e2i,  ind, s1, s2, hs_start, oobh[2], oobm, hi_o, mi_o, ind_m, ind_hs;
	double wh, e2, wm1,  wm,  wh_mat[ne],  u, c,  hp2, ha2, a2, x2,  hsr2, minc, slope_m, slope_hs,   
	       cont, hpt2, m2,  F2, dhp2, hs2,   hsrt2, dhrp2, xrec2, px, eprob[2], temp, hsinc,  y;

	px = rho*(1.0+ lambda);
	eprob[0] = 1.0-rho; eprob[1] = rho; 

	//Indicies for transition probability
	if (zi <nz/2) {
		s1 = 0; s2 = nz/2;
	}else {
		s1 = nz/2; s2 = nz;
	}
	if (own== 1) {
		hp2 = choice[0];
		if (na ==1) ha2 = ah*hp2;
		else ha2 = choice[1]; 
		a2 = choice[2];
		if (rho > 0.99)	x2 =0.0; 
		else x2 = choice[3]; 
		c = m - a2  -ha2/(1.0+eta) - hp2;
		F2 = exp(-theta*pow(ha2/hp2, 1.0/theta));
	} else {
		hsr2 = choice[0]; 
		a2 = choice[1]; 
		x2 =0.0; 
		ha2 =0.0;
		hp2 =0.0;
		c = m - a2;
		F2 = exp(-theta*pow(ar, 1.0/theta));
	}

	if (c<0 || ((ha2/(1.0+eta) + hp2) + a2 < 0.0)) return 1.0e8 + 1000.0*fabs(c);
	else{
		u =pow((pow(c, zeta)*pow(hs, 1.0-zeta)), 1.0-sigma)/(1.0-sigma);
		//hs prime
		for (e2i=0; e2i<ne; e2i++){
			e2 = (double) e2i; 
			if (own ==1){
				dhp2 = e2*Omega*F2*hp2; 
				hpt2 = hp2 -dhp2;
				hs2 = Ah*hpt2; 
			}
			else {
				dhrp2 =  e2*Omega*F2*hsr2/Ah; 
				hsrt2 =Ah*(hsr2/Ah - dhrp2);
				if (qi ==0) hs2 = phi*hsrt2;
				else hs2 = hsrt2;
			}

			if (hs2 > hsgrid[nhs-1]){
				oobh[e2i] =2;   //too high
				hi_o = nhs-1; 
			}
			else if (hs2 < hsgrid[0]){
				oobh[e2i] =1; //too low
				hi_o = 0; 
			}else oobh[e2i] =0; //in bounds
		}	
		cont =0.0; 
		temp = 0.0; 
		for (z2i = s1; z2i< s2; z2i++){
			for (e2i =0; e2i< ne; e2i++){
				e2 = (double) e2i; 

				//Find mprime
				if (own==1){
					dhp2 = e2*Omega*F2*hp2; 
					hpt2 = hp2 -dhp2;
					xrec2 = x2 > (1.0-psi_h)*dhp2 ? (1.0 - psi_h)*dhp2 : x2;
					m2 = (1.0 - tau)* w*zgrid[z2i] + (1.0 +r)*a2 + 
						(1.0 -deltah)*hpt2 +(1.0-deltaa)*ha2/(1.0 + eta) +psi_h*dhp2+ e2*xrec2 - px*x2;
				} else {
					dhrp2 = e2*Omega*F2*hsr2/Ah; 
					hsrt2 = Ah*(hsr2/Ah - dhrp2); 
					m2 =  (1.0 - tau)*w*zgrid[z2i] + (1.0 +r)*a2 -ph*hsrt2 + psi_r*dhrp2;
				}

				if (m2 > mgrid[nm-1]){
					oobm = 2; //m is too high 
					mi_o = nm-1; 
				}
				else if (m2 < mgrid[0]){
					oobm =1; // m is too low
					mi_o = 0; 
				} else oobm =0; //m is perfect

				if (oobh[e2i] >0 && oobm ==0 ){ //h is out of bounds, m is not, find spot on the mgrid
					mi_h = nm-1; 
					for (i =1; i < nm;  i++){
						if (m2 <= mgrid[i]){
							mi_h = i;
							break;
						}
					}
					mi_l = mi_h -1;
					minc =mgrid[mi_h] - mgrid[mi_l]; 
					wm1 = (m2 - mgrid[mi_l])/minc; 
					if (wm1 < 0.0) wm1  = 0.0; 
					if (wm1 > 1.0) wm1 = 1.0;
				} 
				if (oobh[e2i] ==0 && oobm >0){ //m is out of bounds, but h is not, find spot on hgrid
					hi_h = nhs-1; 
					for (i =1; i < nhs;  i++){
						if (hs2 <=hsgrid[i]){
							hi_h = i;
							break;
						}
					}
					hi_l = hi_h -1;
					hsinc = hsgrid[hi_h] - hsgrid[hi_l];
					wh = (hs2 -hsgrid[hi_l])/hsinc; 
					hil_mat[e2i] = hi_l; 
					hih_mat[e2i] = hi_h; 
					wh_mat[e2i] = wh; 
					if (wh < 0.0) wh  = 0.0; 
					if (wh > 1.0) wh = 1.0;

				}
				//Interpolate value function
				if (oobh[e2i]>0 && oobm >0) { //m and h are out of bounds
					ind = qi*nhs*nm*nz + z2i*nhs*nm + hi_o*nm + mi_o; 
					if (oobm ==1) { //too low
						ind_m = qi*nhs*nm*nz+ z2i*nhs*nm + hi_o*nm + mi_o +1;
						slope_m = (v0[ind] - v0[ind_m])/(mgrid[mi_o] - mgrid[mi_o +1]);
					}
					else {
						ind_m = qi*nhs*nm*nz+ z2i*nhs*nm + hi_o*nm + mi_o -1;
						slope_m = (v0[ind] - v0[ind_m])/(mgrid[mi_o] - mgrid[mi_o -1]);
					}
					if (oobh[e2i] ==1) {
						ind_hs = qi*nhs*nm*nz+ z2i*nhs*nm + (hi_o+1)*nm + mi_o;
						slope_hs = (v0[ind] - v0[ind_hs])/(hsgrid[hi_o] - hsgrid[hi_o +1]);
					}
					else { 
						ind_hs = qi*nhs*nm*nz+ z2i*nhs*nm + (hi_o-1)*nm + mi_o;
						slope_hs = (v0[ind] - v0[ind_hs])/(hsgrid[hi_o] - hsgrid[hi_o -1]);
						//		printf("%i %i %i %f %f\n", hi_o, hi_o - 1, mi_o, v0[ind], v0[ind_hs]); 
					}
					if ((v0[ind] + slope_m*(m2 - mgrid[mi_o]) + slope_hs*(hs2 - hsgrid[hi_o])) < 0) {
						cont +=  tmat[zi][z2i]*eprob[e2i]*(v0[ind] + slope_m*(m2 - mgrid[mi_o]) 
							+ slope_hs*(hs2 - hsgrid[hi_o]));
					} else cont +=0.0; 
				} 
				else if (oobh[e2i] >0 && oobm ==0) { //hs is out of bounds, m is in bounds
					for (j = 0; j < 2; j++){

						wm = (j==0) ? 1.0 - wm1 : wm1; 
						mi = (j==0) ? mi_l : mi_h;
						ind = qi*nhs*nm*nz+z2i*nhs*nm + hi_o*nm + mi; 

						if (oobh[e2i] ==1) {
							ind_hs = qi*nhs*nm*nz + z2i*nhs*nm + (hi_o+1)*nm + mi;
							slope_hs = (v0[ind] - v0[ind_hs])/(hsgrid[hi_o] - hsgrid[hi_o +1]);
						}
						else { 
							ind_hs = qi*nhs*nm*nz+z2i*nhs*nm + (hi_o-1)*nm + mi;
							slope_hs = (v0[ind] - v0[ind_hs])/(hsgrid[hi_o] - hsgrid[hi_o -1]);
						}

						if (v0[ind] + slope_hs*(hs2 - hsgrid[hi_o]) < 0)
							cont+=tmat[zi][z2i]*eprob[e2i]*wm*(v0[ind] + slope_hs*(hs2 - hsgrid[hi_o]));
						else  cont+=0.0;
					}
				}	
				else if (oobh[e2i] ==0 && oobm >0) { //hs is in  bounds, m is out of bounds
					for (j = 0; j < 2; j++){
						wh = (j ==0) ? 1.0 - wh_mat[e2i] : wh_mat[e2i]; 
						hi = (j ==0)? hil_mat[e2i]: hih_mat[e2i];
						ind = qi*nhs*nm*nz + z2i*nhs*nm + hi*nm + mi_o; 

						if (oobm ==1) {
							ind_m = qi*nhs*nm*nz + z2i*nhs*nm + hi*nm + mi_o+1;
							slope_m = (v0[ind] - v0[ind_m])/(mgrid[mi_o] - mgrid[mi_o +1]);
						}
						else { 
							ind_m = qi*nhs*nm*nz+z2i*nhs*nm + hi*nm + mi_o-1;
							slope_m = (v0[ind] - v0[ind_m])/(mgrid[mi_o] - mgrid[mi_o -1]);
						}

						if ((v0[ind] + slope_m*(m2 - mgrid[mi_o])) < 0.0 ) {
							cont+=tmat[zi][z2i]*eprob[e2i]*wh*(v0[ind] + slope_m*(m2 - mgrid[mi_o]));
						} else  cont +=0.0;  

					}
				}	
				else { //in bounds: use cubic spline interpolation

					hs_start = qi*nz*nhs + z2i*nhs;
					splin2(hsgrid, mgrid, ya, y2a, nhs, nm, hs_start, hs2, m2, &y);
					if (y >0) {
						printf("IN UCONT: %i %f %f %f %i %f\n", hs_start, hs2, m2, ah,  na,  y);  
						y =0; 
					}
					cont += tmat[zi][z2i]*eprob[e2i]*y; 
					temp += tmat[zi][z2i]*eprob[e2i]; 
				}

			}//end of e2i loop
		}//end of z2i loop
		
		if (hp2 < 0.00001 && a2 < 0.00001 && own ==1) printf("%f\n", cont); 


	//	if (cont >0) printf("1)%i %i %i %f %f %f %f %f %f \n", oobm, oobh[0], oobh[1], hs2, slope_hs, m, hs, c, cont);
	//	if (u >0) printf("2)%f %f %f\n", c, hs, u);

		return -(u + beta*cont);
		

	}//End of else clause
}
