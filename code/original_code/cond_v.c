/*
Program: cond_v.c
Purpose: Computes the contiution value conditional on the realization of the storm shock
*/


#include "header.h"

extern double r,beta, zeta, lambda, sigma, eta, tau, theta,  *hsgrid, *mgrid, zgrid[10],
	psi_h,  psi_r, deltah, deltaa,  *v0, tmat[10][10], phi, Ah,  **ya, **y2a,  *hppol, *hapol, *apol, *xpol, *hsrpol,
	*cont_e0, *cont_e1, *f1; 
extern int  nhs, nm, nz, ne, nq, *qind,  *zind, *hsind,  *mind; 


void cond_v(double rho, double  Omega, double  w, double ar, double ph) 
{

	int i, j,  hi_l, hi_h, mi_l, mi_h, hi, mi,  hil_mat, hih_mat, qi, zi,n, 
	    z2i,e2i,  ind, s1, s2, hs_start, oobh, oobm, hi_o, mi_o, ind_m, ind_hs;
	double wh, e2, wm1,  wm,  wh_mat, u, c,  hp2, ha2, a2, x2, hrp2,  hsr2, minc, slope_m, slope_hs,   
	       cont, hpt2, m2,  Fh2,  Fhr2, dhp2, hs2,   hsrt2, dhrp2, xrec2, px, temp, hsinc, y ;


	n = nm*nz*nhs*nq;
	px = rho*(1.0 + lambda); 

	//solve for conditional value function
	for (e2i = 0; e2i < 2; e2i ++) {
		e2 = (double) e2i; 
		for (j =0; j < n; j++) {
			if (f1[j] > 1.0e-8){
			//	printf("%i %i\n", j, n); 
			qi = qind[j]; 
			zi = zind[j]; 
			hp2 = hppol[j]; 
			ha2 = hapol[j]; 
			a2 = apol[j]; 
			x2 = xpol[j]; 
			hsr2 = hsrpol[j]; 
			hrp2 = hsr2/Ah;

			//Indicies for transition probability
			if (zi <nz/2) {
				s1 = 0; s2 = nz/2;
			}else {
				s1 = nz/2; s2 = nz;
			}

			if (hp2 >0.0) Fh2 = exp(-theta*pow(ha2/hp2, 1.0/theta));
			else Fh2 =0.0;
			Fhr2 = exp(-theta*pow(ar, 1.0/theta)); 

			//hs prime
			if (hp2 >0){
				dhp2 = e2*Omega*Fh2*hp2; 
				hpt2 = hp2 -dhp2;
				hs2 = Ah*hpt2; 
			}
			else {
				dhrp2 =  e2*Omega*Fhr2*hsr2/Ah; 
				hsrt2 =Ah*(hsr2/Ah - dhrp2);
				if (qi ==0) hs2 = phi*hsrt2;
				else hs2 = hsrt2;
			}

			if (hs2 > hsgrid[nhs-1]){
				oobh =2;   //too high
				hi_o = nhs-1; 
			}
			else if (hs2 < hsgrid[0]){
				oobh =1; //too low
				hi_o = 0; 
			}else oobh =0; //in bounds

			cont =0.0; 
			temp = 0.0; 
			for (z2i = s1; z2i< s2; z2i++){
				//Find mprime
				if (hp2 >0){
					dhp2 = e2*Omega*Fh2*hp2; 
					hpt2 = hp2 -dhp2;
					xrec2 = x2 > (1.0-psi_h)*dhp2 ? (1.0 - psi_h)*dhp2 : x2;
					m2 = (1.0 - tau)* w*zgrid[z2i] + (1.0 +r)*a2 + 
						(1.0 -deltah)*hpt2 +(1.0-deltaa)*ha2/(1.0 + eta) +psi_h*dhp2+ e2*xrec2 - px*x2;
				} else {
					dhrp2 = e2*Omega*Fhr2*hsr2/Ah; 
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

				if (oobh >0 && oobm ==0 ){ //h is out of bounds, m is not, find spot on the mgrid
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
				if (oobh==0 && oobm >0){ //m is out of bounds, but h is not, find spot on hgrid
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
					hil_mat= hi_l; 
					hih_mat = hi_h; 
					wh_mat = wh; 
					if (wh < 0.0) wh  = 0.0; 
					if (wh > 1.0) wh = 1.0;

				}
				//Interpolate value function
				if (oobh>0 && oobm >0) { //m and h are out of bounds
					ind = qi*nhs*nm*nz + z2i*nhs*nm + hi_o*nm + mi_o; 
					if (oobm ==1) { //too low
						ind_m = qi*nhs*nm*nz+ z2i*nhs*nm + hi_o*nm + mi_o +1;
						slope_m = (v0[ind] - v0[ind_m])/(mgrid[mi_o] - mgrid[mi_o +1]);
					}
					else {
						ind_m = qi*nhs*nm*nz+ z2i*nhs*nm + hi_o*nm + mi_o -1;
						slope_m = (v0[ind] - v0[ind_m])/(mgrid[mi_o] - mgrid[mi_o -1]);
					}
					if (oobh==1) {
						ind_hs = qi*nhs*nm*nz+ z2i*nhs*nm + (hi_o+1)*nm + mi_o;
						slope_hs = (v0[ind] - v0[ind_hs])/(hsgrid[hi_o] - hsgrid[hi_o +1]);
					}
					else { 
						ind_hs = qi*nhs*nm*nz+ z2i*nhs*nm + (hi_o-1)*nm + mi_o;
						slope_hs = (v0[ind] - v0[ind_hs])/(hsgrid[hi_o] - hsgrid[hi_o -1]);
						//		printf("%i %i %i %f %f\n", hi_o, hi_o - 1, mi_o, v0[ind], v0[ind_hs]); 
					}
					if ((v0[ind] + slope_m*(m2 - mgrid[mi_o]) + slope_hs*(hs2 - hsgrid[hi_o])) < 0) {
						cont +=  tmat[zi][z2i]*(v0[ind] + slope_m*(m2 - mgrid[mi_o]) 
								+ slope_hs*(hs2 - hsgrid[hi_o]));
					} else cont +=0.0; 
				} 
				else if (oobh >0 && oobm ==0) { //hs is out of bounds, m is in bounds
					for (j = 0; j < 2; j++){

						wm = (j==0) ? 1.0 - wm1 : wm1; 
						mi = (j==0) ? mi_l : mi_h;
						ind = qi*nhs*nm*nz+z2i*nhs*nm + hi_o*nm + mi; 

						if (oobh ==1) {
							ind_hs = qi*nhs*nm*nz + z2i*nhs*nm + (hi_o+1)*nm + mi;
							slope_hs = (v0[ind] - v0[ind_hs])/(hsgrid[hi_o] - hsgrid[hi_o +1]);
						}
						else { 
							ind_hs = qi*nhs*nm*nz+z2i*nhs*nm + (hi_o-1)*nm + mi;
							slope_hs = (v0[ind] - v0[ind_hs])/(hsgrid[hi_o] - hsgrid[hi_o -1]);
						}

						if (v0[ind] + slope_hs*(hs2 - hsgrid[hi_o]) < 0)
							cont+=tmat[zi][z2i]*wm*(v0[ind] + slope_hs*(hs2 - hsgrid[hi_o]));
						else  cont+=0.0;
					}
				}	
				else if (oobh ==0 && oobm >0) { //hs is in  bounds, m is out of bounds
					for (j = 0; j < 2; j++){
						wh = (j ==0) ? 1.0 - wh_mat : wh_mat; 
						hi = (j ==0)? hil_mat: hih_mat;
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
							cont+=tmat[zi][z2i]*wh*(v0[ind] + slope_m*(m2 - mgrid[mi_o]));
						} else  cont +=0.0;  

					}
				}	
				else { //in bounds: use cubic spline interpolation

					hs_start = qi*nz*nhs + z2i*nhs;
					splin2(hsgrid, mgrid, ya, y2a, nhs, nm, hs_start, hs2, m2, &y);
					if (y >0) {
						y =0; 
					}
					cont += tmat[zi][z2i]*y; 
				}

			}//end of z2i loop
			if (e2i ==0) cont_e0[j]  = cont; // in state space point j, what is the continuation value if not hit by a storm
			else cont_e1[j] = cont; //in state space point j, what is the continuation value if hit by a storm
			}
			else{
				cont_e0[j]=0; 
				cont_e1[j] =0; 
			}
			}//end of j loop
	}//end of e2i loop

}
