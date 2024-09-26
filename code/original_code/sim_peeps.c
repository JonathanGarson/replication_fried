/*
Program: sim_peep.c
Purpose: uses approximate aggregation to construct the distribution of households across the state space. 
*/


#include "header.h"

extern double *v0, *f1, *f2,  *hsgrid, *mgrid,  zgrid[10], *hppol, *hapol, *apol, *xpol, *hsrpol, *wealth, *earns, *finc, 
       tmat[10][10], dtol, r, theta, psi_h, psi_r, deltah,deltaa,  lambda, tau, phi, Ah, *inc, chi, mdist_sim_peeps; 
extern int  nhs, nm, nz, ne, nq, *qind,  *zind, *hsind,  *mind, *ooh, iter_sim_peeps; 

void sim_peeps(double rho, double Omega, double eta, double  w, double ar, double ph, 
		double *hp_rag, double *ha_rag, double *a_rag, 
		double *dhp_rag, double *m_rag, double *x_dhp_rag, double *hs_rag, double *hpt_rag, double *hsr_rag, double *frac_own)
{

	int i, j,k,mi, hi,ind, e2i, z2i, n, hi_l, hi_h, mi_l, mi_h, mih_mat[nz][ne], mil_mat[nz][ne], hil_mat[ne], 
	    hih_mat[ne], iter, zi, qi, s1, s2, counter;

	n = nm*nz*nhs*nq;
	
	double hp2, ha2, a2, x2, hpt2, hsr2, hsrt2, hs2, Fh2,minc, e2,
	       dhp2, Fhr2, dhrp2,  wh, wm, wh_mat[ne], wm_mat[nz][ne],  dist, mdist, *temp2, frac_own_temp,
		m, hpt,  xrec2, m2, hs,  px,  hp, ha, hsr, a, z,x,  tot, eprob[2],dhp, hrp2, hsinc, temp, x_dhp_temp;


	px = rho*(1.0 + lambda); 
	eprob[0] = 1.0-rho; eprob[1] = rho; 	
	//initialize distributionhs 
	for (j =0; j<n; j++) {
		if (qind[j]==0)  f1[j] = 1.0/(n/(2.0*chi)); //discount renting
		 else f1[j] = 1.0/(n/(2.0*(1.0 - chi))); //Don't discount renting (always rent)
	       // f1[j] = 1.0/n;  
		f2[j] =0.0; 
	}

	//solve for invariant distribution
	iter = -1;
	mdist = 1.0;
	while (mdist > dtol) {
		iter +=1;
		tot =0.; 
		for (i = 0; i < n; i++){
			 tot+=f1[i];
		}
	//	if (iter %1 ==0) printf("iter total n: %i %f %i\n",iter,  tot, n); 

		for (j =0; j < n; j++) {
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

			//Find mprime and grid bounds 
			for (z2i =s1; z2i< s2; z2i++){
				for (e2i =0; e2i< ne; e2i++){
					e2 = (double) e2i;
					dhp2 = e2*Omega*Fh2*hp2; 
					hpt2 = hp2 -dhp2; 
					dhrp2 = e2*Omega*Fhr2*hrp2; 
					hsrt2 = Ah*(hrp2 - dhrp2); 
					
					//xrec2 = x2 > (1.0-deltah)*dhp2 - psi_h*dhp2 ? (1.0-deltah)*dhp2 - psi_h*dhp2  : x2;
					xrec2 = x2 > (1.0-psi_h)*dhp2 ? (1.0 - psi_h)*dhp2 : x2;

					m2 = (1.0 - tau)* w*zgrid[z2i] + (1.0 +r)*a2 + (1.0 -deltah)*hpt2
						+(1.0-deltaa)*ha2/(1.0 + eta) + psi_h*dhp2+ e2*xrec2 - px*x2 - ph*hsrt2 
						+ psi_r*dhrp2; 

					mi_h = nm-1; 
					for (i =1; i < nm;  i++){
						if (m2 <= mgrid[i]){
							mi_h = i;
							break;
						}
					}
					mi_l = mi_h -1;
					minc = mgrid[mi_h] - mgrid[mi_l];
					wm = (m2 - mgrid[mi_l])/minc; 
					if (wm < 0.0) wm  = 0.0; 
					if (wm > 1.0) wm = 1.0;
					wm_mat[z2i][e2i] = wm; //weight on mi_h; 
					mih_mat[z2i][e2i] = mi_h; 
					mil_mat[z2i][e2i] = mi_l;
				}
			}

			//Find hsprime and grid bounds
			for (e2i=0; e2i<ne; e2i++){
				e2 = (double) e2i; 
				dhp2 =e2*Omega*Fh2*hp2; 
				hpt2 = hp2 -dhp2; 

				dhrp2 =  e2*Omega*Fhr2*hrp2; 
				hsrt2 = Ah*(hrp2 - dhrp2);

				if (qi ==0) hs2 = Ah*hpt2 + phi*hsrt2; 
				else hs2 = Ah*hpt2 + hsrt2; 
				
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
				if (wh < 0.0) wh  = 0.0; 
				if (wh > 1.0) wh = 1.0;
				wh_mat[e2i] = wh; 
				hil_mat[e2i] = hi_l; 
				hih_mat[e2i] = hi_h;
			}

			


			//Find new distribution shares
			for (e2i =0; e2i< ne; e2i++){
				for(i =0; i<2; i++){
					wh = (i ==0) ? 1.0 - wh_mat[e2i] : wh_mat[e2i]; 
					hi = (i ==0)? hil_mat[e2i]: hih_mat[e2i]; 
					for (z2i = s1; z2i < s2; z2i++){
						for(k =0; k < 2; k++){
							wm = (k==0) ? 1.0 - wm_mat[z2i][e2i]: wm_mat[z2i][e2i]; 
							mi = (k==0) ? mil_mat[z2i][e2i] : mih_mat[z2i][e2i]; 
							ind = qi*nz*nhs*nm + z2i*nhs*nm + hi*nm + mi; 
							f2[ind] += tmat[zi][z2i]*eprob[e2i]*wh*wm*f1[j];
						}

					}
				}
			}

		}

		//find max distance between distributionhs
		mdist =-1.0; 
		for (j =0; j< n; j++) {
			dist = fabs(f2[j] - f1[j]);
			if (dist > mdist) mdist = dist;
		}

		//update distributions
		temp2 = f2;
		f2 = f1;
		f1 = temp2;
		//initialize f2 to zero
		for (j = 0; j< n; j++) f2[j] =0.0; 
		if (iter > 1500) mdist =0.0; 
	//	printf("sim peeps: iter  = %i mdist = %12.10f\n", iter, mdist); 
	}

	//Calculate aggregates

	*a_rag = 0.0; *ha_rag =0.0; *hp_rag = 0.0; *dhp_rag =0.0; *m_rag = 0.0; x_dhp_temp =0.0;  *hpt_rag =0.0; *hsr_rag =0.0; 
	*frac_own =0.0; *hs_rag =0.0; frac_own_temp =0.0;

	for (j =0; j<n; j++){
		z = zgrid[zind[j]]; 
		m = mgrid[mind[j]]; 
		hs = hsgrid[hsind[j]];
		ha = hapol[j]; 
		hp = hppol[j];
		a = apol[j]; 
		x = xpol[j];
		hsr = hsrpol[j];
		*m_rag +=m*f1[j]; 
		*hs_rag +=hs*f1[j]; 
		*hp_rag += hp*f1[j]; 
		*ha_rag += ha*f1[j]; 
		*a_rag += a*f1[j];
		*hsr_rag +=hsr*f1[j];
		if (hp >0.0){
			dhp = rho*Omega*hp*exp(-theta*pow(ha/hp, 1.0/theta));
			if (dhp <0) dhp =0; 
			*hpt_rag += (hp - dhp)*f1[j];
			*dhp_rag += dhp*f1[j];
			*frac_own += f1[j];
			frac_own_temp += f1[j];
			 ooh[j] = 1; 
			x_dhp_temp +=rho*x*f1[j]/dhp;
		//	printf("%i %f %f %f %f %f %f %f %f\n", j,  rho*x/dhp, x, dhp, rho, Omega, hp, ha, theta);
		} 
		else ooh[j] =0;
	}
	
	*x_dhp_rag = x_dhp_temp/frac_own_temp;
	printf("x_dhp_rag, x_dhp_temp, frac_own_temp: %f %f %f\n", *x_dhp_rag, x_dhp_temp, frac_own_temp); 


	//For wealth quintiles

	for (j=0; j<n; j++){ 
		wealth[j] = hapol[j] + hppol[j] + apol[j]; 
		earns[j] = zgrid[zind[j]]*w;
	//	inc[j] =  earns[j] + r*apol[j];
	}

	//Form income quintiles
	temp  =0.0; 
	for (j = 0; j < n; j ++){
		zi = zind[j];
		if (zi <nz/2) {
			s1 = 0; s2 = nz/2;
		}else {
			s1 = nz/2; s2 = nz;
		}
		counter = 0; 
		for (i =s1; i < s2; i ++){
			inc[j*nz/2+counter] = zgrid[i]*w + r*apol[j]; 
			finc[j*nz/2+counter] = f1[j]*tmat[zi][i]; 
			
			if (hapol[j]  + hppol[j] > 0.0) ooh[j*nz/2 + counter] = 1; 
			else ooh[j*nz/2 + counter] =0;
			temp += f1[j]*tmat[zi][i]; 
			
			counter +=1; 
		}
	}



printf("Fraction of homeowners: %f \n", *frac_own); 
//printf("Temp = %f\n", temp); 














}


