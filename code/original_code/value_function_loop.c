/*
Program: value_function_loop.c
Purpose: Iterates value function   
Calls: nelder.c nelder_nr.c
*/

#include "header.h"

extern int nm, nhs, nz, nq, *qind, *zind, *hsind, *mind, iter; 
extern double  psi_h, sigma, zeta, tol, *v0, *v1, *hppol, *hapol, *apol, *xpol,  *hsgrid, *mgrid, *ahpol,  
       *hsrpol, hpmin, Ah, *hppol_own, *hapol_own, *apol_own, **ya, **y2a, *xpol_own, theta,mdist;
extern char *folder;

double value_function_loop( double rho, double Omega, double eta, double w, double ar, double ph,  int na, int initial, int region)
{

	//Define local variables
	int hsi, mi,  zi,q,  n,  j, i,  k, jprob,  ind2, pr, own, ind3, howard, qi, hs_start, hsind1, mind1;
	double dist,  vdist, cont, *temp2, max, hs, hsprime, m, c, ah, 
	       ytemp,  hpprime1, hpprime2, aprime1, aprime2, haprime1, haprime2, xprime1, xprime2, scale, c2,
	       ystar, prime1[4], prime2[4],  scale_vec[3], scale2_vec[3], c1, hpprime, haprime, aprime, xprime, ubxprime1, 
	       hsrprime1, hsrprime2, hsrprime, ubxprime2, gprime1, gprime2, ytemp2, val_prob, diff1, choice[4], mdist1; 
	FILE *rs;
	char filename[100];

	n = nz*nhs*nm*nq; 
	val_prob =0; 	
	howard =20; 

	//initial =1 => initialize the value function
	//Initial guess for value function
	for (j = 0; j < n; j++){
		hsi = hsind[j];
		mi  = mind[j];
		hs = hsgrid[hsi];
		m = mgrid[mi];
		zi = zind[j]; 
		c =0.5*m; 
		if (initial ==1 || (initial !=1 && v0[j] < -1.0e4 )){ 
			v0[j] = pow((pow(c, zeta)*pow(hs, 1.0-zeta)), 1.0-sigma)/(1.0-sigma);
		}

		v1[j] = v0[j];
		hppol[j] = 0.0;
		hapol[j] = 0.0; 
		apol[j] = 0.0; 
		xpol[j] =0.0; 
		hsrpol[j]  = 0.0; 
	}

	dist = 1.0;
	mdist = tol * 10;
	iter = -1;

	c1 =5.0; 
	scale_vec[0] = 1.0; scale_vec[1] = 0.75;  scale_vec[2] = 1.25; 
	scale2_vec[0] = 1.0; scale2_vec[1] = 0.75; scale2_vec[2] = 1.25;

#pragma omp parallel num_threads (88) default(none) private(j,i, hsi, mi, qi,  hs, \
		max, vdist, m, k, hpprime1, hpprime2, haprime1, haprime2, aprime1, aprime2,xprime1, xprime2,hsind1, mind1, \
		ytemp, ystar, c2, scale, prime1, prime2,  zi,  pr,  ind2, hpprime, haprime, aprime, xprime, hs_start, ah, \
		ubxprime1, ubxprime2, gprime1, gprime2, ytemp2, hsrprime1, hsrprime2, hsrprime, own, ind3, diff1, choice)\
	shared(mdist,tol, iter,  mind, zind, hsind, qind,  hsgrid, mgrid,  c1, w, val_prob, psi_h,  eta, y2a, xpol_own, theta, \
			v1, v0,temp2, hppol, hapol, apol, xpol, hsrpol,  n, nhs, nm, nz,  rho, Omega,  jprob, ahpol, \
			scale_vec, scale2_vec, hpmin, ph, ar, Ah, howard, hppol_own, hapol_own, apol_own, ya, nq, na, region, mdist1)	
	{	
		while (mdist > tol) {
#pragma omp barrier
#pragma omp single
			{
				iter += 1;
				if (iter %howard ==0) mdist = -1.0;
			}
			
			//Find second derivative matrix for cubic splines

# pragma omp for
			for(j =0; j<n; j++){
				qi = qind[j]; 
				zi = zind[j]; 
				hsi =hsind[j]; 
				mi = mind[j];
				hsind1 = qi*nz*nhs + zi*nhs + hsi; 
				mind1 =  mi;
				ya[hsind1][mind1] = v0[j];
			}
#pragma omp barrier
#pragma omp single 
			{
				for(j =0; j <nq*nz; j++){
					qi = ((int)j / nz) % nq;
					zi =j % nz;  
					hs_start = qi*nz*nhs + zi*nhs; 
					splie2(hsgrid, mgrid, ya, nhs, nm,hs_start, y2a);
				}

			}

#pragma omp barrier
			
				if (iter %howard !=0) {

# pragma omp for
				for (j =0; j<n ; j++){
					pr =0;
					qi = qind[j]; 
					zi  = zind[j]; 
					hsi = hsind[j];
					mi = mind[j];
					hs = hsgrid[hsi];
					m = mgrid[mi];

					hpprime = hppol[j];
					if (na ==1){
						if (region ==0) ah = 0.0151; 
						else ah = 0.0365; 
						haprime = ah*hpprime; 
					} else 
					{ 	haprime =hapol[j];
						ah =0.0; 
					}
					aprime = apol[j];
					hsrprime = hsrpol[j];
					xprime = xpol[j]; 
					
					if ((hpprime > 0.0)) {
						own = 1; 
						choice[0] = hpprime; 
						choice[1] = haprime; 
						choice[2] = aprime;
						choice[3]  = xprime; 
					} else {
						own =0; 
						choice[0] = hsrprime; 
						choice[1] = aprime; 
					}

					v1[j]= -ucont(choice, m, hs, ah, zi,  rho,  Omega, w,  ar,  ph, own, qi, na,  pr);
				
				} //end of j loop
			} else {


# pragma omp for
				for (j = 0; j < n; j++) {  // z, hstilde, m
					pr =0;
					qi = qind[j]; 
					zi  = zind[j]; 
					hsi = hsind[j];
					mi = mind[j];
					hs = hsgrid[hsi];
					m = mgrid[mi];
					ytemp = 1.0e30;

					if (region ==0) ah = 0.0151; 
					else ah = 0.0365; 
					
					if (qi ==0) {
						own =1;	
						for (i =0; i<3; i++){
							scale = scale_vec[i]; 
							ind2 = zi*nm*nhs + hsi*nm + mi;

							if (iter ==0){
								hpprime1 = scale*2.0; 
								aprime1 = scale*10.0; 
								haprime1 = scale*0.05; 
								xprime1 = 0.57*Omega*
								exp(-theta*pow(haprime1/hpprime1, 1.0/theta))*hpprime1; 
							}else{
								hpprime1 = scale*hppol_own[ind2]; 
								aprime1 = scale*apol_own[ind2]; 
								xprime1 =scale*xpol_own[ind2]; 
								haprime1 = scale*hapol_own[ind2]; 
							}
							
							if (rho > 0.99) xprime1 = 0.0; //no risk 

							if (hpprime1 < hpmin) hpprime1 = scale/0.9*hpmin;

							prime1[0] = hpprime1; 
							prime1[1] = haprime1; 
							prime1[2] = aprime1;
							prime1[3] = xprime1; 

							for (k =0; k< 3;  k++){
								scale = scale_vec[k];

								hpprime2 =scale*hpprime1 + 0.5; 
								haprime2 = scale*haprime1 + 0.005; 
								aprime2 = scale*aprime1+ 1.0; 
								xprime2 = 0.55*Omega*
									exp(-theta*pow(haprime2/hpprime2, 1.0/theta))*hpprime2; 

								if (hpprime2< hpmin) hpprime2 = scale/0.9*hpmin+ 0.1;

								//if (rho > 0.99) xprime2 = 0.000001; 
								
								prime2[0] = hpprime2; 
								prime2[1] = haprime2; 
								prime2[2] = aprime2;
								prime2[3] = xprime2; 

								if (rho < 0.99){
								nelder(prime1, prime2, ytemp, 4,  &hpprime, &haprime, &aprime,
										&xprime, &hsrprime, &ystar, m, hs, zi, rho,
										Omega, w, ar, ph, ah, own, qi, na,  pr); 
								}
								else {
								nelder_nr(prime1, prime2, ytemp, 3,  &hpprime, &haprime, &aprime,
										&xprime, &hsrprime, &ystar, m, hs, zi, rho,
										Omega, w, ar, ph, ah, own, qi, na,  pr); 
								}
								ytemp = ystar;


							}
						}
						if (iter > 2){
							hppol_own[j] = (hpprime + hppol_own[j])/2.0; 
							hapol_own[j] = (haprime + hapol_own[j])/2.0; 
							apol_own[j] = (aprime + apol_own[j])/2.0;
							xpol_own[j] = (xprime + xpol_own[j])/2.0;
						}else {
							hppol_own[j] = hpprime; 
							hapol_own[j] = haprime; 
							apol_own[j] = aprime;
							xpol_own[j] = xprime; 
						}
					} //end of if (qi==0)
					
					own =0; 
					for (i =0; i<3; i++){
						scale = scale2_vec[i]; 

						if (qi ==0){
							if (na ==0) aprime1 = scale*(aprime +hpprime + haprime); 
							else aprime1 = scale*(aprime + (1.0+ ah)*hpprime); 
							hsrprime1 = scale*Ah*hpprime;
						}else{
							hsrprime1 = scale*0.2655; 
							aprime1 = scale*15.0; 
						}
						prime1[0] = hsrprime1; 
						prime1[1] = aprime1; 

						for (k =0; k< 3;  k++){
							scale = scale_vec[k];

							hsrprime2 =scale*hsrprime1 + 0.1; 
							aprime2 = scale*aprime1+ 1.0; 

							prime2[0] = hsrprime2; 
							prime2[1] = aprime2; 

							nelder(prime1, prime2, ytemp, 2,  &hpprime, &haprime, &aprime, 
									&xprime, &hsrprime, &ystar, m, hs, zi, rho, 
									Omega, w, ar, ph, ah,  own, qi, na,  pr); 
							ytemp = ystar; 

						}
					}
					v1[j] = -ystar; 
					hppol[j] =hpprime;
					if (na ==1) hapol[j] = ah*hpprime; 
					else hapol[j] = haprime; 
					apol[j] = aprime;
					if (rho > 0.99) xpol[j] =0.0; 
					else xpol[j] = xprime;
					hsrpol[j] = hsrprime;

					if (v1[j] >-0.0001) 
						printf("1) %f %f %f %f %f %f %i\n", v1[j], hppol[j], hapol[j], apol[j], xpol[j], hsrpol[j], j); 


					vdist = fabs(v1[j] - v0[j])/fabs(v0[j]);
#pragma omp critical 
					{
						if (vdist > mdist)
						{
							mdist = vdist;
							jprob = j;
						}

					}
				}//end of j loop

			} //end of else clause (for howard)

#pragma omp barrier
#pragma omp single
			{
				if(iter %howard==0) {
					
			printf(" jprob, qi, mi, hsi,  zi, hppol, apol, hapol, xpol, hsrpol, v1, v0: %i %i %i %i %i %f %f %f %f %f %f %f\n",
						jprob,	qind[jprob], mind[jprob], hsind[jprob], 
							zind[jprob], hppol[jprob], apol[jprob], hapol[jprob], xpol[jprob],
							hsrpol[jprob], v1[jprob], v0[jprob]); 
					printf("iter = %i, mdist = %4.8f\n", iter, mdist);
				}

				temp2 = v1;
				v1 = v0;
				v0 = temp2;
				if (iter/howard > 30){
					val_prob = mdist;
					mdist1 =mdist; 
					mdist =0;
				}
			}


		} //end of while loop 

	} // end of parallel region		
	if (iter/howard > 30) mdist = mdist1; 

	return val_prob;

} 

 /*
 int id = omp_get_thread_num();
 int data = id;
 int total = omp_get_num_threads();
 
 printf("Greetings from process %d out of %d with Data %d\n", id, total, data);
 */
