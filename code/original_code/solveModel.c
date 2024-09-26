/*
Program: solveModel.c
Purpose: Solves model, computes aggregrate quantities, reports targets and moments  
Calls: err_lb.c err_kakp.c err_hrahrp.c value_function_loop.c sim_peeps.c save_results.c
*/

#include "header.h"

extern double r, beta,  theta, lambda, zeta, psi_ky, psi_kr,  psi_h, psi_r,deltaa,  deltah, deltay, eta, tfp_l, tfp_h, alpha,  sigma, 
       rho_l, rho_h, tau,  *v0, *v1, *hppol, *hapol, *apol, *xpol, Omega_hh, hpmin,phi, Ah, Omegah_Omegal, chi,
        *v_h_base, *v_l_base, *v_h, *v_l, ttol,  *hapol_l_base, *hapol_h_base, Omegak_Omegah, *f1;
extern double rel_damage_target, subs_output_target, aid_output_target,  wealth_output_target,mdist,
        housing_damage_total_damage_target, rel_gdp_target, 
 	 owner_frac_target, rent_owned_target, owned_total_target, renter_fema_damage_target, 
       fema_damage_target, owner_fema_damage_target, renter_aid_total_aid_target, owner_aid_total_aid_target, 
       rel_damage_firms_target, rel_damage_renters_target, housing_capital_target, rental_housing_aid_total_aid_target,
       *earns, *wealth, *fe_pop, *fw_pop,  *earns_pop, *wealth_pop, *inc_pop, *inc, *fi_pop, mdist_sim_peeps,
       owner_renter_inc_target, *renter_inc_vec, *owner_inc_vec, *fr, *fo, *fr1, *fo1, *finc, rel_aid_target, *cont_e0, *cont_e1; 

extern int nm, nhs,  ne, nz, nq, cal_iter, use_saved, *ooh, ccexp, iter, iter_sim_peeps; 	
extern char *folder;

void solveModel(char *exper, FILE *flog)
{
	FILE *rs;
	int i, j, n,  initial, tau_iter, tau_prob, rcounter, ocounter, na;
	double time_spent,  hp_rag ,ha_rag, a_rag, m_rag, hpt_rag, kp_vec[2],
	ka_vec[2],  kptilde_vec[2], w_vec[2], dk_vec[2], ika_vec[2],
	 hp_vec[2], ha_vec[2], a_vec[2], m_vec[2], hpt_vec[2], d_agg, dhp_agg, xrec_agg,dhp_rag, 
	dh_vec[2],  aid_agg, subs_agg, gdp, d_vec[2], ho_agg, rel_damage, wealth_output, aid_output, 
	subs_output, housing_capital, tfp_vec[2], rho_vec[2], Omegah_vec[2], Omegak_vec[2], Omega_lh, Omega_hk,
	Omega_lk, k_agg, kp_k, mid, hr_agg, ihr_agg, hs_vec[2], hs_rag, ar_vec[2], ay_vec[2], Fy_vec[2], Fr_vec[2], F, Fprime,  
	 iha_rag, iha_vec[2], a_agg, rel_gdp, x_dhp_rag,  cons_p, cons_a, 
	bot, top, ph_vec[2],  hsr_rag, hrp_vec[2], dhr_vec[2],  hra_vec[2], spend, rev, ytilde_vec[2], y_vec[2], 
	tdiff, lb, ik_agg,  x_agg,  housing_damage_total_damage, x_dhp_vec[2], 
	insurance_house_damage,  val_prob, owner_frac, rent_owned, owner_frac_vec[2], owner_frac_rag, owned_total, 
	fema_damage,  rel_damage2, renter_aid_total_aid, owner_aid_total_aid, foc_hrp[2], foc_kp[2], foc_ka[2],  
	foc_hra[2],  aid_vec[2], dhp_vec[2], rel_fema_aid1, rel_fema_aid2, rel_fema_aid3,rel_fema_aid4,  gdp_vec[2], 
	e1, e2, e3, e4, e5, w1, w2, w3, w4, w5, tmp, csw, cse, w_tot, e_tot, i1, i2, i3, i4, i5, i_tot, csi, rel_aid,  
	owner_inc, renter_inc, owners, renters, owner_med, renter_med, rental_housing_aid_total_aid, temp[4], temp3[4];  
	clock_t tstart, tend; 
	char filename[100];

	int ishares =1; //set equal to 1 to correctly calculate income shares by quintile and mean owner and renter income 

	n = nz*nhs*nm*nq;
	Omega_lh = Omega_hh/Omegah_Omegal; //storm seerity for households in the low risk region 
	Omega_hk = Omegak_Omegah*Omega_hh; 
	Omega_lk = Omegak_Omegah*Omega_lh;

	if (strcmp(exper, "e")==0) {
		Omega_lh = (Omega_hh + Omega_lh)/2.0; 
		Omega_hh = Omega_lh; 
		Omega_hk =(Omega_hk +Omega_lk)/2.0; 
		Omega_lk = Omega_hk; 
	}
	printf("Omegas: %f %f %f %f\n", Omega_lh, Omega_hh, Omega_lk, Omega_hk);
	printf("rhos: %f %f\n", rho_l, rho_h); 

	//Region vectors
	tfp_vec[0] = tfp_l; tfp_vec[1] = tfp_h; 
	rho_vec[0] = rho_l; rho_vec[1] = rho_h; 
	Omegah_vec[0] = Omega_lh; Omegah_vec[1] = Omega_hh; 
	Omegak_vec[0] = Omega_lk; Omegak_vec[1] = Omega_hk; 
	
	if (strcmp(exper, "ccna") !=0){

		//Solve for fraction of adaptation capital
		for (i =0; i <2; i++) {

			lb  =  rtbis(err_lb, 1.0, 0.0,  1.0e-20, Omegak_vec[i], rho_vec[i], tfp_vec[i]);
			if (lb <0) lb =0;
			printf("Omega, rho, rho*Omega: %f %f %f\n", Omegak_vec[i], rho_vec[i], rho_vec[i]*Omegak_vec[i]); 
			ay_vec[i] = rtbis(err_kakp, 1.0, lb,  1.0e-20, Omegak_vec[i], rho_vec[i], tfp_vec[i]);

			Fy_vec[i]  = exp(-theta*pow(ay_vec[i], 1.0/theta)); 
			Fprime = -exp(-theta*pow(ay_vec[i], 1.0/theta))*pow(ay_vec[i], 1.0/theta-1.0); 
			cons_p = -ay_vec[i]; //krp*dar/dkrp, productive capital foc
			cons_a = 1.0; //adaptive capital foc

			//kp
			top = 1.0 - rho_vec[i]*Omegak_vec[i]*(Fy_vec[i] + Fprime*cons_p); 
			bot = r + deltay + rho_vec[i]*Omegak_vec[i]*(1.0 - deltay - psi_ky)*(Fy_vec[i] + Fprime*cons_p); 
			kp_vec[i] = pow(top/bot, 1.0/(1.0 - alpha))*(1.0/(1.0 - rho_vec[i]*Omegak_vec[i]*Fy_vec[i]))
				*pow(alpha*tfp_vec[i], 1.0/(1.0-alpha));

			w_vec[i] = tfp_vec[i]*(1.0-alpha)*pow(1.0-rho_vec[i]*Omegak_vec[i]*Fy_vec[i], alpha)*pow(kp_vec[i], alpha); 
			
			foc_kp[i] = alpha*pow(tfp_vec[i], 1.0/alpha)*pow((1.0-alpha)/w_vec[i], (1.0-alpha)/alpha) 
				-rho_vec[i]*Omegak_vec[i]*(alpha *pow(tfp_vec[i], 1.0/alpha)*pow((1.0-alpha)/w_vec[i], (1.0-alpha)/alpha)
						+ 1.0 - deltay - psi_ky)*(Fy_vec[i] + Fprime*cons_p) - r - deltay; 
	
			foc_ka[i] = -rho_vec[i]*Omegak_vec[i]*(1.0+eta)*
				( alpha*pow(tfp_vec[i], 1.0/alpha)*pow((1.0-alpha)/w_vec[i], (1.0-alpha)/alpha) + 1.0 - deltay - psi_ky)*
				Fprime*cons_a - r- deltaa; 

			temp3[i+2] = Omegak_vec[i]*Fy_vec[i]; 

			//Solve for fraction  of rental adaptation  housing
			ar_vec[i] =  rtbis(err_hrahrp, 1.0, 0.0, 1.0e-20, Omegah_vec[i], rho_vec[i], tfp_vec[i]);
			
			//Find ph
			Fr_vec[i]  = exp(-theta*pow(ar_vec[i], 1.0/theta)); 
			Fprime = -exp(-theta*pow(ar_vec[i], 1.0/theta))*pow(ar_vec[i], 1.0/theta-1.0); 
			cons_p = -ar_vec[i]; //krp*dar/dkrp, productive capital foc
			cons_a = 1.0; //adaptive capital foc
			
			top = r+deltah + rho_vec[i]*Omegah_vec[i]*(1.0 - deltah - psi_kr)*(Fr_vec[i] + Fprime*cons_p); 
			bot = Ah*(1.0 - rho_vec[i]*Omegah_vec[i]*(Fr_vec[i] + Fprime*cons_p)); 
			ph_vec[i]  = top/bot;
			
			foc_hra[i] = -rho_vec[i]*(1.0 +eta)*(ph_vec[i]*Ah + 1.0- deltah - psi_kr)*Omegah_vec[i]*Fprime*cons_a - r - deltaa;  
			foc_hrp[i] = ph_vec[i]*Ah - rho_vec[i]*Omegah_vec[i]*(ph_vec[i]*Ah + 1.0 - deltah - psi_kr)*
				(Fr_vec[i] + Fprime*cons_p) - r - deltah; 
			
			temp3[i] = Omegah_vec[i]*Fr_vec[i];

		}
	} else {
		//Read in baseline values of ay and ar
		sprintf(filename, "./%s/adapt_ratios.dat", folder);
		rs = fopen(filename, "r");
		if(rs == NULL) printf("Error opening file"); 
		j =0; 
		while(j < 5){
			fscanf(rs, "%lf", &temp[j]);
			j +=1; 
		}
		fclose(rs);
		ar_vec[0] = temp[0]; ar_vec[1] = temp[1];
		ay_vec[0] = temp[2]; ay_vec[1] = temp[3];

		for(i = 0; i < 2; i++){
			Fy_vec[i]  = exp(-theta*pow(ay_vec[i], 1.0/theta)); 
			Fr_vec[i]  = exp(-theta*pow(ar_vec[i], 1.0/theta));
		}

		//Find ph and Kp
		for (i =0; i < 2; i++){
			top = alpha*tfp_vec[i]*pow((1.0 - rho_vec[i]*Omegak_vec[i]*Fy_vec[i]), alpha);
			bot = ay_vec[i]/(1.0+eta)*(r + deltaa) + r + deltay + rho_vec[i]*Omegak_vec[i]*Fy_vec[i]*
				(1.0 - deltay - psi_ky); 
			kp_vec[i] = pow(top/bot, 1.0/(1.0-alpha));
			
			top = r+deltah + (r+deltaa)*(ar_vec[i]/(1.0+eta)) + rho_vec[i]*Omegah_vec[i]*Fr_vec[i]*(1.0 - deltah- psi_kr); 
			bot = Ah*(1.0 - rho_vec[i]*Omegah_vec[i]*Fr_vec[i]); 
			ph_vec[i] = top/bot;
		}

	}
	

	if (strcmp(exper, "b") ==0){
		
		//Check FOCS
		printf("Final-good Focs: %f %f %f %f\n", foc_kp[0], foc_kp[1], foc_ka[0], foc_ka[1]); 
		printf("Rental housing Focs: %f %f %f %f\n", foc_hrp[0], foc_hrp[1], foc_hra[0], foc_hra[1]); 
		printf("Rental housing damage: %f %f %f\n", temp3[0], temp3[1], temp3[1]/temp3[0]); 
		printf("final good  damage: %f %f %f\n", temp3[2], temp3[3], temp3[3]/temp3[2]); 
		
		// Save baseline values  of rental housing and final-good fractions of adaptation capital
		sprintf(filename, "./%s/adapt_ratios.dat", folder);
		rs = fopen(filename, "w"); 
		if (rs==NULL) printf("The file did not open\n");
		fprintf(rs, "%12.8f\n %12.8f\n %12.8f\n  %12.8f\n", ar_vec[0], ar_vec[1], ay_vec[0], ay_vec[1]); 
		fclose(rs); 
	
	}

	if (strcmp(exper, "c") !=0){
		printf("ar: %f %f\n", ar_vec[0], ar_vec[1]); 
		printf("ph: %f %f\n", ph_vec[0], ph_vec[1]); 
		printf("ay: %f %f\n", ay_vec[0], ay_vec[1]);
		printf("kp: %f %f\n", kp_vec[0], kp_vec[1]);
	}
	//Solve for regional aggregates
	for (i =0; i<2; i++){

		w_vec[i] = tfp_vec[i]*(1.0-alpha)*pow(1.0-rho_vec[i]*Omegak_vec[i]*Fy_vec[i], alpha)*pow(kp_vec[i], alpha); 
		ka_vec[i] = ay_vec[i]*kp_vec[i]; 
		kptilde_vec[i]= (1.0 - rho_vec[i]*Omegak_vec[i]*Fy_vec[i])*kp_vec[i];  
		ytilde_vec[i] = (1.0 - rho_vec[i])*tfp_vec[i]*pow(kp_vec[i], alpha) + rho_vec[i]*tfp_vec[i]*pow(kp_vec[i], alpha)*
			pow(1.0 - Omegak_vec[i]*Fy_vec[i], alpha); 
		y_vec[i] = tfp_vec[i]*pow(kp_vec[i], alpha); 
		dk_vec[i] = rho_vec[i]*Omegak_vec[i]*Fy_vec[i]*kp_vec[i]; 
		ika_vec[i] = deltaa*ka_vec[i]; 

		if (strcmp(exper, "c") !=0){
			printf("Region aggregates:  (KP, KA, KPtilde, w) (%4.5f, %4.5f, %4.5f, %4.5f)\n", 
					kp_vec[i], ka_vec[i], kptilde_vec[i], w_vec[i]); 

		}
	}	

	if (strcmp(exper, "c")!=0) cal_iter =1;

	time(&tstart);
	//Solve_model in each region
	tau_iter=-1;
	tdiff = 1.0; 
	tau_prob =0; 
	while (tdiff > ttol){
		tau_iter +=1;
		for (i =0; i <2;  i++) {
			if (tau_iter ==0 && cal_iter ==1) initial =1; 
			else initial =0; 

			//Use previous value as the inittal guess
			if (tau_iter > 0){
				for(j =0; j < n; j++){
					if(i ==0){
						v0[j] = v_l[j]; 
					}else{
						v0[j] = v_h[j]; 
					}
				}
			}
			
			if (use_saved==1 && tau_iter ==0){
				sprintf(filename, "./%s/v0_b%i.dat", folder, i);
				rs = fopen(filename, "r");
				if(rs == NULL) printf("Error opening file"); 
				j =0; 
				while(j < n){
				//if (i==0) fscanf(rs, "%lf", &v0[j]);
				fscanf(rs, "%lf", &v0[j]);
				j +=1; 
				}
				fclose(rs);
				
				if (i ==0){
					rs = fopen("tau_saved", "r"); 
					fscanf(rs, "%lf", &tau);
					fclose(rs);
					
				}
				initial =0;
			}
			printf("tau = %f\n", tau); 

			//Use low risk region as the initial guess for the high-risk region no-risk experiments
			if ((strcmp(exper, "ccnr") ==0 ||  strcmp(exper, "nr") ==0)  && i ==1 && tau_iter ==0){
				for (j =0; j < n; j++){
					v0[j] = v_l[j]; 
				}
				initial =0; 
			}



			if (strcmp(exper, "c") !=0){
				printf("***Solving  model in region %i tau_iter %i ***\n", i, tau_iter);
			}

			if (strcmp(exper, "ccna") !=0) na =0; //if not the no-adaptation experiment:
			else na =1; 
			
			val_prob = value_function_loop(rho_vec[i], Omegah_vec[i],  eta, w_vec[i],  ar_vec[i], ph_vec[i], na,  initial, i);
			sim_peeps(rho_vec[i], Omegah_vec[i], eta,  w_vec[i],  ar_vec[i], ph_vec[i], 
					&hp_rag, &ha_rag, &a_rag, &dhp_rag, &m_rag, &x_dhp_rag, &hs_rag, &hpt_rag, &hsr_rag, 
					&owner_frac_rag);
			
			cond_v(rho_vec[i],Omegah_vec[i], w_vec[i], ar_vec[i],  ph_vec[i]); 

			hp_vec[i] = hp_rag; 
			hpt_vec[i] = hpt_rag; 
			m_vec[i] = m_rag; 
			ha_vec[i] = ha_rag; 
			a_vec[i] = a_rag; 
			dhp_vec[i] = dhp_rag; 
			x_dhp_vec[i] = x_dhp_rag; //Insurance purchases, not claims
			hs_vec[i] = hs_rag;
			hrp_vec[i] = hsr_rag/Ah; 
			hra_vec[i] = ar_vec[i]*hrp_vec[i]; 
			dhr_vec[i] = rho_vec[i]*Omegah_vec[i]*Fr_vec[i]*hrp_vec[i]; 
			owner_frac_vec[i] = owner_frac_rag;
			dh_vec[i] = dhr_vec[i] + dhp_vec[i]; 
		
			aid_vec[i] = psi_ky*dk_vec[i] +  psi_kr*dhr_vec[i] +psi_h* (dhp_vec[i]) + psi_r*(dhr_vec[i]); 
			
			//set intitial guess for next iteration of tau loop
			for (j =0; j < n; j++){
				if(i ==0) {
					v_l[j] = v0[j];
				}
				else{
					v_h[j] = v0[j];
				}
			}
			//Fill in earnings, wealth, and income and distribution population vectors
			for (j=0; j <n; j++){
				earns_pop[i*n+j] = earns[j]; 
				wealth_pop[i*n+j] = wealth[j];
				fe_pop[i*n+j] = 0.5*f1[j]; 
				fw_pop[i*n+j] = 0.5*f1[j]; 
				if (ishares ==0) fi_pop[i*n+j] = 0.5*f1[j]; 
			}
			
			if (ishares ==1){	
				for (j=0; j < n*nz/2; j ++){
					fi_pop[i*(n*nz/2)+j] = 0.5*finc[j]; 
					inc_pop[i*(n*nz/2)+j] = inc[j];
				}
			}

			if (i ==0){ //initialize values to be zero for first region
				renter_inc =0.0; 
				owner_inc = 0.0;
				renters =0.0; 
				owners= 0.0; 
				rcounter =0; 
				ocounter =0;
			}

			//Find average and median  income for homeowners and renters
			if (ishares ==0){
				for (j=0; j < n; j++){
					if (ooh[j] ==0) {
						renter_inc += inc[j]*f1[j]*0.5; 
						renters +=0.5*f1[j];
						renter_inc_vec[rcounter] = inc[j]; 
						fr1[rcounter] = 0.5*f1[j];
						rcounter +=1; 
					}
					else {
						owner_inc += inc[j]*f1[j]*0.5;
						owners += 0.5*f1[j];  
						owner_inc_vec[ocounter] = inc[j]; 
						fo1[ocounter] = 0.5*f1[j];
						ocounter+=1;
					}
				}
			}
			else{
				for (j=0; j < n*nz/2; j++){
					if (ooh[j] ==0) {
						renter_inc += inc[j]*finc[j]*0.5; 
						renters +=0.5*finc[j];
						renter_inc_vec[rcounter] = inc[j]; 
						fr1[rcounter] = 0.5*finc[j];
						rcounter +=1; 
					}
					else {
						owner_inc += inc[j]*finc[j]*0.5;
						owners += 0.5*finc[j];  
						owner_inc_vec[ocounter] = inc[j]; 
						fo1[ocounter] = 0.5*finc[j];
						ocounter+=1;
					}
				}

			}

			if (strcmp(exper, "c")!=0){
				printf("Aggregates (hp, ha, a, x_dhp, hrp, m, hs): %f %f %f %f %f %f %f\n", 
						hp_vec[i], ha_vec[i], a_vec[i], x_dhp_vec[i], hrp_vec[i], m_rag, hs_rag);
				save_results(i, exper);
				printf("Region %i tau-loop %i experiment %c finished!\n", i, tau_iter, *exper);
			}
		}

		//Government budget constraint
		
		aid_agg = 0.5*psi_ky*(dk_vec[0] + dk_vec[1])  +psi_kr*(dhr_vec[0] + dhr_vec[1]) 
			+0.5*psi_h* (dhp_vec[0] + dhp_vec[1]) + 0.5*psi_r*(dhr_vec[0] + dhr_vec[1]); 
		subs_agg =  0.5*eta*(ika_vec[0]+ ika_vec[1] + deltaa*(hra_vec[0]  + hra_vec[1] +ha_vec[0] + ha_vec[1])); 	
		spend = aid_agg + subs_agg; 
		rev = (tau*w_vec[0] + tau*w_vec[1])/2.0;
		if (spend > 1.0e-5) tdiff = fabs(spend - rev)/spend;
		else tdiff = fabs(spend - rev); 
		
		if (strcmp(exper, "c")!=0) printf("(tdiff, tau, tau_new): (%f %f %f)\n", tdiff, tau, spend/((w_vec[0] + w_vec[1])/2.0));
		tau = spend/((w_vec[0] + w_vec[1])/2.0);
		if (tau_iter == 5){
			tau_prob = 1; 
			tdiff =0; 
		}

	//	if ( (i ==0) && (strcmp(exper, "ns")==0)) break; //don't solve no storm experiment twice
	}

	//Distribution of owners and renters
	for (i = 0 ; i< rcounter; i ++) fr[i] = fr1[i]/renters; 
	for (i=0; i < ocounter; i++) fo[i] = fo1[i]/owners; 

	//Economy-wide aggregate (aggregates in per capita terms (measure 2 of people))
	d_vec[0] = dk_vec[0] +dh_vec[0];  
	d_vec[1] = dk_vec[1] + dh_vec[1]; 
	gdp = 0.5*(ytilde_vec[0]  + ytilde_vec[1] + ph_vec[0]*Ah*(hrp_vec[0] - dhr_vec[0] + hpt_vec[0]) +
			ph_vec[1]*Ah*(hrp_vec[1]- dhr_vec[1] + hpt_vec[1])); 
	gdp_vec[0] = ytilde_vec[0] + ph_vec[0]*Ah*(hrp_vec[0] - dhr_vec[0] + hpt_vec[0]); 
	gdp_vec[1] = ytilde_vec[1] + ph_vec[1]*Ah*(hrp_vec[1] - dhr_vec[1] + hpt_vec[1]); 
	ho_agg =  0.5*(hpt_vec[0] + hpt_vec[1] + ha_vec[0] + ha_vec[1]);
	hr_agg =  0.5*(hrp_vec[0] + hrp_vec[1] - dhr_vec[0] - dhr_vec[1] + hra_vec[0] +hra_vec[1]);
	k_agg  =0.5*(kp_vec[0] + kp_vec[1] + ka_vec[0] + ka_vec[1]- dk_vec[0] - dk_vec[1]); 
	a_agg = 0.5*(a_vec[0] + a_vec[1]); 
	d_agg = 0.5*(d_vec[0] + d_vec[1]);
	
	//Save baseline value of tax
	if (strcmp(exper, "b") ==0){
		//baseline value of tau
		rs = fopen("tau_saved", "w"); 
		if (rs==NULL) printf("The file did not open\n");
		fprintf(rs, "%12.8f\n", tau); 
		fclose(rs); 
	}

	if (strcmp(exper, "c")!=0){
		for(i=0; i <2; i++){	
			//Save aggregates
			if ((strcmp(exper, "cc") ==0) || (strcmp(exper, "ccna") ==0)|| (strcmp(exper, "ccnr") ==0)) {
				sprintf(filename, "./%s/agg_%s%i_%i.dat", folder, exper,i, ccexp );
			}
			else sprintf(filename, "./%s/agg_%s%i.dat", folder, exper,i );
			rs = fopen(filename, "w");
			if (rs==NULL) printf("The file did not open\n");
			fprintf(rs, "%12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n   %12.8f\n %12.8f\n  %12.8f\n",  
					kp_vec[i], ka_vec[i], kptilde_vec[i], dk_vec[i],  w_vec[i], hp_vec[i], 
					hrp_vec[i], ha_vec[i]);
			fprintf(rs, "%12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n %12.8f\n %12.8f\n  %12.8f\n", 
					a_vec[i], x_dhp_vec[i], owner_frac_vec[i],  hra_vec[i], tau, d_vec[i], ika_vec[i]);
			fprintf(rs, " %12.8f\n  %12.8f\n %12.8f\n  %12.8f\n %12.8f\n %12.8f\n %12.8f\n  %12.8f\n %12.8f\n", 
					  theta, Omegah_vec[i],psi_ky, eta, tdiff, ph_vec[i],
					rho_vec[i], dhr_vec[i], psi_kr);
			fprintf(rs, "%12.8f\n  %12.8f\n %12.8f\n %12.8f\n  %12.8f\n %12.8f\n %12.8f\n %12.8f\n %i\n %12.8f\n", 
					ar_vec[i], ay_vec[i], Omegak_vec[i], dhp_vec[i], ytilde_vec[i], gdp_vec[i], psi_h,psi_r,iter, mdist);
			
			fprintf(rs, "%i\n  %12.9f\n", iter_sim_peeps, mdist_sim_peeps);  
			fclose(rs);
		}
	}
	if (strcmp(exper, "e") ==0){
		top = (d_vec[1]/(rho_vec[1]*gdp_vec[1]))/(d_vec[0]/(rho_vec[0]*gdp_vec[0])) -1;
		bot = rho_vec[1]/rho_vec[0] -1; 
		printf("Elasticity of damage = %f\n", top/bot); 
		top = (aid_vec[1]/(rho_vec[1]*gdp_vec[1]))/(aid_vec[0]/(rho_vec[0]*gdp_vec[0])) -1;
		printf("Elasticity of aid = %f\n", top/bot); 
		
		top = (d_vec[1]/(rho_vec[1]*(hp_vec[1] + hrp_vec[1] + kp_vec[1])))/
				(d_vec[0]/(rho_vec[0]*(hp_vec[0] + hrp_vec[0] + kp_vec[0]))) -1;
		bot = rho_vec[1]/rho_vec[0] -1; 
		printf("Elasticity of damage = %f\n", top/bot); 
		top = (aid_vec[1]/(rho_vec[1]*(hp_vec[1] + hrp_vec[1] + kp_vec[1])))/
				(aid_vec[0]/(rho_vec[0]*(hp_vec[0] + hrp_vec[0] + kp_vec[0]))) -1;
		printf("Elasticity of aid = %f\n", top/bot); 
	}
	
	if (strcmp(exper, "b")==0 || strcmp(exper, "ht") ==0){
		
		//Moments
		wealth_output = (a_agg + ho_agg)/gdp; 
		aid_output = aid_agg/gdp; 
		subs_output = subs_agg/gdp; 
		housing_capital = (ho_agg + hr_agg)/k_agg;
		owned_total = ho_agg/(ho_agg+ hr_agg);
		rel_gdp = gdp_vec[1]/gdp_vec[0]; 
		insurance_house_damage = (owner_frac_vec[0]*x_dhp_vec[0] + owner_frac_vec[1]*x_dhp_vec[1])/
			(owner_frac_vec[0] + owner_frac_vec[1]);
		housing_damage_total_damage = 0.5*(dh_vec[0] + dh_vec[1])/
			(0.5*(d_vec[0] + d_vec[1])); 
		owner_frac = 0.5*(owner_frac_vec[0] + owner_frac_vec[1]);
		fema_damage = aid_agg/(0.5*(d_vec[0] + d_vec[1]));
		renter_aid_total_aid = psi_r*0.5*(dhr_vec[0] + dhr_vec[1])/aid_agg;
		owner_aid_total_aid = (psi_h*0.5*(dhp_vec[0] + dhp_vec[1]))/aid_agg; 
		rental_housing_aid_total_aid = (psi_kr*0.5*(dhr_vec[0] + dhr_vec[1]))/aid_agg; 
		rel_fema_aid1 = (psi_h*dhp_vec[1]/(hp_vec[1]*rho_vec[1]))/(psi_h*dhp_vec[0]/(hp_vec[0]*rho_vec[0])); 
		rel_fema_aid2 = (psi_r*dhr_vec[1]/(hrp_vec[1]*rho_vec[1]))/(psi_r*dhr_vec[0]/(hrp_vec[0]*rho_vec[0])); 
		rel_fema_aid3 = (psi_ky*dk_vec[1]/(kp_vec[1]*rho_vec[1]))/(psi_ky*dk_vec[0]/(kp_vec[0]*rho_vec[0]));
		rel_aid = (aid_vec[1]/(rho_vec[1]*(hp_vec[1] + hrp_vec[1] + kp_vec[1])))/
						(aid_vec[0]/(rho_vec[0]*(hp_vec[0] + hrp_vec[0] + kp_vec[0])));
		

		//Earnings, wealth, and income,  quintiles 		
		e1 =0.0; e2 =0.0; e3 =0.0; e4 = 0.0; e5 =0.0; 
		w1 =0.0; w2 =0.0; w3 =0.0; w4 = 0.0; w5 =0.0; 
		i1 =0.0; i2 =0.0; i3 =0.0; i4 = 0.0; i5 =0.0; 
		cse =0.0; csw =0.0; csi=0.0; 

		//Sort 
		
		for (i=0; i <2*n; i++){
			for (j =0; j <2*n; j++){
				if (wealth_pop[j] > wealth_pop[i]){
					tmp = wealth_pop[i]; 
					wealth_pop[i] = wealth_pop[j]; 
					wealth_pop[j] = tmp; 

					tmp = fw_pop[i]; 
					fw_pop[i] = fw_pop[j]; 
					fw_pop[j]  = tmp;
				}
				if (earns_pop[j] > earns_pop[i]){
					tmp = earns_pop[i]; 
					earns_pop[i] = earns_pop[j]; 
					earns_pop[j] = tmp; 

					tmp = fe_pop[i]; 
					fe_pop[i] = fe_pop[j]; 
					fe_pop[j]  = tmp;
				}
			/*	if (inc_pop[j] > inc_pop[i]){
					tmp = inc_pop[i]; 
					inc_pop[i] = inc_pop[j]; 
					inc_pop[j] = tmp; 

					tmp = fi_pop[i]; 
					fi_pop[i] = fi_pop[j]; 
					fi_pop[j]  = tmp;
				} 
			*/	
			}
		}
		if (ishares ==1){
			for (i=0; i <2*(n*nz/2); i++){
				for (j =0; j <2*(n*nz/2); j++){
					if (inc_pop[j] > inc_pop[i]){
						tmp = inc_pop[i]; 
						inc_pop[i] = inc_pop[j]; 
						inc_pop[j] = tmp; 

						tmp = fi_pop[i]; 
						fi_pop[i] = fi_pop[j]; 
						fi_pop[j]  = tmp;
					}
				}
			}
		}
		
		printf("here2\n"); 
		//Find quintiles

		for (i=0; i <2*n; i++){
			csw +=fw_pop[i]; 
			if (csw <=0.2) w1 += fw_pop[i]*wealth_pop[i];
			else if (csw > 0.2 & (csw - fw_pop[i]) < 0.2){
				tmp  = 0.2 - (csw- fw_pop[i]); 
				w1 += tmp*wealth_pop[i]; 
				w2 += (fw_pop[i] - tmp)*wealth_pop[i]; 
			}
			else if (csw <=0.4) w2 += fw_pop[i]*wealth_pop[i]; 
			else if (csw > 0.4 & (csw - fw_pop[i]) < 0.4){
				tmp  = 0.4 - (csw- fw_pop[i]); 
				w2 += tmp*wealth_pop[i]; 
				w3 += (fw_pop[i] - tmp)*wealth_pop[i]; 
			}
			else if (csw <=0.6) w3 += fw_pop[i]*wealth_pop[i]; 
			else if (csw > 0.6 & (csw - fw_pop[i]) < 0.6){
				tmp  = 0.6 - (csw- fw_pop[i]); 
				w3 += tmp*wealth_pop[i]; 
				w4 += (fw_pop[i] - tmp)*wealth_pop[i]; 
			}
			else if (csw <=0.8) w4 += fw_pop[i]*wealth_pop[i]; 
			else if (csw > 0.8 & (csw - fw_pop[i]) < 0.8){
				tmp  = 0.8 - (csw- fw_pop[i]); 
				w4 += tmp*wealth_pop[i]; 
				w5 += (fw_pop[i] - tmp)*wealth_pop[i]; 
			}
			else  w5 += fw_pop[i]*wealth_pop[i]; 
			
			cse +=fe_pop[i]; 
			if (cse <=0.2) e1 += fe_pop[i]*earns_pop[i]; 
			else if (cse> 0.2 & (cse - fe_pop[i]) < 0.2){
				tmp  = 0.2 - (cse- fe_pop[i]); 
				e1 += tmp*earns_pop[i]; 
				e2 += (fe_pop[i] - tmp)*earns_pop[i]; 
			}
			else if (cse <=0.4) e2 += fe_pop[i]*earns_pop[i]; 
			else if (cse> 0.4 & (cse - fe_pop[i]) < 0.4){
				tmp  = 0.4 - (cse- fe_pop[i]); 
				e2 += tmp*earns_pop[i]; 
				e3 += (fe_pop[i] - tmp)*earns_pop[i]; 
			}
			else if (cse <=0.6) e3 += fe_pop[i]*earns_pop[i]; 
			else if (cse> 0.6 & (cse - fe_pop[i]) < 0.6){
				tmp  = 0.6 - (cse- fe_pop[i]); 
				e3 += tmp*earns_pop[i]; 
				e4 += (fe_pop[i] - tmp)*earns_pop[i]; 
			}
			else if (cse <=0.8) e4 += fe_pop[i]*earns_pop[i]; 
			else if (cse> 0.8 & (cse - fe_pop[i]) < 0.8){
				tmp  = 0.8 - (cse- fe_pop[i]); 
				e4 += tmp*earns_pop[i]; 
				e5 += (fe_pop[i] - tmp)*earns_pop[i]; 
			}
			else  e5 += fe_pop[i]*earns_pop[i]; 
			
		/*	csi +=fi_pop[i]; 
			if (csi <=0.2) i1 += fi_pop[i]*inc_pop[i]; 
			else if (csi> 0.2 & (csi - fi_pop[i]) < 0.2){
				tmp  = 0.2 - (csi- fi_pop[i]); 
				i1 += tmp*inc_pop[i]; 
				i2 += (fi_pop[i] - tmp)*inc_pop[i]; 
			}
			else if (csi <=0.4) i2 += fi_pop[i]*inc_pop[i]; 
			else if (csi> 0.4 & (csi - fi_pop[i]) < 0.4){
				tmp  = 0.4 - (csi- fi_pop[i]); 
				i2 += tmp*inc_pop[i]; 
				i3 += (fi_pop[i] - tmp)*inc_pop[i]; 
			}
			else if (csi <=0.6) i3 += fi_pop[i]*inc_pop[i]; 
			else if (csi> 0.6 & (csi - fi_pop[i]) < 0.6){
				tmp  = 0.6 - (csi- fi_pop[i]); 
				i3 += tmp*inc_pop[i]; 
				i4 += (fi_pop[i] - tmp)*inc_pop[i]; 
			}
			else if (csi <=0.8) i4 += fi_pop[i]*inc_pop[i]; 
			else if (csi> 0.8 & (csi - fi_pop[i]) < 0.8){
				tmp  = 0.8 - (csi- fi_pop[i]); 
				i4 += tmp*inc_pop[i]; 
				i5 += (fi_pop[i] - tmp)*inc_pop[i]; 
			}
			else  i5 += fi_pop[i]*inc_pop[i];
			*/
		}

		if (ishares ==1){
			for (i=0; i <2*(n*nz/2); i++){
				csi +=fi_pop[i]; 
				if (csi <=0.2) i1 += fi_pop[i]*inc_pop[i]; 
				else if (csi> 0.2 & (csi - fi_pop[i]) < 0.2){
					tmp  = 0.2 - (csi- fi_pop[i]); 
					i1 += tmp*inc_pop[i]; 
					i2 += (fi_pop[i] - tmp)*inc_pop[i]; 
				}
				else if (csi <=0.4) i2 += fi_pop[i]*inc_pop[i]; 
				else if (csi> 0.4 & (csi - fi_pop[i]) < 0.4){
					tmp  = 0.4 - (csi- fi_pop[i]); 
					i2 += tmp*inc_pop[i]; 
					i3 += (fi_pop[i] - tmp)*inc_pop[i]; 
				}
				else if (csi <=0.6) i3 += fi_pop[i]*inc_pop[i]; 
				else if (csi> 0.6 & (csi - fi_pop[i]) < 0.6){
					tmp  = 0.6 - (csi- fi_pop[i]); 
					i3 += tmp*inc_pop[i]; 
					i4 += (fi_pop[i] - tmp)*inc_pop[i]; 
				}
				else if (csi <=0.8) i4 += fi_pop[i]*inc_pop[i]; 
				else if (csi> 0.8 & (csi - fi_pop[i]) < 0.8){
					tmp  = 0.8 - (csi- fi_pop[i]); 
					i4 += tmp*inc_pop[i]; 
					i5 += (fi_pop[i] - tmp)*inc_pop[i]; 
				}
				else  i5 += fi_pop[i]*inc_pop[i]; 
			}
		}
		
		w_tot = w1 + w2 + w3+ w4+ w5; 
		e_tot = e1+ e2 +  e3+ e4 + e5;
		i_tot = i1 + i2 + i3 + i4+ i5;

		
		//Median owner and renter income
		//sort
		for (i=0; i <rcounter; i++){
			for (j =0; j <rcounter; j++){
				if (renter_inc_vec[j] > renter_inc_vec[i]){
					tmp = renter_inc_vec[i]; 
					renter_inc_vec[i] = renter_inc_vec[j]; 
					renter_inc_vec[j] = tmp; 

					tmp = fr[i]; 
					fr[i] = fr[j]; 
					fr[j]  = tmp;
				}
			}
		}
		
		for (i=0; i <ocounter; i++){
			for (j =0; j <ocounter; j++){
				if (owner_inc_vec[j] > owner_inc_vec[i]){
					tmp = owner_inc_vec[i]; 
					owner_inc_vec[i] = owner_inc_vec[j]; 
					owner_inc_vec[j] = tmp; 

					tmp = fo[i]; 
					fo[i] = fo[j]; 
					fo[j]  = tmp;
				}
			}
		}
		
		printf("here4\n"); 
		//Find median
		csi =0.0; 
		for (i=0; i < rcounter; i++){
			csi += fr[i];
			if (csi <= 0.5) renter_med = renter_inc_vec[i]; 
			else break; 
		} 

		csi =0.0; 
		for (i=0; i < ocounter; i++){
			csi += fo[i];
			if (csi <= 0.5) owner_med = owner_inc_vec[i]; 
			else break; 
		} 
		
		if (strcmp(exper, "b")==0 || strcmp(exper, "ht") ==0){
			//Save moments and targets
			sprintf(filename, "./%s/targets_%s.dat", folder, exper );
			rs = fopen(filename, "w");
			if (rs==NULL) printf("The file did not open\n");
			fprintf(rs, "%12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n   %12.8f\n %12.8f\n  %12.8f\n",  
					wealth_output, aid_output, subs_output, housing_capital, owned_total, rel_gdp, 
					insurance_house_damage, housing_damage_total_damage);
			fprintf(rs, "%12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n   %12.8f\n %12.8f\n  %12.8f\n  %12.8f\n",  
					owner_frac, fema_damage, renter_aid_total_aid, owner_aid_total_aid, rental_housing_aid_total_aid, 
					rel_fema_aid1, rel_fema_aid2, rel_fema_aid3);

			fprintf(rs, "%12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n   %12.8f\n %12.8f\n  %12.8f\n",  
					wealth_output_target, aid_output_target, subs_output_target, housing_capital_target,
					owned_total_target, rel_gdp_target, 0.57, housing_damage_total_damage_target);
			fprintf(rs, "%12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n   %12.8f\n %12.8f\n  %12.8f\n  %12.8f\n",  
					owner_frac_target, rel_damage_target, renter_aid_total_aid_target, owner_aid_total_aid_target, 
					rental_housing_aid_total_aid_target, fema_damage_target,
					rel_damage_renters_target, rel_damage_firms_target, owner_renter_inc_target);
			fprintf(rs, "%12.8f\n%12.8f\n%12.8f\n%12.8f\n%12.8f\n%12.8f\n%12.8f\n %12.8f\n %12.8f\n %12.8f\n %12.8f\n",  
					(owner_inc/owners)/(renter_inc/renters),  
					i1/i_tot*100, i2/i_tot*100, i3/i_tot*100, i4/i_tot*100, i5/i_tot*100,  
					3.1, 8.3, 14.1, 22.6, 52.0);
			fprintf(rs, "%12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n   %12.8f\n %12.8f\n  %12.8f\n  %12.8f\n %12.8f\n",  
					w1/w_tot*100, w2/w_tot*100, w3/w_tot*100, w4/w_tot*100, w5/w_tot*100, -0.7, 0.6, 3.2, 9.8, 87.0);

			fprintf(rs, "%12.8f\n  %12.8f\n ", rel_aid_target, rel_aid);  
			fclose(rs); 
		}

		if (strcmp(exper, "b") ==0){

			printf("\n\n************Targeted Moments***********\n"); 
			printf("Relative aid  %4.5f %4.5f\n", rel_aid, rel_aid_target); 
			printf("wealth_output: %4.5f %4.5f\n", wealth_output, wealth_output_target); 
			printf("aid_output: %4.5f %4.5f\n",aid_output, aid_output_target); 
			printf("subs_output: %4.5f %4.5f\n", subs_output, subs_output_target); 
			printf("insurance_house_damage: %4.5f %4.5f\n", insurance_house_damage, 0.57); 
			printf("fema_damage: %4.5f %4.5f\n", fema_damage,fema_damage_target); 
			printf("renter_aid_total_aid: %4.5f %4.5f\n", renter_aid_total_aid, renter_aid_total_aid_target); 
			printf("owner_aid_total_aid: %4.5f %4.5f\n", owner_aid_total_aid, owner_aid_total_aid_target); 
			printf("rental_housing_aid_total_aid: %4.5f %4.5f\n", rental_housing_aid_total_aid, 
					rental_housing_aid_total_aid_target); 
			printf("Housing_damage_total_damage: %4.5f %4.5f\n", housing_damage_total_damage,
					housing_damage_total_damage_target);
			printf("housing_capital: %4.5f %4.5f\n", housing_capital, housing_capital_target); 
			printf("owned to total housing: %4.5f %4.5f\n",  owned_total, owned_total_target);
			printf("owner fraction: %4.5f %4.5f\n", owner_frac, owner_frac_target); 
			printf("Mean: owner/renter income %4.5f %4.5f\n", 
					(owner_inc/owners)/(renter_inc/renters), owner_renter_inc_target); 
			
			printf("\n\n************UnTargeted Moments***********\n"); 
			printf("rel_fema_aid_renters: %4.5f %4.5f\n", rel_fema_aid2, rel_damage_renters_target); 
			printf("rel_fema_aid_firms: %4.5f %4.5f\n", rel_fema_aid3, rel_damage_firms_target); 
			printf("rel_fema_aid_owners: %4.5f %4.5f\n", rel_fema_aid1, rel_damage_target); 
			printf("rel_gdp: %4.5f %4.5f\n", rel_gdp, rel_gdp_target); 
			printf("Earnings shares: (%4.2f %4.2f %4.2f %4.2f %4.2f) (3.1, 8.3, 14.1, 22.6, 52.0)\n", 
					e1/e_tot*100, e2/e_tot*100, e3/e_tot*100, e4/e_tot*100, e5/e_tot*100); 
			printf("Income shares: (%4.2f %4.2f %4.2f %4.2f %4.2f) (3.1, 8.3, 14.1, 22.6, 52.0)\n", 
					i1/i_tot*100, i2/i_tot*100, i3/i_tot*100, i4/i_tot*100, i5/i_tot*100); 
			printf("Wealth  shares: (%4.2f %4.2f %4.2f %4.2f %4.2f) (-0.7, 0.6, 3.2, 9.8, 87.0)\n", 
					w1/w_tot*100, w2/w_tot*100, w3/w_tot*100, w4/w_tot*100, w5/w_tot*100); 
			printf("Wealth shares: (%4.2f %4.2f)(3.3, 96.7)\n", (w1+ w2+w3)/w_tot*100, (w4+w5)/w_tot*100); 
			printf("Median: owner/renter income %4.5f %4.5f\n", owner_med/renter_med, owner_renter_inc_target);

		}

		if (strcmp(exper, "c") ==0){
			fprintf(flog, "%f, %f, %f,  %f, %f, %f,  %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", 
					theta, Omegak_Omegah, tfp_h, beta, lambda, Omega_hh, eta, hpmin, phi,chi,  
					rel_damage, aid_output, subs_output, owned_total, owner_frac, 
					housing_damage_total_damage,  wealth_output, insurance_house_damage, 
					rel_gdp,  fema_damage, renter_aid_total_aid, owner_aid_total_aid, (owner_inc/owners)/(renter_inc/renters),
					housing_capital);  
		}
	}
	
			time(&tend);
			time_spent = (double)(tend - tstart);
			printf("Time spent: %12.2f seconds, %12.2f minutes, %12.2f hours\n", 
					 time_spent, time_spent/60, time_spent/60/60);
}
