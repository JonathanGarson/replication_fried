/*
Program: shell.c
Purpose: Runs different specifications of the model and saves results. 
Calls: solveModel.c
*/


#include "header.h"

int  *eind, *zind, *hsind, *mind, *qind, *ooh,  cal_iter, ccexp, iter, iter_sim_peeps; 
double dist, tol =0.00005,  ftol = 1.0e-14, gtol = 1.0e-12, dtol = 1.0e-9, ttol = 0.0001;
double *v0, *v1,  *f1, *f2,  *hppol, *hapol, *apol, *xpol, *hsrpol,  hsinc,  *hsgrid, *mgrid, zgrid[10], tmat[10][10], 
       *v_h_base, *v_l_base, *v_h, *v_l, *hapol_l_base, *hapol_h_base, Omega_h_base, mdist, mdist_sim_peeps, 
       Omega_hh_base,  ph_vec_base[2], *earns, *wealth,  *fe_pop, *fw_pop,  *earns_pop, *wealth_pop, *inc_pop, *inc, *fi_pop, *xfgrid, 
       *hppol_own, *hapol_own, *apol_own, *xpol_own, *owner_inc_vec, *renter_inc_vec, *fr1, *fo1, *fo, *fr, *ahpol, *finc, **ya, **y2a, 
       *cont_e0, *cont_e1; 

char *folder = "results_1109";  

//Grid parameters
double hsmax =30.0, hsmin= 0.000001, mmax=200.0, mmin = 0.001; 
int nhs =68, nm =60, nz =10, ne=2, nq = 2;  

// External Parameter values
double tau= 0.00105;
double deltay = 0.068;
double deltah = 0.023; 
double deltaa = 0.025; 
double sigma =2.0;
double rho_l =0.0021; 
double rho_h = 0.022; 
double r = 0.04; 
double tfp_l = 1.0; 
double alpha =0.26; 
double Omegah_Omegal = 1.98; 
double Ah =0.14; 
double tfp_h = 1.00;

//Internal parameter values 
double zeta = 0.8491; 
double psi_ky = 0.23;
double psi_kr = 0.015;
double psi_h = 0.1;
double psi_r = 0.034;
double theta= 3.78;
double beta = 0.944;
double lambda = 0.00512;
double Omega_hh = 0.3; //severity of storms for households in the high-risk region 
double eta = 0.086;
double Omegak_Omegah = 0.7; // Difference in storm severity between housing and non-housing capital 
double hpmin = 1.15; 
double phi = 0.79; //discount from renting
double chi = 0.789; //Fraction that discount renting (phi <1) increase chi => increase fraction of homeowners

int use_saved =1; //set to 1 if use existing starting guess for  value function

//Targeted moments
double rel_aid_target = 1.16; 
double rel_damage_target = 1.19; 
double rel_damage_renters_target = 0.97; 
double rel_damage_firms_target = 1.22; 
double subs_output_target = .0001037, aid_output_target = .00058;
double renter_aid_total_aid_target = 0.03,  owner_aid_total_aid_target =0.17, rental_housing_aid_total_aid_target = 0.01;
double rel_gdp_target = 0.98;
double fema_damage_target =0.17; 
double wealth_output_target = 3.08; 
double owner_frac_target =0.67, owned_total_target = 0.78; //housing target
double housing_damage_total_damage_target = 0.43;
double housing_capital_target = 0.88;
double owner_renter_inc_target = 1.5; 
double med_owner_renter_inc_target = 1.38; 

int main() {
	//Flags (set equal to 1 to run particular experiment)
	int baseline =1; //=1: solves baseline model and calculates moments
	int no_policy =0;  //Sets fema and and adaptation subsidy to zero
	int no_subsidy =0; //Sets adaptation subsidty to zerom
	int cc =0; //runs climate change experiments
	int no_risk =0; //runs no risk experiment 
	int higher_rho = 0;  //Changes probability of a storm to compute elasticities
	
	//Calculate Omega_h
	//Define local variables
	FILE *rs, *flog, *fcal;
	int i, n, n2, ri, ci, j, nrows, ncols;
	double theta_base, eta_base, psiky_base, psikr_base,  psih_base, tau_base, rho_l_base, rho_h_base,psir_base, 
	       minc1, minc2, minc, temp, 
	       zgrid1[nz/2], tmat1[nz/2][nz/2], zgrid2[nz/2], Omegah_Omegal_base, sfrho_l[6], sfrho_h[6], sfomega_l[6],
	       sfomega_h[6]; 
	char filename[100], *exper;
	mkdir(folder, 0700);

	//Save basleine values of parameters 
	theta_base = theta; rho_l_base  =rho_l; rho_h_base = rho_h; Omega_hh_base =Omega_hh;  
	eta_base = eta; psiky_base =psi_ky; psih_base = psi_h; psikr_base = psi_kr;  psir_base = psi_r; Omegah_Omegal_base = Omegah_Omegal;  

	n = nz*nhs*nm*nq;
	n2 = 2*n;

	//open log file
	sprintf(filename, "./%s/results1.doc", folder);
	flog = fopen(filename, "w");

	//Assign memory 
	hsgrid = (double *)malloc(sizeof(double)*nhs);
	mgrid = (double *)malloc(sizeof(double)*nm); 
	hppol = (double *)malloc(sizeof(double)*n);
	hapol = (double *)malloc(sizeof(double)*n);
	apol = (double *)malloc(sizeof(double)*n);
	xpol = (double *)malloc(sizeof(double)*n);
	hsrpol = (double *)malloc(sizeof(double)*n);
	hppol_own = (double *)malloc(sizeof(double)*n);
	hapol_own = (double *)malloc(sizeof(double)*n);
	xpol_own = (double *)malloc(sizeof(double)*n);
	apol_own = (double *)malloc(sizeof(double)*n);
	v0 = (double *)malloc(sizeof(double)*n);
	v1 = (double *)malloc(sizeof(double)*n);
	v_h_base = (double *)malloc(sizeof(double)*n);
	v_l_base = (double *)malloc(sizeof(double)*n);
	v_h = (double *)malloc(sizeof(double)*n);
	v_l = (double *)malloc(sizeof(double)*n);
	f1 = (double *)malloc(sizeof(double)*n);
	f2 = (double *)malloc(sizeof(double)*n);
	earns = (double *)malloc(sizeof(double)*n);
	wealth = (double *)malloc(sizeof(double)*n);
	ahpol = (double *)malloc(sizeof(double)*n);
	hsind = (int *)malloc(sizeof(int)*n);
	mind = (int *)malloc(sizeof(int)*n);
	zind = (int *)malloc(sizeof(int)*n); 
	qind = (int *)malloc(sizeof(int)*n); 
	ooh = (int *)malloc(sizeof(int)*n*nz/2); 
	finc = (double *)malloc(sizeof(double)*n*nz/2); 
	inc = (double *)malloc(sizeof(double)*n*nz/2);
	cont_e0 = (double *)malloc(sizeof(double)*n);
	cont_e1 = (double *)malloc(sizeof(double)*n);
	
	hapol_l_base = (double *)malloc(sizeof(double)*n);
	hapol_h_base = (double *)malloc(sizeof(double)*n);

	//Vector for the whole population
	fe_pop = (double *)malloc(sizeof(double)*n2);
	fi_pop = (double *)malloc(sizeof(double)*n2*nz/2);
	fw_pop = (double *)malloc(sizeof(double)*n2);
	fr1 = (double *)malloc(sizeof(double)*n2*nz/2);
	fr = (double *)malloc(sizeof(double)*n2*nz/2);
	fo1 = (double *)malloc(sizeof(double)*n2*nz/2);
	fo = (double *)malloc(sizeof(double)*n2*nz/2);
	renter_inc_vec = (double *)malloc(sizeof(double)*n2*nz/2);
	owner_inc_vec = (double *)malloc(sizeof(double)*n2*nz/2);
	earns_pop = (double *)malloc(sizeof(double)*n2);
	wealth_pop = (double *)malloc(sizeof(double)*n2);
	inc_pop = (double *)malloc(sizeof(double)*n2*nz/2);
	
	//allocate array for cubic splines
	nrows = nhs*nq*nz; ncols = nm*nq*nz;
	ya = malloc(nrows* sizeof(double*));
	y2a = malloc(nrows* sizeof(double*));
	if (ya==NULL) printf("out of memory");  
	if (y2a==NULL) printf("out of memory");  
	for (i =0; i <nrows; i++){
		ya[i]  = malloc(ncols * sizeof(double)); 
		y2a[i]  = malloc(ncols * sizeof(double)); 
		if (ya[i] == NULL) printf("out of memory"); 
		if (y2a[i] == NULL) printf("out of memory"); 
	}


	minc = (log(mmax) - log(mmin))/(nm-1.0); 
	hsinc = (log(hsmax) - log(hsmin))/(nhs-1.0);

	for (i = 0; i < nhs; i++) {
		temp= log(hsmin) + i *hsinc;
		hsgrid[i] = exp(temp); 
	}
	
	for (i =0; i < nm; i ++) {
		temp= log(mmin) + i*minc; 
		mgrid[i] =exp(temp); 
	}

	//Create indicies
	for (i = 0; i < n; i++) { //(q, z,  hs, m)
		qind[i] = ((int)i / (nhs*nm*nz));
		zind[i] = ((int)i / (nhs*nm)) %nz;
		hsind[i] = ((int)i / nm) % nhs;
		mind[i] = (int) i % nm; 
	}
	
	
	//Earnings shocks, Kaplan estimates (1)
	zgrid[0] = 0.1017; zgrid[1] = 0.1711; zgrid[2] = 0.2879; zgrid[3] = 0.4843; zgrid[4] = 0.8149; zgrid[5] = 0.5163; 
	zgrid[6] = 0.8687; zgrid[7]= 1.4616; zgrid[8] = 2.4592; zgrid[9] = 4.1377; 

	//Transition matrix (rho =0.97) Kaplan estimates(1)
	tmat1[0][0] = 0.94133655; tmat1[0][1] = 0.0573403; tmat1[0][2] = 0.0013098; tmat1[0][3] = 0.0000133; tmat1[0][4] = 0.00000005;
	tmat1[1][0] = 0.01433507; tmat1[1][1] = 0.94199145; tmat1[1][2]= 0.0430152; tmat1[1][3] = 0.00065495; tmat1[1][4] = 0.00000332; 
	tmat1[2][0] = 0.0002183; tmat1[2][1] = 0.0286768; tmat1[2][2] = 0.9422098; tmat1[2][3] = 0.0286768;tmat1[2][4] = 0.0002183; 
	tmat1[3][0] =  0.00000332; tmat1[3][1] = 0.00065495;  tmat1[3][2] = 0.0430152; tmat1[3][3] = 0.94199145; tmat1[3][4] = 0.01433507;
	tmat1[4][0] = 0.00000005; tmat1[4][1]= 0.0000133; tmat1[4][2]= 0.0013098; tmat1[4][3] = 0.0573403; tmat1[4][4] = 0.94133655; 


	
	for (ri =0; ri< nz; ri++){ //current state
		for(ci=0; ci < nz; ci++){ //future state
			if (ri <nz/2 & ci < nz/2) tmat[ri][ci] = tmat1[ri][ci]; 
			else if (ri >=nz/2 & ci >=nz/2) tmat[ri][ci] = tmat1[ri-nz/2][ci-nz/2]; 
			else tmat[ri][ci] = 0.0; 
		}
	}

/*
	//Check at columns sum to one 
	for (i = 0; i < nz; i++){
		printf("%12.10f\n", tmat[i][4] +  tmat[i][0] + tmat[i][1] + tmat[i][2] + tmat[i][3]+ tmat[i][5] + tmat[i][6] + tmat[i][7]
				+ tmat[i][8] + tmat[i][9]); //+ tmat[i][10] +tmat[i][11] + tmat[i][12] + tmat[i][13]); 
	}
	
*/	//Form region vectors (low risk, high risk)

	fprintf(flog, "*********External values***********\n\n"); 
	fprintf(flog, "VFI tolerance %f\n", tol); 

	fprintf(flog, "\n\n\n"); 
	fprintf(flog, "*********Internally calibrated parameter values***********\n\n"); 
	fprintf(flog, "beta %f\n", beta); 
	fprintf(flog, "lambda %f\n", lambda); 
	fprintf(flog, "theta %f\n", theta); 
	fprintf(flog, "eta %f\n", eta); 
	fprintf(flog, "zeta %f\n", zeta); 
	fprintf(flog, "Omega_hh %f\n", Omega_hh); 
	fprintf(flog, "Omegak_Omegah %f\n", Omegak_Omegah); 
	fprintf(flog, "\n\n\n"); 

	//save grid and parameter values results

	//Grid values
	sprintf(filename, "./%s/grids.dat", folder);
	rs = fopen(filename, "w");
	if (rs==NULL) printf("The file did not open\n");
	fprintf(rs, "%i\n %i\n  %i\n  ", nhs, nm,   nz); 
	fprintf(rs, "%12.8f\n  %12.8f\n %12.8f\n  %12.8f\n ", hsmax, hsmin, mmax, mmin ); 

	fclose(rs);

	//Parameter values
	sprintf(filename, "./%s/params.dat", folder);
	rs = fopen(filename, "w");
	if (rs==NULL) printf("The file did not open\n");
	fprintf(rs, "%12.8f\n %12.8f\n   %12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n  %12.8f\n %12.8f %12.8f\n %12.8f\n\n ",
			deltay, zeta, beta, theta, Omega_hh,Omegak_Omegah,  sigma, psi_ky, psi_h,  rho_l, rho_h);
	fprintf(rs, "%12.8f\n %12.8f\n %12.8f\n %12.8f\n  %12.10f\n %12.10f\n %12.10f\n %12.10f\n", 
			eta, lambda, deltah,  r,   tol, tfp_l, tfp_h,  alpha);
	fprintf(rs, "%12.10f\n %12.10f\n %12.10f\n %12.10f\n %12.10f\n %12.10f\n %12.10f\n %12.10f\n  %12.10f\n %12.10f\n" 
			, dtol, gtol, ftol,phi, hpmin, psi_r, Ah, deltaa, chi, psi_kr);
	fclose(rs);

	
	//tmat 
	sprintf(filename, "./%s/tmat.dat", folder);
	rs = fopen(filename, "w");
	if (rs == NULL)	printf("The file did not open\n");
	for (i = 0; i < nz; i++) { 
		for (j =0; j < nz; j++){
			fprintf(rs, "%12.8f\n ", tmat[i][j]);
		}
	}
	fclose(rs);
	
	//zgrid
	sprintf(filename, "./%s/zgrid.dat", folder);
	rs = fopen(filename, "w");
	if (rs == NULL)	printf("The file did not open\n");
	for (i = 0; i < nz; i++) fprintf(rs, "%12.8f\n ", zgrid[i]);
	fclose(rs);
	
	//hsgrid
	sprintf(filename, "./%s/hsgrid.dat", folder);
	rs = fopen(filename, "w");
	if (rs == NULL)	printf("The file did not open\n");
	for (i = 0; i < nhs; i++) fprintf(rs, "%12.8f\n ", hsgrid[i]);
	fclose(rs);
	
	//mgrid
	sprintf(filename, "./%s/mgrid.dat", folder);
	rs = fopen(filename, "w");
	if (rs == NULL)	printf("The file did not open\n");
	for (i = 0; i < nm; i++) fprintf(rs, "%12.8f\n ", mgrid[i]);
	fclose(rs);


	//Experiments 
	if (baseline == 1)
	{
		printf("\n\n**********Baseline Experiment***********\n\n"); 
		fprintf(flog, "\n\n*********Baseline Experiment***********\n\n"); 
		solveModel("b", flog);
		tau_base =tau; 
	}
	
	if (no_policy ==1){
		tau = 0.0; 
		rho_l = rho_l_base; 
		rho_h = rho_h_base; 
		Omega_hh = Omega_hh_base; 
		Omegah_Omegal = Omegah_Omegal_base; 
		psi_ky =0.0; psi_kr = 0.0; psi_h = 0.0; psi_r =0.0;  eta = 0.0; 
		printf("\n\n********* No Policy Experiment   ***********\n\n"); 
		fprintf(flog, "\n\n*********No policy  experiment  ***********\n\n"); 
		solveModel("np", flog);
	}

	if (no_subsidy ==1){
		rho_l = rho_l_base; 
		rho_h = rho_h_base; 
		Omega_hh = Omega_hh_base; 
		Omegah_Omegal = Omegah_Omegal_base; 
		psi_ky =psiky_base; psi_kr = psikr_base; psi_h = psih_base; psi_r = psir_base; 
		eta = 0.0;
		tau = tau_base; 
		printf("\n\n********* No Subsidy Experiment  ***********\n\n"); 
		fprintf(flog, "\n\n*********No Subsidy Experiment  ***********\n\n"); 
		solveModel("ns", flog);
	}

	sfrho_l[0]= 1.0; sfrho_l[1] = 1.0; sfrho_l[2] = 1.77; sfrho_l[3] = 2.44; sfrho_l[4] = 1.5; sfrho_l[5] = 2.5; 
	sfrho_h[0]= 1.0; sfrho_h[1] = 1.0; sfrho_h[2] = 1.57; sfrho_h[3] = 2.03; sfrho_h[4] = 1.5; sfrho_h[5] = 2.5; 
	sfomega_l[0]= 1.5; sfomega_l[1] = 2.0; sfomega_l[2] = 1.09; sfomega_l[3] = 1.15; sfomega_l[4] = 1.0; sfomega_l[5] = 1.0; 
	sfomega_h[0]= 1.5; sfomega_h[1] = 2.0; sfomega_h[2] = 1.09; sfomega_h[3] = 1.15; sfomega_h[4] = 1.0; sfomega_h[5] = 1.0; 
	
	if (cc ==1){
		printf("Baseline values: %f %f %f\n", Omega_hh_base, rho_l_base, rho_h_base); 
		//Climate change experiments: cyclones, extreme-precip, floods
		printf("\n\n*********Constant  adaptataion baseline experiment ***********\n\n"); 
		ccexp = 9; 
		solveModel("ccna", flog);
		
		for (i = 0; i < 6; i ++){
			if (i ==0) tau = tau_base;
			tau = 0.00140397; 
			rho_l = rho_l_base*sfrho_l[i]; 
			rho_h = rho_h_base*sfrho_h[i]; 
			Omega_hh = Omega_hh_base*sfomega_h[i]; 
			Omegah_Omegal = Omegah_Omegal_base*sfomega_h[i]/sfomega_l[i];
			psi_ky =psiky_base; psi_kr = psikr_base;  psi_h= psih_base; psi_r = psir_base; eta =eta_base; 
			ccexp = i; 
			printf("\n\n*********CC No adaptataion  experiment %i ***********\n\n", i); 
			solveModel("ccna", flog);
			printf("\n\n*********CC  experiment %i***********\n\n", i); 
			fprintf(flog, "\n\n*********CC experiment %i ***********\n\n", i); 
			solveModel("cc", flog); 
		}
	}
	if (no_risk ==1){
		rho_l = 1.0; 
		rho_h = 1.0; 
		Omega_hh = rho_h_base*Omega_hh_base; 
		Omegah_Omegal = rho_h_base/rho_l_base*Omegah_Omegal_base; 
		tau = tau_base;
		tau = 0.0014; 
		eta = eta_base; 
		psi_ky =psiky_base; psi_kr = psikr_base; psi_h = psih_base; psi_r = psir_base; 
		printf("\n\n********* No Risk Experiment  ***********\n\n"); 
		fprintf(flog, "\n\n*********No Risk Experiment  ***********\n\n"); 
		solveModel("nr", flog);

		for (i =0; i < 6; i++){
			rho_l = rho_l_base*sfrho_l[i]; 
			rho_h = rho_h_base*sfrho_h[i]; 
			Omega_hh = Omega_hh_base*sfomega_h[i]*rho_h; 
			Omegah_Omegal = Omegah_Omegal_base*sfomega_h[i]/sfomega_l[i]*rho_h/rho_l;
			rho_h = 1.0; 
			rho_l = 1.0;
			psi_ky =psiky_base; psi_kr = psikr_base;  psi_h= psih_base; psi_r = psir_base; eta =eta_base; 
			ccexp = i; 
			printf("\n\n*********CC No risk  experiment %i***********\n\n", i); 
			fprintf(flog, "\n\n*********CC No risk experiment %i ***********\n\n", i); 
			solveModel("ccnr", flog); 
		}
	}

	if (higher_rho ==1){
		rho_l = 1.01*rho_l_base; 
		rho_h = 1.01*rho_h_base; 
		Omega_hh = Omega_hh_base; 
		Omegah_Omegal = Omegah_Omegal_base; 
		psi_ky =psiky_base; psi_kr = psikr_base; psi_h = psih_base; psi_r = psir_base; 
		eta = eta_base;
		tau = tau_base;
		printf("\n\n********* Higher rho experiment  ***********\n\n"); 
		solveModel("hr", flog);
	}
	fclose(flog); 
	return 0;
}
