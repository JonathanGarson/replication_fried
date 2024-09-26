#include "header.h"

extern int nloc, *eivec_h, *eivec_l;
extern double delta, sigma, rho, tau, s, psi, zeta, lambda, lambdag, gdpshare_h, gdpshare_l, p_h, p_l, A_h, A_l;
extern double fema_l_targ, fema_h_targ, capital_output_target, insurance_target, subsidy_gdp_target, frac_subs_target, fema_gdp_target, gdppc_h_l_target, ratio_target;
extern double *iavec_h, *hvec_h, *dvec_h, *qvec_h, *iclaimvec_h, *iavec_l, *hvec_l, *dvec_l, *qvec_l, *iclaimvec_l, *xvec_h, *xvec_l, *ipvec_h, *ipvec_l;
extern double  *kavec_h, *kavec_l, *kptildevec_h, *kptildevec_l, *dkvec_h, *dkvec_l, *kpvec_h, *kpvec_l;


void moments(FILE *flog)
{

	double fema_y_sum_h, fema_y_sum_l, k_gdp, k_sum_h, k_sum_l, ytilde_sum_h, ytilde_sum_l, x_d_avg, k, xsum_h, xsum_l, gdp, prod_sum_h, prod_sum_l, ip_sum_h, ip_sum_l;
	double x_d_sum_l, x_d_sum_h, kptilde, ytilde, fema, femay_avg_h, femay_avg_l, ka_sum_h, ka_sum_l, ia_sum_h, ia_sum_l, ka_agg, x_sum_h_2, x_sum_l_2, d_sum_h, d_sum_l;
	double kp_sum_l, kp_sum_h, subs_sum_h, subs_sum_l, fema_sum_h, fema_sum_l, fema_sum, subs_sum, frac_fema, fema_gdp, subs_gdp, I_gdp, ratio;
	int num_storm_h, num_storm_l, loc;


	k_sum_h = 0.0;
	ytilde_sum_h = 0.0;
	x_d_sum_h = (double) 0.0;
	num_storm_h = 0;
	fema_y_sum_h = (double) 0.0;
	k_sum_l = (double) 0.0;
	ytilde_sum_l = (double) 0.0;
	x_d_sum_l = (double) 0.0;
	num_storm_l = 0;
	fema_y_sum_l = (double) 0.0;
	ka_sum_h = 0.0;
	ka_sum_l = 0.0;
	xsum_h = 0.0;
	xsum_l = 0.0;
	kp_sum_l = 0.0;
	kp_sum_h = 0.0;
	ia_sum_h = 0.0;
	ia_sum_l = 0.0;
	ip_sum_h = 0.0;
	ip_sum_l = 0.0;
	prod_sum_l = 0.0;
	prod_sum_h = 0.0;
	fema_sum_h = 0.0; 
	subs_sum_h = 0.0; 
	fema_sum_l = 0.0;
	subs_sum_l = 0.0;
	d_sum_l = 0.0; 
	x_sum_l_2 = 0.0;
	d_sum_h = 0.0; 
	x_sum_h_2 = 0.0;

	for (loc = 0; loc < nloc; loc++) {
		//high risk localities
		xsum_h += xvec_h[loc];
		ka_sum_h += kavec_h[loc];
		kp_sum_h += kpvec_h[loc];
		k = kpvec_h[loc] + kavec_h[loc];
		kptilde = kpvec_h[loc] - dkvec_h[loc];
		prod_sum_h += A_h*pow(kptilde, zeta);
		ytilde = A_h*pow(kptilde, zeta) + dvec_h[loc] * psi; //locality-level income
		k_sum_h += k;
		ytilde_sum_h += ytilde;
		ia_sum_h += iavec_h[loc];
		ip_sum_h += ipvec_h[loc];
		fema_sum_h += psi * dvec_h[loc]*(1+lambdag);
		subs_sum_h += s * iavec_h[loc];

		if (eivec_h[loc] > 0) {
			fema = psi * dvec_h[loc]*(1+lambdag);
			fema_y_sum_h += (fema / ytilde);
			num_storm_h += 1;
			x_sum_h_2 += xvec_h[loc];
			d_sum_h += dvec_h[loc];
			x_d_sum_h += xvec_h[loc] / dvec_h[loc];
			// printf("highrisk: %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n", dvec_h[loc], ytilde, fema, fema / ytilde, xvec_h[loc]/dvec_h[loc], ytilde -xvec_h[loc]);
		}


		//Low risk localities
		xsum_l += xvec_l[loc];
		ka_sum_l += kavec_l[loc];
		kp_sum_l += kpvec_l[loc];
		k = kpvec_l[loc] + kavec_l[loc];
		kptilde = kpvec_l[loc] - dkvec_l[loc];
		prod_sum_l += A_l*pow(kptilde, zeta);
		ytilde = A_l*pow(kptilde, zeta) + dvec_l[loc] * psi;
		k_sum_l += k;
		ytilde_sum_l += ytilde;
		ia_sum_l += iavec_l[loc];
		ip_sum_l += ipvec_l[loc];
		fema_sum_l += psi * dvec_l[loc]*(1+lambdag);
		subs_sum_l += s * iavec_l[loc];

		if (eivec_l[loc] > 0) {
			fema = psi * dvec_l[loc]*(1+lambdag);
			fema_y_sum_l += (fema / ytilde);
			num_storm_l += 1;
			x_sum_l_2 += xvec_l[loc];
			d_sum_l += dvec_l[loc];
			x_d_sum_l += xvec_l[loc] / dvec_l[loc];
		}
	}

	femay_avg_h = fema_y_sum_h / num_storm_h;
	femay_avg_l = fema_y_sum_l / num_storm_l;
	ratio = femay_avg_h / femay_avg_l;
	fema_sum = fema_sum_h + fema_sum_l;
	subs_sum = subs_sum_h + subs_sum_l;
	frac_fema = fema_sum/ (fema_sum + subs_sum);
	gdp = prod_sum_h + prod_sum_l;
	k_gdp = (k_sum_h + k_sum_l) / gdp;
	fema_gdp = fema_sum*(1+lambdag)/ gdp;
	subs_gdp = subs_sum / gdp;
	x_d_avg = (x_d_sum_h / num_storm_h +  x_d_sum_l / num_storm_l)/2;
	I_gdp = (ip_sum_h + ia_sum_h + ip_sum_l + ia_sum_l) / gdp; 

	printf("\n");
	printf("**********Targeted moments***************\n");
	printf("ratio_target=%2.4f, ratio = %2.4f\n", ratio_target, ratio);
	printf("k_gdp_target = %2.3f k_gdp = %2.3f\n", capital_output_target, k_gdp);
	printf("I_gdp_target = 0.255 I_gdp = %2.3f\n",  I_gdp);
	printf("x_target =  %2.3f. x_d_avg = %1.3f\n", insurance_target, x_d_avg);
	printf("gdppc_h_l_target =  %2.3f, gdp_h_l = %2.3f\n", gdppc_h_l_target, prod_sum_h/prod_sum_l);
	printf("s_gdp_target = %2.6f, s_gdp = %2.6f\n", subsidy_gdp_target, subs_gdp);
	printf("fema_gdp_target = %2.6f, fema_gdp = %2.6f\n", fema_gdp_target, fema_gdp);
	printf("\n");
	printf("**********Not Targeted moments***************\n");
	printf("frac_subsidy_target =  %2.3f frac_subsidy = %2.3f\n", frac_subs_target, 1-frac_fema);
	printf("femay_avg_l_target=%2.4f, femay_avg_l = %2.4f\n", fema_l_targ, femay_avg_l);
	printf("fema_avg_h_target = %2.4f, femay_avg_h = %2.4f\n", fema_h_targ, femay_avg_h);
	printf("\n");
	printf("*********Aggregates***************\n");
	printf("x/d: h = %12.8f  \n", x_sum_h_2 / d_sum_h);
	printf("x/d: l = %12.8f\n", x_sum_l_2 / d_sum_l);
	printf("x_avg_h = %12.8f\n", xsum_h / nloc);
	printf("x_avg_l = %12.8f\n", xsum_l / nloc);
	printf("ka_avg_h = %12.10f\n", ka_sum_h / nloc);
	printf("ka_avg_l = %12.10f\n", ka_sum_l / nloc);
	printf("kp_avg_h = %12.10f\n", kp_sum_h / nloc);
	printf("kp_avg_l = %12.10f\n", kp_sum_l / nloc);
	printf("k_avg_h = %12.10f\n", k_sum_h / nloc);
	printf("k_avg_l = %12.10f\n", k_sum_l / nloc);
	printf("num_storm_l = %i, num_storm_h = %i\n", num_storm_l, num_storm_h);
	printf("\n");


	fprintf(flog, "femay_avg_l_target=%12.8f, femay_avg_l = %12.8f\n", fema_l_targ, femay_avg_l);
	fprintf(flog, "fema_avg_h_target = %12.8f, femay_avg_h = %12.8f\n", fema_h_targ, femay_avg_h);
	fprintf(flog, "k_gdp_target = %12.8f k_gdp = %12.8f\n", capital_output_target, k_gdp);
	fprintf(flog, "I_gdp_target = 0.255 I_gdp = %12.3f\n", I_gdp);
	fprintf(flog, "x_target =  %12.8f. x_d_avg = %12.8f\n", insurance_target, x_d_avg);
	fprintf(flog, "frac_fema_target =  %12.8f frac_fema = %12.8f\n", frac_subs_target, 1- frac_fema);
	fprintf(flog, "s_gdp_target = %12.8f, s_gdp = %12.8f\n", subsidy_gdp_target, subs_gdp);
	fprintf(flog, "fema_gdp_target = %12.8f, fema_gdp = %12.8f\n", fema_gdp_target, fema_gdp);
	fprintf(flog, "x_d_avg_h = %12.8f\n", x_d_sum_h / nloc);
	fprintf(flog, "x_d_avg_l = %12.8f\n", x_d_sum_l / nloc);
	fprintf(flog, "x_avg_h = %12.8f\n", xsum_h / nloc);
	fprintf(flog, "x_avg_l = %12.8f\n", xsum_l / nloc);
	fprintf(flog, "ka_avg_h = %12.10f\n", ka_sum_h / nloc);
	fprintf(flog, "ka_avg_l = %12.10f\n", ka_sum_l / nloc);
	fprintf(flog, "kp_avg_h = %12.10f\n", kp_sum_h / nloc);
	fprintf(flog, "kp_avg_l = %12.10f\n", kp_sum_l / nloc);
	fprintf(flog, "k_avg_h = %12.10f\n", k_sum_h / nloc);
	fprintf(flog, "k_avg_l = %12.10f\n", k_sum_l / nloc);
	fprintf(flog, "num_storm_l = %i, num_storm_h = %i\n", num_storm_l, num_storm_h);
	fprintf(flog, "\n");

}
