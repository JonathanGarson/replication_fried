/*
Program: err_kakp.c 
Purpose: solve the final-good firm's optimization problem from the first order conditions. 
*/

#include "header.h"

extern double eta, theta, deltay, deltaa,  psi_ky, r, alpha; 

double err_kakp(double ay, double Omega, double rho, double tfp) 
{

	double err,  kp, w,focka, F, Fprime, cons_p, cons_a,  top, bot; 

	F = exp(-theta*pow(ay, 1.0/theta)); 
	Fprime = -exp(-theta*pow(ay, 1.0/theta))*pow(ay, 1.0/theta-1.0); 
	cons_p = -ay; //krp*dar/dkrp, productive capital foc
	cons_a = 1.0; //adaptive capital foc 

	//kp
	top = 1.0 - rho*Omega*(F + Fprime*cons_p); 
	bot = r + deltay + rho*Omega*(1.0 - deltay - psi_ky)*(F + Fprime*cons_p); 
	kp = pow(top/bot, 1.0/(1.0 - alpha))*(1.0/(1.0 - rho*Omega*F))*pow(alpha*tfp, 1.0/(1.0-alpha)); 
	
	//wages
	w = tfp*(1.0-alpha)*pow(1.0-rho*Omega*F, alpha)*pow(kp, alpha); 

	//Foc Ka
	focka = -(1.0 + eta)*rho*Omega*(alpha*pow(tfp, 1.0/alpha)*pow((1.0-alpha)/w, (1.0-alpha)/alpha) + 1.0 - deltay - psi_ky)
		*Fprime*cons_a - r - deltaa; 

	err = focka; 
	return err;
}

