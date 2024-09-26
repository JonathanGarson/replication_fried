/*
Program: err_hrahrp.c
Purpose: Solve the rental-houseing firm's optimization problem from the first orderconditions
*/

#include "header.h"

extern double eta, theta, deltah,  deltaa, psi_kr,  r, alpha, Ah; 

double err_hrahrp(double ar, double Omega, double rho, double tfp) 
//need to pass tfp so that you can call same rtbis function
{

	double top, bot, ph, F, Fprime, cons_p, cons_a,  err, fochra; 
	
	F = exp(-theta*pow(ar, 1.0/theta)); 
	Fprime = -exp(-theta*pow(ar, 1.0/theta))*pow(ar, 1.0/theta-1.0); 
	cons_p = -ar; //krp*dar/dkrp, productive capital foc
	cons_a = 1.0; //adaptive capital foc

	//ph
	top = r+deltah + rho*Omega*(1.0 - deltah - psi_kr)*(F + Fprime*cons_p); 
	bot = Ah*(1.0 - rho*Omega*(F + Fprime*cons_p)); 
	ph = top/bot; 

	//foc for adaptive  rental housing capital
	fochra = - rho*(1.0+eta)*(ph*Ah + 1.0 - deltah - psi_kr)*Omega*Fprime*cons_a - r - deltaa;  
	
	err = fochra;

	return err;
}
