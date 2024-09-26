#include "header.h"

extern double eta, theta, deltaa,  deltay, psi_ky, r, alpha; 

double err_lb(double ay, double Omega, double rho, double tfp) 
{

	double err,  top,  bot, cons_p, F, Fprime; 

	F = exp(-theta*pow(ay, 1.0/theta)); 
	Fprime = -exp(-theta*pow(ay, 1.0/theta))*pow(ay, 1.0/theta-1.0); 
	cons_p = -ay; //da/dkp*kp 
	
	err = 1.0 - rho*Omega*(F + Fprime*cons_p); 

}
