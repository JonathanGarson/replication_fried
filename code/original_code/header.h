#pragma once

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

//Prototype functions
double err_kakp(double ka_k, double Omega, double rho, double tfp);
double rtbis(double (*func)(double ka_kp, double Omega, double rho, double tfp), 
		double x1, double x2, double  xacc, double Omega, double rho, double tfp);
void save_results(int region, char *exper);
void  golden (double rho, double ax, double bx, double cx, 
	double(*f)(double rho, double x, int hpi, int hai,  int ai, int zi, int ub), 
	double tol, double *xmin, int hpi, int hai, int ai, int zi);
double err_hrahrp(double ph, double Omega, double rho, double tfp);
void solveModel(char *exper, FILE *flog);
void sim_peeps(double rho, double Omega, double eta, double w, double hra_hr, double ph, 
		double *hp_rag, double *ha_rag, double *a_rag, 
		double *dhp_rag, double *m_rag, double *x_dhp_rag, double *hs_rag, double *hpt_rag, double *hrp_rag, double *frac_own); 
void nelder(double prime1[], double prime2[],  double ytemp, int ndim, double *hpprime, double *haprime, double *aprime, 
		double *xprime, double *hrpprime, double *ystar, double m, double hs, int zi, double rho, double Omega, 
		double  w, double hra_hr, double ph, double ah,  int own,  int qi, int na, int pr) ;
double value_function_loop( double rho, double Omega, double eta, double w, double hra_hr, double ph, int na,  int initial, int region);
double ucont(double choice[],  double m, double hs, double ah,  int zi, double rho, double  Omega, double w, double hra_hr, 
		double ph, int own,  int qi, int na, int pr);
void amoeba(double p[5][4], double y[], int ndim, double ftol, 
		double (*funk)(double[], double m, double hs, double ah, int zi, double rho, double Omega, double w,
			double hra_hr, double ph, int own,  int qi, int na, int pr),
	 int *nfunk,  double m, double hs, int zi, double rho, double Omega, double w, double hra_hr, double ph, double  ah,
	 int own, int qi, int na, int pr);
double amotry(double p[5][4], double y[], double psum[], int ndim, 
		double (*funk)(double[], double m, double hs, double ah,  int zi, double rho, double Omega, double w,
			double hra_hr, double ph, int own, int qi, int na,  int pr),
		int ihi, double fac, double m, double hs,  int zi, double rho, double Omega, double w, double hra_hr,
		double ph, double ah,  int own, int qi, int na, int pr);
double amotry_r(double p[3][2], double y[], double psum[], int ndim, 
		double (*funk)(double[], double m, double hs, double ah,   int zi, double rho, double Omega, double w,
			double hra_hr, double ph, int own, int qi, int na, int pr),
		int ihi, double fac, double m, double hs,  int zi, double rho, double Omega, double w, double hra_hr,
		double ph, double ah,  int own, int qi, int na,  int pr);
void amoeba_r(double p[3][2], double y[], int ndim, double ftol, 
		double (*funk)(double[], double m, double hs, double ah,  int zi, double rho, double Omega, double w,
			double hra_hr, double ph, int own, int qi,  int na,  int pr),
		int *nfunk,  double m, double hs, int zi, double rho, double Omega, double w, double hra_hr, double ph, double ah, 
		int own, int qi, int na,  int pr);
double err_lb(double ka_k, double Omega, double rho, double tfp); 
double ucont_na(double choice[],  double m, double hs, int zi, double rho, double  Omega, double  w, double ar, 
		double ph, double ah,  int own,  int qi, int pr);
void amoeba_na(double p[4][3], double y[], int ndim, double ftol, 
		double (*funk)(double[], double m, double hs,  int zi, double rho, double Omega,  double w,
			double hra_hr, double ph, double ah,  int own,   int qi, int pr),
		int *nfunk,  double m, double hs, int zi, double rho, double Omega,  double w, double hra_hr, double ph, double ah,  
		int own, int qi, int pr);
double amotry_na(double p[4][3], double y[], double psum[], int ndim, 
		double (*funk)(double[], double m, double hs,  int zi, double rho, double Omega, double w,
			double hra_hr, double ph, double ah,  int own,  int qi, int pr ),
		int ihi, double fac, double m, double hs,  int zi, double rho, double Omega,  double w, double hra_hr,
		double ph, double ah, int own, int qi, int pr);
void nelder_na(double prime1[], double prime2[],  double ytemp, int ndim, double *hpprime, double *aprime, 
		double *xprime, double *hsrprime, double *ystar, double m, double hs, int zi, double rho, double Omega, 
		double w, double hra_hr, double ph, double ah,  int own, int qi,  int pr);
double value_function_loop_na( double rho, double Omega, double eta, double w, double ar, double ph,  int initial, int region);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splie2(double x1a[], double x2a[], double **ya, int m, int n, int hs_start,  double **y2a);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, int hs_start,  double x1, double x2, double *y);
void amoeba_nr(double p[4][3], double y[], int ndim, double ftol, 
		double (*funk)(double[], double m, double hs, double ah, int zi, double rho, double Omega,  double w,
			double ar, double ph, int own,  int qi, int na, int pr),
		int *nfunk,  double m, double hs, int zi, double rho, double Omega,  double w, double ar, double ph, 
	       double ah, int own, int qi, int na,  int pr);
double amotry_nr(double p[4][3], double y[], double psum[], int ndim, 
		double (*funk)(double[], double m, double hs, double ah,  int zi, double rho, double Omega, double w,
			double hra_hr, double ph, int own,  int qi, int na, int pr ),
		int ihi, double fac, double m, double hs,  int zi, double rho, double Omega,  double w, double hra_hr,
		double ph, double ah, int own, int qi, int na, int pr);
void nelder_nr(double prime1[], double prime2[],  double ytemp, int ndim, double *hpprime, double *haprime, double *aprime, 
		double *xprime, double *hsrprime, double *ystar, double m, double hs, int zi, double rho, double Omega, 
		double w, double ar, double ph, double ah,  int own, int qi, int na, int pr); 
void cond_v(double rho, double  Omega, double  w, double ar, double ph);
