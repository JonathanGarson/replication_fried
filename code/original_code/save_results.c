/*
Program: save_results.c
Purpose: saves model results as dat files to be read into matlab
*/

#include "header.h"

extern double *v0, *v1, *hppol, *hapol, *apol, *xpol, *hsrpol,   *f1, *f2, *mgrid, *cont_e0, *cont_e1; 
extern char *exper, *folder; 
extern int  nm, nhs,  nz, nq, ccexp;

void save_results(int region, char *exper)
{
	int n,  i;
	char filename[100];
	FILE *rs; 

	n = nhs*nm*nz*nq;
	

	if ((strcmp(exper, "cc") ==0) || (strcmp(exper, "ccna") ==0)|| (strcmp(exper, "ccnr") ==0)) {
		//v0
		sprintf(filename, "./%s/v0_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", v0[i]);
		fclose(rs);

		//v1
		sprintf(filename, "./%s/v1_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", v1[i]);
		fclose(rs);

		//hppol
		sprintf(filename, "./%s/hppol_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", hppol[i]);
		fclose(rs);

		//hapol
		sprintf(filename, "./%s/hapol_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL) printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", hapol[i]);
		fclose(rs);

		//apol
		sprintf(filename,  "./%s/apol_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", apol[i]);
		fclose(rs);

		//distribution
		sprintf(filename,  "./%s/f1_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", f1[i]);
		fclose(rs);

		//xpol
		sprintf(filename,  "./%s/xpol_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", xpol[i]);
		fclose(rs);

		//hsrpol
		sprintf(filename, "./%s/hsrpol_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", hsrpol[i]);
		fclose(rs);
		
		//cont e0
		sprintf(filename, "./%s/conte0_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", cont_e0[i]);
		fclose(rs);

		//cont e1
		sprintf(filename, "./%s/conte1_%s%i_%i.dat", folder, exper, region, ccexp);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", cont_e1[i]);
		fclose(rs);

	} else{

		//v0
		sprintf(filename, "./%s/v0_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", v0[i]);
		fclose(rs);

		//v1
		sprintf(filename, "./%s/v1_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", v1[i]);
		fclose(rs);

		//hppol
		sprintf(filename, "./%s/hppol_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", hppol[i]);
		fclose(rs);

		//hapol
		sprintf(filename, "./%s/hapol_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL) printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", hapol[i]);
		fclose(rs);

		//apol
		sprintf(filename,  "./%s/apol_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", apol[i]);
		fclose(rs);

		//distribution
		sprintf(filename,  "./%s/f1_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", f1[i]);
		fclose(rs);

		//xpol
		sprintf(filename,  "./%s/xpol_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", xpol[i]);
		fclose(rs);

		//hsrpol
		sprintf(filename, "./%s/hsrpol_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", hsrpol[i]);
		fclose(rs);
		
		//cont e0
		sprintf(filename, "./%s/conte0_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", cont_e0[i]);
		fclose(rs);
	
		//cont e1
		sprintf(filename, "./%s/conte1_%s%i.dat", folder, exper, region);
		rs = fopen(filename, "w");
		if (rs == NULL)	printf("The file did not open\n");
		for (i = 0; i < n; i++) fprintf(rs, "%12.8f\n ", cont_e1[i]);
		fclose(rs);
	}
}
