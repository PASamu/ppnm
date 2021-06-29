#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>

void driver12(int n, void f(double t, double y[], double dydt[]),
		double a, double y[], double b, double h, double acc, double eps);

void SIR(double T_c, double t,double y[], double dydt[]){

	double N=6*1e6;	 	//Population size
	double T_r=7;		//Recovery time

	dydt[0]=-y[0]*y[1]/N/T_c;
	dydt[1]=y[0]*y[1]/N/T_c-y[1]/T_r;
	dydt[2]=y[1]/T_r;
}

static double contact_time;

void SIR_with_Tc(double t, double y[],double dydt[]){
	SIR(contact_time,t,y,dydt);
}

int main(){
	int a=0;
	int b=100;
	
	double h=0.1,acc=1e-3,eps=1e-3;
	double T_c[]={2,4,7}; 	//Time between contacts
	int n = sizeof(T_c)/sizeof(T_c[0]);
	double y[n];

	FILE* sus = fopen ("sus.txt", "w");	//File containing data on susceptible
	FILE* inf = fopen ("inf.txt", "w");	//File containing data on infectious
	FILE* rem = fopen ("rem.txt", "w");	//File containing data on removed

	for(int t=1;t<b;t++){
		//Initial values
		y[0]=6*1e6;
		y[1]=100;
		y[2]=0;

		fprintf (sus,"%i ", t);
		fprintf (inf,"%i ", t);
		fprintf (rem,"%i ", t);

		for(int i=0;i<n;i++){
			
			contact_time = T_c[i];
			
			driver12(n,SIR_with_Tc,(double) a,y,(double) t,h,acc,eps);
			
			fprintf (sus, "%15g ", y[0]);
			fprintf (inf, "%15g ", y[1]);
			fprintf (rem, "%15g ", y[2]);
		}
		fprintf (sus, "\n");
		fprintf (inf, "\n");
		fprintf (rem, "\n");
	}
	fclose (sus);
	fclose (inf);
	fclose (rem);
	
	return 0;
}
