#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
#include<gsl/gsl_interp.h>
#include<assert.h>
#include<math.h>

double linterp(int n, double* x, double* y, double z);
double linterp_integ(int n, double x[],double y[],double z);

int main(){
	int n=10;
	double x[n],y[n];
	for(int i=0;i<n;i++){
		x[i]=0.5+i;
		y[i]=0.3*i+RND;
		printf("%g %g\n",x[i],y[i]);
	}
	printf("\n\n");

	gsl_interp* gsl_linterp = gsl_interp_alloc(gsl_interp_linear,n);
	gsl_interp_init(gsl_linterp,x,y,n);

	double dz=0.1;
	for(double z=x[0];z<=x[n-1];z+=dz){
		double yi=linterp(n,x,y,z);
		double yi_i=linterp_integ(n,x,y,z);

		double yiGSL=gsl_interp_eval(gsl_linterp,x,y,z,NULL);
		double yiGSL_i=gsl_interp_eval_integ(gsl_linterp,x,y,x[0],z,NULL);

		printf("%g %g %g %g %g\n",z,yi,yi_i,yiGSL,yiGSL_i);
	}
	gsl_interp_free(gsl_linterp);
return 0;
}

