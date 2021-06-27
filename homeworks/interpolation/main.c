#include<stdio.h>
#include<stdlib.h>
#define RND (double)rand()/RAND_MAX
#include<gsl/gsl_interp.h>
#include<math.h>

double linterp(int n, double* x, double* y, double z);
double linterp_integ(int n, double x[],double y[],double z);

int main(){
	int n=10;
	double x[n],y[n];
	for(int i=0;i<n;i++){
		x[i]=0.5+i;
		y[i]=0.2*i+RND;
		printf("%g %g\n",x[i],y[i]);
	}



	

return 0;
}

