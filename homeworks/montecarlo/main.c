#include<complex.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>

complex plainmonte(int dim, double f(double* x), double* a, double *b, int N);

int main(){

	double f(double* x){
	return 1./(1.-cos(x[0])*cos(x[1])*cos(x[2]))/M_PI/M_PI/M_PI;
	}

	//Boundaries
	double a[]={0,0,0};
	double b[]={M_PI,M_PI,M_PI};

	//Number of points
	int N=1e6;
	int dim=3;

	complex integral=plainmonte(dim,f,a,b,N);
	printf("Value of integral (plain Montecarlo) =%g +- %g\n",creal(integral), cimag(integral));
	printf("Value of integral =1.39 (approx)\n");

return 0;
}
