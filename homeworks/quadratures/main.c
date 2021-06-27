#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integ.h"

int main(){
	double a=0,b=1;
	double acc=1e-3,eps=1e-3;

	double f(double x){
		return sqrt(x);
	}
	double I=adapt_integ(f,a,b,acc,eps);

	printf("Integral of sqrt(x) from 0 to 1 (adaptive)=%g\n",I);
	printf("Integral of sqrt(x) from 0 to 1 (analytical)=%g\n",2./3);
	
	double g(double x){
		return 4*sqrt(1-x*x);
	}
	double I2=adapt_integ(g,a,b,acc,eps);
	printf("Integral of 4sqrt(1-x^2) from 0 to 1 (adaptive)=%g\n",I2);
	printf("Integral of 4sqrt(1-x^2) from 0 to 1 (analytical)=pi\n");
}
