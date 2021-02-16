#include "stdio.h"
#include "math.h"
#include "complex.h"
int main(void){
	int n=5;
	double m=0.5;
	double complex o=csqrt(-2.0);
	double p=cexp(I*M_PI);
	double complex q=cexp(I);
	double complex r=cpow(I,M_E);
	double complex s=cpow(I,I);
	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	printf("gammafunction of 5 = %g \n",tgamma(n));
	printf("bessel function of 0.5 = %g \n",j0(m));
	printf("sqrt(-2)=%f%+fi \n", crealf(o), cimagf(o));
	printf("exp(i*pi)=%f%+fi \n", crealf(p), cimagf(p));
	printf("exp(i)=%f%+fi \n", crealf(q), cimagf(q));
	printf("i^e=%f%+fi \n", crealf(r),cimagf(r));
	printf("i^i=%f%+fi \n", crealf(s),cimagf(s));
	printf("1/9 (float) = %.25g \n", x_float);
	printf("1/9 (double) = %.25lg \n", x_double);
	printf("1/9 (long double) = %.25Lg \n", x_long_double);
return 0;
}
