#include "stdio.h"
#include "math.h"
#include "complex.h"
int main(void){
	int m=5;
	double n=0.5;
	double complex o=csqrt(-2.0);
	double complex p=cexp(I*M_PI);
	double complex r=cpow(I,M_E);
	double s=cpow(I,I);
	float float_x=1.f/9;
	double x_double=1./9;
	long double x_long_double=1.L/9;
	printf("gammafunction of 5 = %g \n", tgamma(m));
	printf("bessel function of 0.5 = %g \n", j0(n));
	printf("sqrt(-2)=%f%+fi \n", crealf(o), cimagf(o));
	printf("exp(i*pi)=%f%+fi \n", crealf(p), cimagf(p));
	printf("i^e=%f%+fi \n", crealf(r), cimagf(r));
	printf("i^i=%f%+fi \n", crealf(s), cimagf(s));
	printf("1/9 (float) = %.25g \n", float_x);
	printf("1/9 (double) = %.25lg \n", x_double);
	printf("1/9 (long double) = %.25Lg \n", x_long_double);
return 0;
}
