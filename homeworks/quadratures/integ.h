#ifndef HAVE_INTEG_H
#define HAVE_INTEG_H
double adapt24(double f(double),double a, double b, 
		double acc, double eps, double f2, double f3,int lim);
double adapt_integ(double f(double),double a, double b, double acc, double eps);
#endif
