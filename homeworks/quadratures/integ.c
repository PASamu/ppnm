#include<math.h>
#include<assert.h>
#include<stdio.h>

double adapt24(double f(double x), double a, double b, double acc,
		double eps, double f2, double f3, int lim){
	double f1=f(a+(b-a)/6),f4=f(a+5*(b-a)/6);
	double Q=(2*f1+f2+f3+2*f4)/6*(b-a), q=(f1+f2+f3+f4)/4*(b-a);
	double tol=acc+eps*fabs(Q),err=fabs(Q-q);
	if(lim==0){
		fprintf(stderr,"adapt24: recursion limit reached\n");
		return Q;
	}
	if(err<tol) return Q;
	else{
		double Q1=adapt24(f,a,(a+b)/2,acc/sqrt(2),eps,f1,f2,lim-1);
		double Q2=adapt24(f,(a+b)/2,b,acc/sqrt(2),eps,f3,f4,lim-1);
		return Q1+Q2;}
}

double adapt_integ(double f(double x),double a, double b,
		double acc, double eps){
	double f2=f(a+2*(b-a)/6),f3=f(a+4*(b-a)/6);
	int lim=999;
	return adapt24(f,a,b,acc,eps,f2,f3,lim);
}
	
