#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f (double x, void* params){
	double f = log(x)/sqrt(x);
	return f;
}

int main(){
	gsl_function F;
	F.function=&f;
	F.params=0;
	int limit=999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0,b=1,result,error,acc=1e-6,eps=1e-6;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	printf("result of integral = %g, error = %g \n",result,error);
	printf("integral = -4 (from wolframalpha)\n");
	return 0;
}
