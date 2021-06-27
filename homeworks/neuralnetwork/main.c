#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include"neurons.h"

double act_func(double x){
	return x*exp(-x*x);
}

double fit_func(double x){
	return cos(5*x-1)*exp(-x*x);
}

int main(){
	int n=4;
	ann* network=ann_alloc(n,act_func);
	double a=-1,b=1;
	int nx=20;
	gsl_vector* vx=gsl_vector_alloc(nx);
	gsl_vector* vy=gsl_vector_alloc(nx);

	for(int i=0;i<nx;i++){
		double x=a+(b-a)*i/(nx-1);
		double f=fit_func(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
	}

	for(int i=0;i<network->n;i++){
		gsl_vector_set(network->params,3*i+0,a+(b-a)*i/(network->n-1));
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,1);
	}

	train(network,vx,vy);

	for(int i=0;i<vx->size;i++){
		double x=gsl_vector_get(vx,i);
		double f=gsl_vector_get(vy,i);
		printf("%g %g \n",x,f);
	}
	printf("\n\n");

	for(double z=a;z<b;z+=1./64){
		double y=response(network,z);
		printf("%g %g\n",z,y);
	}

	gsl_vector_free(vx);
	gsl_vector_free(vy);
	ann_free(network);

return 0;
}
