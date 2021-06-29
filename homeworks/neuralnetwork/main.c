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

ann* network;
gsl_vector* xs;
gsl_vector* ys;

double cost_func(gsl_vector* p){
	gsl_vector_memcpy(network->params,p);
	double sum=0;
	for(int i=0;i<xs->size;i++){
		double xi=gsl_vector_get(xs,i);
		double yi=gsl_vector_get(ys,i);
		double fi=response(network,xi);
		sum+=fabs(fi-yi);
	}
	return sum/xs->size;
}

int main(){
	int n=4;
	ann* network=ann_alloc(n,act_func);
	double a=-1,b=1;
	int m=10;
	gsl_vector* vx=gsl_vector_alloc(m);
	gsl_vector* vy=gsl_vector_alloc(m);

	for(int i=0;i<m;i++){
		double x=a+(b-a)*i/(m-1);
		gsl_vector_set(vx,i,x);
		double f=fit_func(x);
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

	for(double z=a;z<=b;z+=1.0/64){
		double y=response(network,z);
		printf("%g %g\n",z,y);
	}
	
	gsl_vector_free(vx);
	gsl_vector_free(vy);
	ann_free(network);

return 0;
}
