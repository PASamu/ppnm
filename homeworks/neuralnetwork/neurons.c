#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include"neurons.h"

void num_grad(double f(gsl_vector*),gsl_vector* x,gsl_vector* grad);

int qnewt(double f(gsl_vector*),gsl_vector*x,double eps);

ann* ann_alloc(int n, double(*f)(double)){
	ann* network = malloc(sizeof(ann));
	network->n=n;
	network->f=f;
	network->params=gsl_vector_alloc(3*n);
	return network;
}

void ann_free(ann* network){
	gsl_vector_free(network->params);
	free(network);
}

double response(ann* network,double x){
	double s=0;
	for(int i=0;i<network->n;i++){
		double a=gsl_vector_get(network->params,3*i+0);
	 	double b=gsl_vector_get(network->params,3*i+1);
		double w=gsl_vector_get(network->params,3*i+2);
		s+=network->f((x-a)/b)*w;
	}
	return s;
}

void train(ann* network, gsl_vector* xs, gsl_vector* ys){
	double cost_func(gsl_vector* p){
		gsl_vector_memcpy(network->params,p);
		double sum=0;
		for(int i=0;i<xs->size;i++){
			double xi=gsl_vector_get(xs,i);
			double yi=gsl_vector_get(ys,i);
			double fi=response(network,xi);
			sum+=(fi-yi)*(fi-yi);
		}
		return sum/xs->size;
	}
	gsl_vector* p=gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);

	qnewt(cost_func,p,1e-3);
	gsl_vector_memcpy(network->params,p);

	gsl_vector_free(p);
}


