#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
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

double cost_func(gsl_vector *p);

void train(ann* network, gsl_vector* xs, gsl_vector* ys){
	gsl_vector* p=gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);
	qnewt(cost_func,p,1e-3);
	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}


